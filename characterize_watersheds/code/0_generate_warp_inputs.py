import os
import pandas as pd
import numpy as np
import gdal
import math


class WarpInputs(object):
    def __init__(self, region, states, use_table, county_raster, use_raster, nhd_dir, allocation_dir, overwrite=False):
        self.region = region
        self.states = states
        self.county_raster = county_raster
        self.use_raster = use_raster
        self.use_table = use_table
        self.allocation_dir = allocation_dir
        self.overwrite = overwrite

        # Set paths
        nhd_dir = os.path.join(nhd_dir, "NHDPlus{}".format(region))

        # Initialize rasters
        catchment_raster = Raster(os.path.join(nhd_dir, "NHDPlusCatchment", "cat"), no_data=255)
        gridcode_table = os.path.join(nhd_dir, "NHDPlusCatchment", "featureidgridcode.dbf")

        # Perform regional allocation (area of use in each county/catchment)
        regional_allocation = self.allocate_use(county_raster, use_raster, catchment_raster)

        # Calculate use density (use per county / pixels per county) and merge with allocation
        density = self.get_use_density(regional_allocation)

        # Multiply area by density to get load
        data = regional_allocation.join(density, on="FIPS")
        data['Load'] = data['Area'].astype(np.float32) * data['Density'].astype(np.float32)
        data.Gridcode = data.Gridcode.astype(np.int64)
        data = data.groupby("Gridcode")["Load"].sum()

    def get_use_density(self, regional_allocation):
        """Create a table of use density in kg/m2 for each FIPS code and year"""
        all_years = None
        for year in sorted(np.unique(self.use_table.YEAR)):
            annual_use = self.use_table[self.use_table.YEAR == year]
            county_use = annual_use.groupby("FIPS")["UseKG"].sum()  # FIPS, KG
            county_area = regional_allocation.groupby("FIPS")["Area"].sum()  # FIPS, Area
            county_stats = pd.concat([county_use, county_area], axis=1, join='inner')
            density = pd.DataFrame(data=county_stats.UseKG.astype(np.float32) / county_stats.Area.astype(np.float32),
                                   index=county_stats.index, columns=["Density"])
            density['Year'] = year
            all_years = density if all_years is None else all_years.append(density)
        return all_years.groupby(all_years.index)["Density"].mean().to_frame()

    def allocate_use(self, county_raster_dir, use_raster, catchment_raster):
        """Allocate use to counties"""
        regional_allocation = None

        # Iterate through all states in region
        for state in self.states:

            # Initialize the statewide FIPS raster
            state_raster = Raster(os.path.join(county_raster_dir, "{}_cty_30m".format(state.lower())))

            # Allocate use to catchment/county rasters
            gridcodes, fips, _, area = use_raster.allocate([catchment_raster, state_raster], 100000)
            fips = np.array([str(f).zfill(5) for f in fips])

            # Append to a complete table
            allocation = pd.DataFrame(data=np.array([gridcodes, fips, area]).T, columns=["Gridcode", "FIPS", "Area"])
            for field in ("Gridcode", "Area"):
                allocation[field] = allocation[field].map(np.int32)
            regional_allocation = allocation if regional_allocation is None else regional_allocation.append(allocation)

        return regional_allocation


class Raster(object):
    def __init__(self, path, no_data=None):
        self.path = path
        self.no_data = no_data
        self.obj = gdal.Open(path)
        gt = self.obj.GetGeoTransform()
        self.cell_size = int(gt[1])
        self.tl = (gt[0], gt[3])  # top left
        self.size = (self.obj.RasterXSize * self.cell_size, self.obj.RasterYSize * self.cell_size)
        self.shape = Envelope(self.tl[0], self.tl[0] + self.size[0], self.tl[1] - self.size[1], self.tl[1])
        self.max_val = int(self.obj.GetRasterBand(1).GetStatistics(True, True)[1])
        if self.no_data:
            self.max_val = max((self.max_val, self.no_data))
        self.precision = 10 ** int(math.ceil(math.log10(self.max_val)))
        self.envelope = None

        self._array = np.array([])

    def array(self, envelope, zero_min=True, datatype=np.int64):
        if not self._array.size or (envelope != self.envelope):
            offset_x = (envelope.left - self.shape.left)
            offset_y = (self.shape.top - envelope.top)
            x_max = (envelope.right - envelope.left)
            y_max = (envelope.top - envelope.bottom)
            bounds = map(lambda x: int(x / self.cell_size), (offset_x, offset_y, x_max, y_max))
            self._array = datatype(self.obj.ReadAsArray(*bounds))
            if zero_min:
                self._array[self._array < 0] = 0
            if self.no_data:
                self._array[self._array == self.no_data] = 0
            self._array = np.int64(self._array)
            self.envelope = envelope
        return self._array

    def allocate(self, zone_rasters, tile='max'):
        """ Allocates raster classes to a set of overlapping zones """

        # Function uses a list of zone rasters. If a single zone is provided, stick it in a list
        zone_rasters = [zone_rasters] if not type(zone_rasters) in (list, set) else zone_rasters

        # Overlap rasters and create envelope covering common areas
        overlap_area = self.shape
        for raster in zone_rasters:
            overlap_area = raster.shape.overlap(overlap_area)
            assert overlap_area, "Zone and allocation rasters do not overlap"

        # Divide the overlap area into tiles to aid in processing
        tiles = overlap_area.make_tiles(tile)

        # Sort all rasters, zone and allocation, by precision, lowest to highest (the magnitude of the numbers in the array)
        all_rasters = sorted(zone_rasters, key=lambda x: x.precision, reverse=True) + [allocation_raster]

        # Iterate through tiles
        finished = None
        for tile in tiles:

            # Extract tile from the allocation raster. This is the starting point
            combined_array = self.array(tile)

            # Loop through the zone rasters,
            # multiply by precision factor to add zeros such that all overlapping data is preserved
            for i, raster in enumerate(zone_rasters):
                # Multiply the zone raster by an adjustment factor based on precision
                raster.adjust = np.prod([r.precision for r in all_rasters[i + 1:]])

                # Add the adjusted zone raster to the combined array.
                combined_array += (raster.array(tile) * raster.adjust)

            self.adjust = 1

            # Break down combined array to get zones and classes
            values, counts = np.unique(combined_array.flat, return_counts=True)  # Cells in each zone and class

            zones = np.zeros((len(all_rasters), counts.size), dtype=np.int64)
            for i, r in enumerate(all_rasters):
                zones[i] = np.int64(values / r.adjust)
                values -= (zones[i] * r.adjust)

            counts *= (self.cell_size ** 2)  # Convert cell count to area

            final = np.vstack((zones, counts))

            # Filter out combinations with no cells
            final = final[:, (np.prod(zones, axis=0) > 0)]

            # Append to running
            finished = final if finished is not None else np.hstack((finished, final))

        return finished


class Envelope(object):
    """ Object representing a simple bounding rectangle, used primarily to measure raster overlap """

    def __init__(self, left, right, bottom, top):
        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top

    # Returns the rectangle corresponding to the overlap with another Envelope object
    def overlap(self, r2):
        def range_overlap(a_min, a_max, b_min, b_max):
            return (a_min <= b_max) and (b_min <= a_max)

        if not all((range_overlap(self.left, self.right, r2.left, r2.right),
                    range_overlap(self.bottom, self.top, r2.bottom, r2.top))):
            return None
        else:
            left, right = sorted([self.left, self.right, r2.left, r2.right])[1:3]
            bottom, top = sorted([self.bottom, self.top, r2.bottom, r2.top])[1:3]
        return Envelope(left, right, bottom, top)

    @property
    def area(self):
        return abs(self.top - self.bottom) * abs(self.right - self.left)

    def make_tiles(self, tile_size):
        if tile_size == 'max':
            return [self]
        else:
            h = list(range(int(self.left), int(self.right), tile_size)) + [self.right]
            v = list(range(int(self.bottom), int(self.top), tile_size)) + [self.top]
            return [Envelope(h[i], h[i + 1], v[j], v[j + 1]) for i in range(len(h) - 1) for j in range(len(v) - 1)]

    def __repr__(self):
        return "Rectangle(left: {}, right: {}, top: {}, bottom: {}".format(self.left, self.right, self.top, self.bottom)

    def __eq__(self, other):
        return (self.left, self.right, self.bottom, self.top) == (other.left, other.right, other.bottom, other.top)


def read_use(use_dir, level):
    all_use = None
    tables = filter(lambda x: x.split("_")[0] == level, os.listdir(use_dir))
    for f in tables:
        print(f)
        table = pd.read_csv(os.path.join(use_dir, f), sep="\t")
        table["FIPS"] = table.pop("STATE_FIPS_CODE").map(lambda x: str(int(x)).zfill(2)) + \
                        table.pop("COUNTY_FIPS_CODE").map(lambda x: str(int(x)).zfill(3))
        table["UseKG"] = table.pop("KG").astype(np.float32)
        all_use = table if all_use is None else all_use.append(table)
    return all_use


def main():
    from utilities import compounds, nhd_states

    # Initial paths
    nhd_dir = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"
    grid_dir = r"C:\Users\Trip Hook\Documents\BiasFactors\Projected"
    county_raster_dir = r"C:\Users\Trip Hook\Documents\NationalData\County_Grids"
    use_dir = os.path.join("..", "bin", "use_tables")
    pmayjun_grid = os.path.join(grid_dir, "pmj_prism_2")
    rfactor_grid = os.path.join(grid_dir, "rfactor_30m")
    satof48_grid = os.path.join(grid_dir, "satof48_30m")
    srl25ag_grid = os.path.join(grid_dir, "srl25ag_30m")
    use_grid = os.path.join("..", "Projected", "atr_cdl")

    level = 'high'  # or 'low'

    overwrite = True

    # Pull use data from online or local storage into pandas tables
    pesticide_use = read_use(use_dir, level)

    # Initialize use raster
    use_raster = Raster(use_grid, no_data=255)

    # Loop through NHD regions
    for region, states in sorted(nhd_states.items()):
        WarpInputs(region, states, pesticide_use, county_raster_dir, use_raster, nhd_dir)


main()