import os
import pandas as pd
import numpy as np
from urllib.request import urlopen

from warp_lib import allocate, read_dbf, overlaps, Raster
import os
import numpy as np
import pandas as pd


class Navigator(object):
    def __init__(self, region_id, upstream_path):
        self.file = upstream_path.format(region_id)
        self.paths, self.times, self.map, self.alias_to_reach, self.reach_to_alias = self.load()

    def load(self):
        assert os.path.isfile(self.file), "Upstream file {} not found".format(self.file)
        data = np.load(self.file, mmap_mode='r')
        conversion_array = data['alias_index']
        reverse_conversion = dict(zip(conversion_array, np.arange(conversion_array.size)))
        return data['paths'], data['time'], data['path_map'], conversion_array, reverse_conversion

    def upstream_watershed(self, reach_id, mode='reach', return_times=True, verbose=False):

        def unpack(array):
            first_row = [array[start_row][start_col:]]
            remaining_rows = list(array[start_row + 1:end_row])
            return np.concatenate(first_row + remaining_rows)

        # Look up reach ID and fetch address from upstream object
        reach = reach_id if mode == 'alias' else self.reach_to_alias.get(reach_id)

        try:
            start_row, end_row, col = map(int, self.map[reach])
            start_col = list(self.paths[start_row]).index(reach)
        except TypeError:
            if verbose:
                print("Reach {} not found in region".format(reach))
            return [] if not return_times else ([], [])
        except ValueError:
            if verbose:
                print("{} not in upstream lookup".format(reach))
            return [] if not return_times else ([], [])
        except IndexError:
            if verbose:
                print("Reach {} not found in region".format(reach))
            return [] if not return_times else ([], [])
        else:
            # Fetch upstream reaches and times
            aliases = unpack(self.paths)
            reaches = aliases if mode == 'alias' else np.int32(self.alias_to_reach[aliases])
        if not return_times:
            return reaches
        else:
            times = unpack(self.times)
            adjusted_times = np.int32(times - self.times[start_row][start_col])
            return reaches, adjusted_times


def iterregions(upstream_path):
    regions = [f.split("_")[1].rstrip(".npz") for f in os.listdir(upstream_path)]
    regions = ((r, os.path.join(upstream_path, "upstream_{}.npz".format(r))) for r in regions)
    for region, region_file in sorted(regions):
        yield region, region_file


def histogram(allocation, watershed, progress=True):
    upstream = np.zeros(allocation.shape)
    for reach in range(1, watershed.path_map.shape[0]):
        if progress and not reach % 10000:
            print("\t\tAccumulating reach {} of {}".format(reach, watershed.path_map.shape[0]))
        upstream_reaches = watershed.all_upstream(reach, mode='alias')
        if upstream_reaches.any():
            upstream[reach] = np.take(allocation, upstream_reaches, axis=0).sum(axis=0)
    return upstream


def accumulate_histogram(table, nav):
    alias_to_reach = pd.DataFrame(nav.alias_to_reach, columns=["COMID"], dtype=np.int32)
    table = pd.merge(alias_to_reach, table, how='left', on="COMID").set_index("COMID")

    table['cn_nml'] = table.cn * table.LocalArea
    data = table[['LocalArea', 'cn_nml']].as_matrix()

    header = ['AreaSum', 'cn_Ws']
    new_table = np.zeros((table.shape[0], len(header)))
    for i, (index, row) in enumerate(table.iterrows()):
        if not i % 1000:
            print(i)

        upstream_reaches = nav.upstream_watershed(nav.reach_to_alias[index], mode='alias', return_times=False)
        upstream_array = np.nan_to_num(np.take(data, upstream_reaches, axis=0))

        new_table[i] = upstream_array.sum(axis=0)
        new_table[i, 1] /= new_table[i, 0]  # average by area

    old_table = table[[field for field in table.columns if not field.endswith("_nml")]]
    new_table = pd.DataFrame(data=new_table, columns=header, index=table.index)
    full_table = pd.merge(old_table, new_table, how='inner', left_index=True, right_index=True)

    return full_table


def accumulate(warp_table, nhd_dir, region, upstream_dir):
    # Load area table and join to warp table
    area_table = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusAttributes", "PlusFlowlineVAA.dbf")
    area_table = pd.DataFrame(data=read_dbf(area_table, ["ComID", "AreaSqKM", "TotDASqKM"]),
                              columns=["COMID", "LocalArea", "TotalArea"])

    warp_table = warp_table.merge(area_table, on="COMID", how='inner')

    if region == '01':
        print(warp_table)

    # Create navigator object
    nav = Navigator(region, os.path.join(upstream_dir, "upstream_{}.npz".format(region)))

    # Run accumulation
    accumulated = accumulate_histogram(warp_table, nav)

    return accumulated


def zonal_mean(value_raster, watershed_raster, name, gridcode_table):
    allocation = allocate(value_raster, watershed_raster, 100000)
    allocation = pd.DataFrame(data=allocation.T, columns=["Gridcode", "Value", "Area"])
    gridcode_to_feature = pd.Series(dict(read_dbf(gridcode_table, ["GRIDCODE", "FEATUREID"])), name="COMID").to_frame()
    allocation = allocation.groupby("Gridcode").apply(lambda g: (g["Value"] * g["Area"]).sum() / g["Area"].sum())
    allocation.name = name
    allocation = pd.merge(allocation.to_frame(), gridcode_to_feature, left_index=True, right_index=True)

    return allocation


def write_to_file(region, out_dir, all_data, gridcode_table):
    outfile = os.path.join(out_dir, "region_{}.csv".format(region))
    all_data.reset_index().to_csv(outfile, index=False)


def main():
    global overwrite

    # Initial paths
    nhd_dir = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"
    upstream_dir = r"S:\bin\Preprocessed\Upstream"
    out_dir = r"C:\Users\Trip Hook\Documents\BiasFactors\CurveNumber"
    value_raster = \
        r"C:\Users\Trip Hook\Documents\WatershedTools_new\SpecialApplications\CurveNumberMap\Results\national_30.tif"

    overwrite = True

    # Loop through NHD regions
    value_raster = Raster(value_raster, no_data=0)

    for region, states in sorted(overlaps.items()):
        # Set regional paths
        gridcode_table = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusCatchment", "featureidgridcode.dbf")

        catchment_raster = \
            Raster(os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusCatchment", "cat"), no_data=255)

        result = zonal_mean(value_raster, catchment_raster, "cn", gridcode_table)

        result = accumulate(result, nhd_dir, region, upstream_dir)

        write_to_file(region, out_dir, result, gridcode_table)


main()

