import arcpy
import os
from collections import defaultdict


def get_watersheds(watersheds_table):
    out_dict = defaultdict(list)
    with open(watersheds_table) as f:
        for line in f:
            data = line.strip().split(',')
            site_id, comid, region, time = data[:4]
            out_dict[region].append([site_id, comid, time] + data[4:])
    return sorted(out_dict.items())


def get_grid_props(grid_layer, catchment_layer, upstream_reaches):
    temp_file = os.path.join("..", "bin", "temp", "clip.shp")
    if os.path.exists(temp_file):
        arcpy.Delete_management(temp_file)
    sel_string = "\"FEATUREID\" = " + " OR \"FEATUREID\" = ".join(map(str, upstream_reaches))
    arcpy.SelectLayerByAttribute_management(catchment_layer, "NEW_SELECTION", sel_string)
    clip = arcpy.Clip_analysis(grid_layer, catchment_layer, temp_file)
    return [(station_id, area) for station_id, area in arcpy.da.SearchCursor(clip, ["stationID", "SHAPE@AREA"])]


def write_output(output, outfile):
    with open(outfile, 'w') as f:
        for row in output:
            f.write(",".join(map(str, row)) + "\n")


def main():
    # Set paths
    grid_file = r"C:\Users\Trip Hook\Documents\BiasFactors\Data\weather_stations_highres_thiessens_US_alb\weather_stations_highres_thiessens_US_alb.shp"
    nhd_dir = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"
    watersheds_table = os.path.join("..", "bin", "Tables", "amp_watersheds_061418_2.txt")
    outfile = os.path.join("..", "bin", "Tables", "amp_061418_grid_props.txt")

    # Loop through sites, sorted by region
    output = [["Site_ID", "Comid", "Dayshed", "Grid_ID", "Grid_Area"]]
    grid_layer = arcpy.MakeFeatureLayer_management(grid_file, "grid_layer")
    for region, sites in get_watersheds(watersheds_table):
        print(region)
        catchment_shapefile = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusCatchment", "Catchment.shp")
        catchment_layer = arcpy.MakeFeatureLayer_management(catchment_shapefile, "catch_layer")
        for site in sites:
            site_id, comid, time = site[:3]
            if int(time) < 8:
                print(site_id, time)
                upstream_reaches = site[3:]
                for grid_id, area in get_grid_props(grid_layer, catchment_layer, upstream_reaches):
                    output.append([site_id, comid, time, grid_id, area])
            else:
                print("Skipping {}, {} for now".format(site_id, time))
        arcpy.Delete_management(catchment_layer)

    write_output(output, outfile)


main()