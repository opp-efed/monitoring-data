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


def get_grid_props(grid_layer, catchment_layer, upstream_reaches, grid_id_field):
    temp_file = os.path.join("..", "bin", "temp", "clip.shp")
    if os.path.exists(temp_file):
        arcpy.Delete_management(temp_file)
    sel_string = "\"FEATUREID\" = " + " OR \"FEATUREID\" = ".join(map(str, upstream_reaches))
    arcpy.SelectLayerByAttribute_management(catchment_layer, "NEW_SELECTION", sel_string)
    clip = arcpy.Clip_analysis(grid_layer, catchment_layer, temp_file)
    grid_props = \
        [(station_id, area) for station_id, area in arcpy.da.SearchCursor(clip, [grid_id_field, "SHAPE@AREA"])]
    arcpy.Delete_management(temp_file)
    return grid_props


def write_output(output, outfile):
    with open(outfile, 'w') as f:
        for row in output:
            f.write(",".join(map(str, row)) + "\n")


def main():
    from paths import grid_file, watershed_file, catchment_path, grid_props, grid_id_field

    # Loop through sites, sorted by region
    output = [["Site_ID", "Comid", "Dayshed", "Grid_ID", "Grid_Area"]]
    grid_layer = arcpy.MakeFeatureLayer_management(grid_file, "grid_layer")
    for region, sites in get_watersheds(watershed_file):
        print(region)
        catchment_layer = arcpy.MakeFeatureLayer_management(catchment_path.format(region), "catch_layer")
        for site in sites:
            site_id, comid, time = site[:3]
            upstream_reaches = site[3:]
            print(site_id, time, len(upstream_reaches))
            for grid_id, area in get_grid_props(grid_layer, catchment_layer, upstream_reaches, grid_id_field):
                output.append([site_id, comid, time, grid_id, area])
        arcpy.Delete_management(catchment_layer)

    write_output(output, grid_props)


main()
