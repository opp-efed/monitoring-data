import arcpy
import os


def initialize_temp():
    from paths import scratch_workspace
    temp_file = os.path.join(scratch_workspace, "intersect.shp")
    if os.path.exists(temp_file):
        arcpy.Delete_management(temp_file)
    return temp_file


def overlay(paths, id_fields):
    temp_file = initialize_temp()
    layers = []
    for i, path in enumerate(paths):
        layer = 'lyr{}'.format(i)
        arcpy.MakeFeatureLayer_management(path, layer)
        layers.append(layer)
    arcpy.Intersect_analysis(layers, temp_file)
    if arcpy.GetCount_management(temp_file) > 0:
        results = {tuple(row) for row in arcpy.da.SearchCursor(temp_file, id_fields)}
    else:
        results = set()
    for bogey in layers + [temp_file]:
        arcpy.Delete_management(bogey)
    return results


def write_output(sites, outfile):
    with open(outfile, 'w') as f:
        f.write("SITE_ID,COMID,REGION\n")
        for site_id, comid, state in sites:
            f.write("{},{},{}\n".format(site_id, comid, state))


def main():
    from paths import catchment_path, points_file, sites_file, site_id_field, regions_file, region_id_field

    catchment_field = "FEATUREID"

    sites = set()

    # Get the regions which contain points
    active_regions = \
        sorted({r[1] for r in overlay([points_file, regions_file], [site_id_field, region_id_field])})

    for region in active_regions:
        print(region)
        local_sites = overlay([points_file, catchment_path.format(region)], [site_id_field, catchment_field])
        sites |= {(site_id, comid, region) for site_id, comid in local_sites}

    write_output(sites, sites_file)


main()
