import arcpy
import os


def overlay(points_file, catchment_file):
    temp_file = r"..\Data\intersect.shp"
    if os.path.exists(temp_file):
        arcpy.Delete_management(temp_file)
    points = arcpy.MakeFeatureLayer_management(points_file, 'lyr1')
    catchments = arcpy.MakeFeatureLayer_management(catchment_file, 'lyr2')
    arcpy.Intersect_analysis([points, catchments], temp_file)
    if arcpy.GetCount_management(temp_file) > 0:
        results = {field1: comid for field1, comid in arcpy.da.SearchCursor(temp_file, ['WSDA_Site', "FEATUREID"])}
    else:
        results = {}
    for bogey in ("lyr1", "lyr2", temp_file):
        arcpy.Delete_management(bogey)
    return results


def write_output(site_dict, outfile):
    with open(outfile, 'w') as f:
        f.write("SITE_ID,COMID\n")
        for site_id, comid in site_dict.items():
            f.write("{},{}\n".format(site_id, comid))


def main():
    from utilities import nhd_states
    from paths import nhd_dir, points_file, sites_file
    site_dict = {}

    for region in nhd_states:
        catchment_file = os.path.join(nhd_dir, "NHDPlus{}".format(region), "NHDPlusCatchment", "Catchment.shp")
        results = overlay(points_file, catchment_file)
        site_dict.update(results)

    write_output(site_dict, sites_file)


main()
