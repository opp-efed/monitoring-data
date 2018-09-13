import arcpy
import os

vpus_nhd = {'01': 'NE', '02': 'MA', '03N': 'SA', '03S': 'SA', '03W': 'SA', '04': 'GL', '05': 'MS',
            '06': 'MS', '07': 'MS', '08': 'MS', '09': 'SR', '10L': 'MS', '10U': 'MS', '11': 'MS',
            '12': 'TX', '13': 'RG', '14': 'CO', '15': 'CO', '16': 'GB', '17': 'PN', '18': 'CA'}


def overlay(paths, id_fields):
    temp_file = "in_memory/intersect"
    layers = []
    for i, path in enumerate(paths):
        layer = 'lyr{}'.format(i)
        arcpy.MakeFeatureLayer_management(path, layer)
        layers.append(layer)
    arcpy.Intersect_analysis(layers, temp_file)
    try:
        if arcpy.GetCount_management(temp_file) > 0:
            results = {tuple(row) for row in arcpy.da.SearchCursor(temp_file, id_fields)}
        else:
            results = set()
    except Exception as e:
        print(e)
        results = set()
    try:
        for bogey in layers + [temp_file]:
            arcpy.Delete_management(bogey)
    except Exception as e:
        print(e)
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
        local_sites = overlay([points_file, catchment_path.format(vpus_nhd.get(region), region)],
                              [site_id_field, catchment_field])
        sites |= {(site_id, comid, region) for site_id, comid in local_sites}

    write_output(sites, sites_file)


main()
