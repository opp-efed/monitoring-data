import pandas as pd
import numpy as np

from classes import MetfileMatrix


def pull_precip(grid_props, met):
    site_ids = np.unique(grid_props.Site_ID)
    columns = [site_id + "_" + field for site_id in site_ids for field in ("uniform", "tot")]
    all_sites = pd.DataFrame(data=np.zeros((met.n_dates, len(columns))), columns=columns, index=met.dates)
    for (site_id, dayshed), group in grid_props.groupby(['Site_ID', 'Dayshed']):
        for _, row in group.iterrows():
            precip = met.fetch_station(str(row.Grid_ID))[0] * row.Grid_Area
            all_sites["{}_uniform".format(site_id)] += precip
            if dayshed > 0:
                all_sites["{}_tot".format(site_id)][int(dayshed):] += precip[:-int(dayshed)]
            else:
                all_sites["{}_tot".format(site_id)] += precip
    return all_sites


def main():
    from paths import metfile_path, grid_props, outfile

    met = MetfileMatrix(metfile_path)

    grid_props = pd.read_csv(grid_props)

    precip_series = pull_precip(grid_props, met)

    precip_series.to_csv(outfile)


main()
