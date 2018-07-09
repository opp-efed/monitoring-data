import os
import pandas as pd
import numpy as np

from classes import MetfileMatrix


def main():
    sam_dir = r"C:\Users\Trip Hook\Documents\opp-efed\sam"
    metfile_path = os.path.join(sam_dir, "bin", "Tables", "WeatherFiles", "metfile")
    grid_props_table = os.path.join("..", "bin", "Tables", "amp_061418_grid_props.txt")
    outfile = os.path.join("..", "bin", "Tables", "amp_061418_precip.txt")

    met = MetfileMatrix(metfile_path)

    grid_props = pd.read_csv(grid_props_table)
    all_sites = None
    for site_id in np.unique(grid_props.Site_ID):
        print(site_id)
        site_rows = grid_props[grid_props.Site_ID == site_id]
        precip_series = np.zeros((2, met.n_dates))
        site_area = 0
        for dayshed in np.unique(site_rows.Dayshed):
            dayshed_rows = site_rows[site_rows.Dayshed == dayshed]
            for _, row in dayshed_rows.iterrows():
                site_area += row.Grid_Area
                precip = met.fetch_station(str(row.Grid_ID))[0] * row.Grid_Area
                precip_series[0] += precip
                if dayshed > 0:
                    precip_series[1][int(dayshed):] += precip[:-int(dayshed)]
                else:
                    precip_series[1] += precip
        labels = ["{}_uniform".format(site_id), "{}_tot".format(site_id)]
        precip_series = pd.DataFrame(index=labels, data=precip_series)
        if all_sites is None:
            all_sites = precip_series
        else:
            all_sites = pd.concat([all_sites, precip_series])
        print(site_area)
    all_sites.to_csv(outfile)


main()