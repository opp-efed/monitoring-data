import numpy as np
import pandas as pd

from classes import Navigator


def main():
    from paths import sites_file, nav_path, watershed_file

    # Read sites file
    sites = pd.read_csv(sites_file)
    sites.columns = map(str.lower, sites.columns)


    # Delineate watersheds and write to output
    with open(watershed_file, 'w') as f:
        for region, sites in sites.groupby('region'):
            print(region)
            nav = Navigator(region, nav_path)
            for _, site in sites.iterrows():
                upstream_reaches, times, message = nav.upstream_watershed(int(site.comid), return_times=True)
                for time in sorted(np.unique(times)):
                    active_reaches = upstream_reaches[np.where(times == time)[0]]
                    out_data = [site.site_id, site.comid, region, time] + list(set(active_reaches))
                    f.write(",".join(map(str, out_data)) + "\n")

main()
