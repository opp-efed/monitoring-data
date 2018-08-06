import numpy as np
import pandas as pd

from collections import defaultdict

from classes import Navigator, nhd_states


class Sites(object):
    def __init__(self, sites_file, nav_dir):
        self.file = sites_file
        self.nav_dir = nav_dir
        self.sites = self.get_sites()
        self.get_nhd_regions()

    def get_sites(self):
        with open(self.file) as f:
            next(f)
            sites = [line.strip().split(",") for line in f]
        return pd.DataFrame(sites, columns=('site_id', 'comid', 'state'))

    def get_nhd_regions(self):
        # Get the potential list of NHD regions for each site
        regions = defaultdict(list)
        for _, site in self.sites.iterrows():
            for nhd_region, states in nhd_states.items():
                if site.state in states:
                    regions[nhd_region].append((site.site_id, site.comid))

        # Check each region for sites
        self.sites['region'] = None
        for region, sites in regions.items():
            nav = Navigator(region, self.nav_dir)
            for site_id, comid in sites:
                if int(comid) in nav.reach_ids:
                    self.sites.loc[self.sites.site_id == site_id, 'region'] = region
                else:
                    print("{} not in watershed".format(comid))


def main():
    from paths import sites_file, nav_path, watershed_file

    # Read sites file
    sites = Sites(sites_file, nav_path)

    # Delineate watersheds and write to output
    with open(watershed_file, 'w') as f:
        for region, sites in sites.sites.groupby('region'):
            print(region)
            nav = Navigator(region, nav_path)
            for _, site in sites.iterrows():
                upstream_reaches, times, message = nav.upstream_watershed(int(site.comid), return_times=True)
                for time in sorted(np.unique(times)):
                    active_reaches = upstream_reaches[np.where(times == time)[0]]
                    f.write("{},{},{},{},{}\n".format(site.site_id, site.comid, region, time,
                                                      ",".join(map(str, active_reaches))))


main()
