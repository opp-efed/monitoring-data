import os
import pandas as pd


class JoinedTable(object):
    def __init__(self, region, streamcat_dir, warp_dir, curve_number_dir, nhd_dir, output_dir, overwrite=True):
        self.region = region
        self.outfile = output_dir.format(region)

        if overwrite or not os.path.exists(self.outfile):
            self.matrix = self.read_streamcat(streamcat_dir)
            self.read_warp()
            self.read_curve_number()
            self.read_nhd()
            self.matrix.to_csv(self.outfile)

    def read_streamcat(self, table_dir):
        from utilities import streamcat_tables

        master_table = None
        for table in streamcat_tables:
            streamcat_file = os.path.join(table_dir, "{}_Region{}.csv".format(table, self.region))
            if os.path.exists(streamcat_file):
                streamcat_file = pd.read_csv(streamcat_file)
                if master_table is not None:
                    master_table = pd.concat([master_table, streamcat_file], axis=0)
                else:
                    master_table = streamcat_file
            else:
                print("StreamCat table {} not found for Region {}".format(table, self.region))

        return master_table


def main():
    from utilities import nhd_states

    streamcat_dir = os.path.join("..", "bin", "streamcat")
    warp_dir = os.path.join("..", "bin", "warp_inputs")
    curve_number_dir = os.path.join("..", "curve_number")
    nhd_dir = os.path.join(r"C:\"", "Documents", "NationalData", "NHDPlusV2")
    output_dir = os.path.join("..", "bin", "Results")

    for region in nhd_states:
        master_table = JoinedTable(region, streamcat_dir, warp_dir, curve_number_dir, nhd_dir, output_dir)


main()
