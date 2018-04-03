import os
import pandas as pd
from ftplib import FTP


class StreamCat(object):
    def __init__(self, variables, output_dir, overwrite=False):
        self.vars = variables
        self.out_dir = output_dir
        self.overwrite = overwrite

        streamcat_files = self.assemble_tables()

        self.condense_tables(streamcat_files)


    def assemble_tables(self):
        from utilities import nhd_states

        # Make FTP connection.  FTP site for StreamCat is
        # ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/HydroRegions/
        ftp = FTP('newftp.epa.gov')  # connect to host, default port
        ftp.login()  # user anonymous, passwd anonymous@
        ftp.cwd('EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/HydroRegions')

        # Initialize file storage
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # Loop through region and variable and pull files
        file_index = []
        for variable in self.vars:
            for region in nhd_states:
                dataset = "{}_Region{}.zip".format(variable, region)
                out_dataset = os.path.join(self.out_dir, dataset)
                if self.overwrite or not os.path.exists(out_dataset):
                    print("Pulling dataset {} from EPA site...".format(dataset))
                    ftp.retrbinary('RETR {}'.format(dataset), open(out_dataset, 'wb').write)
                file_index.append((variable, region, out_dataset))
        return file_index

    def condense_tables(self, file_index):
        for variable, region, dataset in file_index:
            dataset = pd.read_csv(dataset)
            print(dataset)
            input()


def main():
    from utilities import streamcat_tables

    # Set output directory and variables of interest
    output_dir = os.path.join("..", "bin", "StreamCat")

    # Do it
    StreamCat(streamcat_tables, output_dir)


main()