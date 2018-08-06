import os

""" Adjust run id and input shapefile for each individual run """

run_id = "csw_june18"

# Shapefile containing the locations of the sample sites
points_file = os.path.join("..", "bin", "Sites", "cws_june18.shp")

""" Paths which can stay fixed between runs """

# Path containing NHD Plus dataset
nhd_dir = r"C:\Users\Trip Hook\Documents\NationalData\NHDPlusV2"

# Location of the weather files
metfile_path = os.path.join(r"C:\Users\jhook\Documents\opp-efed\sam", "bin", "Tables", "WeatherFiles", "metfile")

# Shapefile containing weather file grids
grid_file = r"C:\Users\Trip Hook\Documents\BiasFactors\Data\weather_stations_highres_thiessens_US_alb\weather_stations_highres_thiessens_US_alb.shp"

# Location of the Navigator files used to delineate watersheds
nav_path = os.path.join(r"C:\Users\Trip Hook\Documents\opp-efed\sam", "bin", "Preprocessed", "HydroFiles")

# Table containing the NHD reach IDs for each sample site
sites_file = r"..\Data\{}_overlay.txt".format(run_id)

# Table containing all reaches upstream of each sample site
watershed_file = os.path.join("..", "bin", "Tables", "{}_watersheds.csv".format(run_id))

# Table with data on the overlay of the watersheds and weather grids
grid_props = os.path.join("..", "bin", "Tables", "{}_grid_props.csv".format(run_id))

# Output precip file
outfile = os.path.join("..", "bin", "Tables", "{}_precip.txt".format(run_id))
