import os

#################################################################
### Adjust run id and input shapefile for each individual run ###
#################################################################
run_id = "cww_june18"

# Shapefile containing the locations of the sample sites and the field containing the site id
points_file = os.path.join("..", "bin", "Sites", "cws_june18.shp")
site_id_field = "PWSID"

##################################################################


""" Paths which can stay fixed between runs """

""" These paths are used as inputs to the tool and must point to existing content """

# Path containing NHD Plus dataset
nhd_dir = r"C:\Users\jhook\Documents\NationalData\NHDPlusV2"
catchment_path = os.path.join(nhd_dir, "NHDPlus{}", "NHDPlusCatchment", "Catchment.shp")

# NHD region boundaries
regions_file = r"C:\Users\jhook\Documents\NationalData\NHD_Regions\all_regions.shp"
region_id_field = "REGION"

# Location of the weather files
metfile_path = r"C:\Users\Jhook\Documents\opp-efed\sam\bin\Preprocessed\WeatherFiles\metfile"

# Shapefile containing weather file grids
grid_file = r"C:\Users\jhook\Documents\NationalData\WeatherFiles\weather_stations_highres_thiessens_US_alb\weather_stations_highres_thiessens_US_alb.shp"
grid_id_field = "stationID"

# Location of the Navigator files used to delineate watersheds
nav_path = os.path.join(r"C:\Users\jhook\Documents\opp-efed\sam", "bin", "Preprocessed", "HydroFiles")

# Scratch workspace
scratch_workspace = os.path.join("..", "bin", "temp")

""" Files that are generated as the tool runs.  These do not need to exist ahead of time """

# Table containing the NHD reach IDs for each sample site
sites_file = os.path.join("..", "bin", "Tables", "{}_overlay.txt".format(run_id))

# Table containing all reaches upstream of each sample site
watershed_file = os.path.join("..", "bin", "Tables", "{}_watersheds.csv".format(run_id))

# Table with data on the overlay of the watersheds and weather grids
grid_props = os.path.join("..", "bin", "Tables", "{}_grid_props.csv".format(run_id))

# Output precip file
outfile = os.path.join("..", "bin", "Tables", "{}_precip.txt".format(run_id))
