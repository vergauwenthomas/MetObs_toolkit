# VLINDER toolkit

This repo contains all the software for the [vlinder_toolkit](https://test.pypi.org/project/vlinder-toolkit/).
The vlinder_toolt is a package for scientists who make use of the VLINDER and MOCCA observations. 


## Installing the package
First make shure [GDAL](https://gdal.org/) is installed on your machine. Than install the package by:

`pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple vlinder_toolkit --upgrade`


## Exaples
Some examples with comments can be found in [here](https://github.com/vergauwenthomas/vlinder_toolkit/tree/master/exaples)

## Using the package
* Importing the package: `import vlinder_toolkit`
* Create Settings object: `settings_name = vlinder_toolkit.Settings()`
    * (you can see what is in the settings with `settings_name.show()`
    * (you can update the settings by: `settings_name.update_settings(input_file=..., output_data_folder=...)`
    * (You can look for missing settings by: `settings_name.check_settings()`
* Create a Dataset object: `dataset_name = vlinder_toolkit.Dataset()`
    * There are two ways to import data in the dataset: by .csv file (download using the [brian-tool](https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php) or directly from database (user and password needed).
        *csv-import: `dataset_name.import_data_from_file()` (make sure `input_file` is provided to the settings, see above)
        *database-import: (nesesary to import datetime: `from datetime import datetime` 
            * `dataset_name.import_data_from_database(start_datetime=datetime(2022, 6,12), #2022/6/12 00:00:00
                                    end_datetime=datetime(2022,6,19,12,45)) #2022/7/19 12:45:00`
                                   
    * (To get the dataset in a pandas.Dataframe: `dataset_name.df()`)
    * (To extract one station from the dataset: `station_name = dataset_name.get_station('vlinder05')`)
    * (To make a timeseries plot of a station: `station_name.make_plot(args**)`)
    
 ## Extra info:
 ### Observation names
 Here a list of all possible observationtypes (these names are used):
 `['temp','radiation_temp','humidity','precip','precip_sum','wind_speed','wind_gust','wind_direction','pressure','pressure_at_sea_level']`
 
 These are attributes of a stationobject and can be extracted by `station_name.precip`
 
 ### Metadata names
 Here a list of all metadata per  (these names are used):
 `['network', 'name', 'lat', 'lon', 'call_name', 'location', 'units', 'obs_description']`
 
 These are attributes of a stationobject and can be extracted by `station_name.location`
 
### Using the database
In order to use the database for importing data, you need to have an active VPN connection with the UGent network or working from within the UGent network. 
In addition you need a specific **USER** and **PASSWORD** to connect with the database. (Contact thomas.vergauwen@meteo.be for this account).

To give the user and password to the vlinder toolkit, you need to set them as envrionment variables:
(on linux execute in terminal (or better add them in  `.bashrc`:)

 `     export VLINDER_DB_USER_NAME="...."` (no spaces, fill in the username for the Database)
 
 
 `     export VLINDER_DB_USER_PASW="...."` (no spaces, fill in the password for the Database)
