# VLINDER toolkit

This repo contains all the software for the [vlinder_toolkit](https://test.pypi.org/project/vlinder-toolkit/).
The vlinder_toolt is a package for scientists who make use of the VLINDER and MOCCA observations.
## Documentation ##
Documentation can be found [here](https://vlinder-toolkit.readthedocs.io/en/latest/).

## Installing the package
Install the package by:

`pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple vlinder_toolkit --upgrade`


## Exaples
Some examples with comments can be found in [here](https://github.com/vergauwenthomas/vlinder_toolkit/tree/master/exaples)

Start with the [introduction_example.py](https://github.com/vergauwenthomas/vlinder_toolkit/tree/master/exaples/introduction_example.py).

## Using the package
* Importing the package: `import vlinder_toolkit`
* Create Settings object: `settings_name = vlinder_toolkit.Settings()`
   * Make shure to update the templates (see [templates_example](https://github.com/vergauwenthomas/vlinder_toolkit/tree/master/exaples/templates_example.py)) if you use a new datset.
* Create a Dataset object: `dataset_name = vlinder_toolkit.Dataset()`
* Use the analysing/visualisation tools available in the vlinder_toolkit package.
    
 ## Extra info:
 The templates are used to map the data from the input-file (or sql table) to the vlinder_toolkit-space. Here a short overview on the names used in the vlinder_toolkit-space.
 
 ### Observation names
 Here a list of all possible observationtypes (these names are used):
 `['temp','radiation_temp','humidity','precip','precip_sum','wind_speed','wind_gust','wind_direction','pressure','pressure_at_sea_level']`
 
 These are attributes of a Dataset and are stored in the `Dataset.df` attribute.
 
 ### Metadata names
 Here a list of all metadata per  (these names are used):
 `['network', 'name', 'lat', 'lon', 'call_name', 'location', 'units', 'obs_description']`
 
 These are attributes of a Dataset and are stored in the `Dataset.metadf` attribute.
 
### Using the database
In order to use the database for importing data, you need to have an active VPN connection with the UGent network or working from within the UGent network. 
In addition you need a specific **USER** and **PASSWORD** to connect with the database. (Contact thomas.vergauwen@meteo.be for this account).

To give the user and password to the vlinder toolkit, you need to set them as envrionment variables:
(on linux execute in terminal (or better add them in  `.bashrc`:)

 `     export VLINDER_DB_USER_NAME="...."` (no spaces, fill in the username for the Database)
 
 
 `     export VLINDER_DB_USER_PASW="...."` (no spaces, fill in the password for the Database)
 
 
 # Overview
 ![alt text](https://github.com/vergauwenthomas/vlinder_toolkit/blob/master/examples/overview_fig.png?raw=true)
