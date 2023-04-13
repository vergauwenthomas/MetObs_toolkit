# MetObs-toolkit

This repo contains all the software for the [metobs_toolkit](https://test.pypi.org/project/metobs-toolkit/).
The MetObs-toolkit is a package for scientists who make use meteorological observations.
## Documentation ##
Documentation can be found [here](https://vergauwenthomas.github.io/vlinder_toolkit/).

## Installing the package
Install the package by:

`pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple metobs_toolkit --upgrade`


### Using the database
In order to use the database for importing data, you need to have an active VPN connection with the UGent network or working from within the UGent network.
In addition you need a specific **USER** and **PASSWORD** to connect with the database. (Contact thomas.vergauwen@meteo.be for this account).

To give the user and password to the vlinder toolkit, you need to set them as envrionment variables:
(on linux execute in terminal (or better add them in  `.bashrc`:)

 `     export VLINDER_DB_USER_NAME="...."` (no spaces, fill in the username for the Database)


 `     export VLINDER_DB_USER_PASW="...."` (no spaces, fill in the password for the Database)


 # Overview
 ![alt text](https://github.com/vergauwenthomas/vlinder_toolkit/blob/master/examples/overview_fig.png?raw=true)
