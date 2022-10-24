#!/bin/sh

#step 0 update apt
apt update

#step 1 install GEOS
apt-get install libgeos-dev

#step 2 install proj dependencies
apt install cmake
apt install sqlite3
apt install curl && apt-get install libcurl4-openssl-dev

#step 3 install Proj
#Unfortunately, cartopy requires proj v8.0.0 as a minimum, but if you install proj using apt you can only install proj v6.3.1

# Just for reference in case anything changes, this is the command to install proj from apt:
#apt-get install proj-bin
#apt-get install proj-bin libproj-dev proj-data

#Else build proj from source

cd proj
tar -xf proj-9.0.0.tar.gz
cd proj-9.0.0
mkdir build && cd build

cmake ..
cmake --build .
cmake --build . --target install

#make test
ctest

#move binaries toe the /bin dir
cp ./bin/* /bin
cp ./lib/* /lib


