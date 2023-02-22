#!/usr/bin/env bash

#This script will:
 # run

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir


WORKDIR=$(pwd)
DISTDIR=${WORKDIR}/dist
TESTDIR=${WORKDIR}/tests 


#1. Update the pyproject.toml for any depending packages



#2 install the package using poetry

poetry update #to update the poetry.lock with the latest versions of the depending packages 
poetry install 


#3 build the package

echo "Removing old builds before building the package ..."
cd ${DISTDIR}
rm *.whl
rm *.tar.gz



echo "Make build ..."
cd ${WORKDIR}
poetry build


#4 testing
#run tests in the poetry environment
cd ${TESTDIR}
poetry run python package_import_test.py 

cd ${DEPLOY_DIR}
echo "Testing Done"
