#!/usr/bin/env bash

# This script will build the package and updates dependencies by updateing the pyproject.toml file  of the MetObs toolkit

echo " ---- Building and Updating Metobs Toolkit ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)

DISTDIR=${WORKDIR}/dist


#1. cleanup previous builds
rm ${DISTDIR}/*.whl
rm ${DISTDIR}/*.tar.gz

#2. Update the dependencies in the  toml
poetry update

# Toolkit dependencies
poetry add numpy@^1 #v2.x.x conflicting with titanlib
poetry add cartopy@^0.23 #from <= v0.24 requires python >= 3.10 so not valid for py3.9
poetry add earthengine-api@latest
poetry add geemap@latest
poetry add geopandas@latest
poetry add geos@latest
poetry add mapclassify@latest
poetry add matplotlib@latest
poetry add pandas@^2
poetry add shapely@latest


# Toolkit DEV group
poetry add poetry@latest --group dev
poetry add pre-commit@latest --group dev
poetry add poetry-plugin-export --group dev

# Toolkit documentation group
poetry add myst_parser@^3 --group documentation #4.0.x not comp with py3.9
poetry add nbsphinx@latest --group documentation
poetry add pandoc@latest --group documentation
poetry add pydata-sphinx-theme@latest --group documentation
poetry add "sphinx@>=7" --group documentation #v8.x.x not comp with py3.9
poetry add sphinx-copybutton@latest --group documentation
poetry add sphinx-rtd-theme@latest --group documentation
poetry add ipykernel --group documentation #else there is a error when building doc: No such kernel named python3

# Toolkit titan group
# poetry add titanlib@latest --group titan


# 3. update the lock file
#poetry update #to update the poetry.lock with the latest versions of the depending packages
poetry install --all-extras
poetry show


# 4. Make a new build
poetry build
cd ${DEPLOY_DIR}
