#!/usr/bin/env bash

# This script will build the documentation of the MetObs toolkit

echo " ---- Building Documentation ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)



#1. Setup environment
poetry install --all-extras

#2. Cleanup builds
DOCDIR=${WORKDIR}/docs
#clear builds

rm -r ${DOCDIR}/_build/*
rm -r ${DOCDIR}/_autosummary/*

#3. Build documentation
cd ${WORKDIR}
sphinx-build -a -E -v docs/ docs/_build/

#4. Back to deploy dir
cd ${DEPLOY_DIR}
