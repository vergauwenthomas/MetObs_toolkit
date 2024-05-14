#!/usr/bin/env bash

# This script will run the packgage tests

echo " ---- Running the toolkit-DOC-tests ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)


#Run Doctests on all modules
echo "RUNNING DOCTESTS NOW"
cd ${WORKDIR}/metobs_toolkit
modules=`ls ./*.py`
cd ${WORKDIR} #call poetry run from in root?
for t in $modules; do
        module_file=${WORKDIR}/metobs_toolkit/${t}
	poetry run python3 -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE ${module_file}
done

rm ${WORKDIR}/metobs_toolkit/*.pkl #created by doctest
rm ${WORKDIR}/metobs_toolkit/*.csv #created by docstest

cd ${DEPLOY_DIR}
