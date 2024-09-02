#!/usr/bin/env bash

# This script will run the packgage tests

echo " ---- Running the toolkit-DOC-tests ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)
LOGDIR=${DEPLOY_DIR}/logs
BIGLOGFILE=${LOGDIR}/all_log

#clear big log file
rm -f ${BIGLOGFILE}

#clear all specific log files
rm ${LOGDIR}/DOCtest_output_*_log



#Run Doctests on all modules
echo "RUNNING DOCTESTS NOW"
cd ${WORKDIR}/metobs_toolkit
modules=`ls ./*.py`
cd ${WORKDIR} #call poetry run from in root?
for t in $modules; do
        module_file=${WORKDIR}/metobs_toolkit/${t}
        echo "Running doctests in ${module_file}"

	python3 ${module_file} 2>&1 | tee -a ${LOGDIR}/DOCtest_output_${t:2:-3}_log ${BIGLOGFILE}
done

rm ${WORKDIR}/*.pkl #created by doctest
rm ${WORKDIR}/*.csv #created by docstest
rm ${WORKDIR}/*.html #created by doctest



cd ${DEPLOY_DIR}
