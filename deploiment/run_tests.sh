#!/usr/bin/env bash

# This script will run the packgage tests

echo " ---- Running the toolkit-tests ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)
LOGDIR=${DEPLOY_DIR}/logs
TESTDIR=${WORKDIR}/tests


#make logfile for each test and stream  prompt output for the test
make_test_log () {
        testfile="$1"
        specificlogfile="${LOGDIR}/${testfile::-3}_log"
        rm -f ${specificlogfile}
        touch ${specificlogfile}
        echo "$specificlogfile"
}


#1. Run tests
#Run tests
cd ${TESTDIR}/push_test
filenames=`ls ./*.py`
for t in $filenames; do
        push_file=${TESTDIR}/push_test/${t}
        logfile="$(make_test_log ${t})"
        echo Running push tests: ${t}
        poetry run python ${push_file} >> ${logfile} 2>&1
        if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done



cd ${DEPLOY_DIR}
