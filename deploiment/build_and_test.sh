#!/usr/bin/env bash

#This script will:
 # run

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir


WORKDIR=$(pwd)
DISTDIR=${WORKDIR}/dist
TESTDIR=${WORKDIR}/tests 
EXAMPLESDIR=${WORKDIR}/examples

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

#make logfile for each test and stream  prompt output for the test
make_test_log () {
	testfile="$1"
	specificlogfile="${DEPLOY_DIR}/${testfile::-3}_log"
	rm -f ${specificlogfile}
	touch ${specificlogfile}
	echo "$specificlogfile"
}



#Run examples
cd ${EXAMPLESDIR}
filenames=`ls ./*.py`
for t in $filenames; do
	example_file=${EXAMPLESDIR}/${t}
	logfile="$(make_test_log ${t})"
	echo Running ${t} as a test
	poetry run python ${example_file} >> ${logfile} 2>&1
	if [ $? -eq 0 ]; then
    		echo "succeeded !!"
	else
    		echo "FAIL!!"
	fi

done




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
echo "Testing Done"
