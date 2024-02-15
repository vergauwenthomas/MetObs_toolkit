#!/usr/bin/env bash

#This script will:
 # run

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir


WORKDIR=$(pwd)
DOCDIR=${WORKDIR}/docs
DISTDIR=${WORKDIR}/dist
TESTDIR=${WORKDIR}/tests
DOCEXAMPLEDIR=${WORKDIR}/docs/examples


#1 install the package using poetry

poetry update #to update the poetry.lock with the latest versions of the depending packages
poetry install --all-extras



#list all packages installed (for debugging)
poetry show


#2. build documentation
cd ${DOCDIR}
source build_doc
cd ${DEPLOY_DIR}


#3 build the package
echo "Removing old builds before building the package ..."
cd ${DISTDIR}
rm *.whl
rm *.tar.gz


echo "Make build ..."
cd ${WORKDIR}
poetry build



#echo "Export requirements file ..."
#poetry export -f requirements.txt --output requirements.txt



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





#Run examples included in the documentation
echo 'Running the documentation examples as test'
cd ${DOCEXAMPLEDIR}
#delete all .py versions of the examples (rebuild them from the notebooks)
rm ${DOCEXAMPLEDIR}/*.py
#convert nb to python files
jupyter nbconvert --to python *.ipynb

cd ${DOCEXAMPLEDIR}
filenames=`ls ./*.py`
for t in $filenames; do
	example_file=${DOCEXAMPLEDIR}/${t}
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
	echo Running push tests: ${t}g
	poetry run python ${push_file} >> ${logfile} 2>&1
	if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done


#Run Doctests on all modules
echo "RUNNING DOCTESTS NOW"
cd ${WORKDIR}/metobs_toolkit
modules=`ls ./*.py`
for t in $modules; do
	module_file=${WORKDIR}/metobs_toolkit/${t}
	python3 -m doctest -o ELLIPSIS -o NORMALIZE_WHITESPACE ${module_file}
done

rm ${WORKDIR}/metobs_toolkit/*.pkl #created by doctest
rm ${WORKDIR}/metobs_toolkit/*.csv #created by docstest



cd ${DEPLOY_DIR}
echo "Testing Done"
