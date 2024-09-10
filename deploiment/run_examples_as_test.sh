#!/usr/bin/env bash

# This script will run the notebook examples and documentation notebooks as tests

echo " ---- Running notebook examples (from docs) as tests ----"

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}
cd .. #Navigate to workdir
WORKDIR=$(pwd)
LOGDIR=${DEPLOY_DIR}/logs
DOCEXAMPLEDIR=${WORKDIR}/docs/examples
OTHER_NOTEBOOKS=${WORKDIR}/docs/notebook_references
TOPIC_NOTEBOOKS=${WORKDIR}/docs/topics
PAPER_NOTEBOOKS=${WORKDIR}/docs/paper



#make logfile for each test and stream  prompt output for the test
make_test_log () {
        testfile="$1"
        specificlogfile="${LOGDIR}/${testfile::-6}_log"
        rm -f ${specificlogfile}
        touch ${specificlogfile}
        echo "$specificlogfile"
}


#3. Run examples
cd ${DOCEXAMPLEDIR}
filenames=`ls ./*.ipynb`
for t in $filenames; do
        example_file=${DOCEXAMPLEDIR}/${t}
        logfile="$(make_test_log ${t})"
        echo Running ${t} as a test
        #poetry run python ${example_file} >> ${logfile} 2>&1
        #Execute the notebooks, and save the new output (== update them)
        poetry run jupyter nbconvert --execute --to notebook --inplace ${example_file} >> ${logfile} 2>&1
        if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done

#3. Run other notebooks used in the documentation
cd ${OTHER_NOTEBOOKS}
filenames=`ls ./*.ipynb`
for t in $filenames; do
        example_file=${OTHER_NOTEBOOKS}/${t}
        logfile="$(make_test_log ${t})"
        echo Running ${t} as a test
        poetry run jupyter nbconvert --execute --to notebook --inplace ${example_file} >> ${logfile} 2>&1
        if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done


#3. Run other notebooks of topics used in the documentation
cd ${TOPIC_NOTEBOOKS}
filenames=`ls ./*.ipynb`
for t in $filenames; do
        example_file=${TOPIC_NOTEBOOKS}/${t}
        logfile="$(make_test_log ${t})"
        echo Running ${t} as a test
        poetry run jupyter nbconvert --execute --to notebook --inplace ${example_file} >> ${logfile} 2>&1
        if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done


#Run the paper notebook (for figures creation) again with new version
cd ${PAPER_NOTEBOOKS}

filenames=`ls ./*.ipynb`
for t in $filenames; do
        example_file=${PAPER_NOTEBOOKS}/${t}
        logfile="$(make_test_log ${t})"
        echo Running ${t} as a test
        poetry run jupyter nbconvert --execute --to notebook --inplace ${example_file} >> ${logfile} 2>&1
        if [ $? -eq 0 ]; then
                echo "succeeded !!"
        else
                echo "FAIL!!"
        fi

done





# Back to deploy
cd ${DEPLOY_DIR}
