#!/usr/bin/env bash

# Run a full develop pipeline on the toolkit (update + build + doc + tests)

echo " -----------------------------------"
echo " ---- Running DEV pipeline ----"
echo " -----------------------------------"


DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DEP_PIPE_LOG=${DEPLOY_DIR}/dev_pipeline_log

cd ${DEPLOY_DIR}

rm -f ${DEP_PIPE_LOG} #clean start
touch ${DEP_PIPE_LOG}

#1. Update and built package
./build_package.sh 2>&1 | tee -a ${DEP_PIPE_LOG}


#remove previous logs
LOGDIR=${DEPLOY_DIR}/logs
rm ${LOGDIR}/*_log

#2. Run examples as test
./run_examples_as_test.sh 2>&1 | tee -a ${DEP_PIPE_LOG}

#3. Build documentation
./build_documentation.sh 2>&1 | tee -a ${DEP_PIPE_LOG}


#4. Run toolkit tests
./run_tests.sh 2>&1 | tee -a ${DEP_PIPE_LOG}

#5. Run Docstring tests
./run_doctests.sh 2>&1 | tee -a ${DEP_PIPE_LOG}

echo " -----------------------------------"
echo " ---- DONE RUNNING ----"
echo " -----------------------------------"
