#!/usr/bin/env bash

# Run a full develop pipeline on the toolkit (update + build + doc + tests)

echo " -----------------------------------"
echo " ---- Running DEV pipeline ----"
echo " -----------------------------------"


DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd ${DEPLOY_DIR}

#1. Update and built package
./build_package.sh


#2. Build documentation
./build_documentation.sh

#remove previous logs
LOGDIR=${DEPLOY_DIR}/logs
rm ${LOGDIR}/*_log

#3. Run examples as test
./run_examples_as_test.sh

#4. Run toolkit tests
./run_tests.sh

#5. Run Docstring tests
./run_doctests.sh
