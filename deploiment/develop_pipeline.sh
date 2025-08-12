#!/usr/bin/env bash

# Run a full develop pipeline on the toolkit (update + build + doc + tests + code cleanup)

echo " -----------------------------------"
echo " ---- Running DEV pipeline ----"
echo " -----------------------------------"


DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPODIR=${DEPLOY_DIR}/..

cd ${DEPLOY_DIR}



#2. Update and built package
echo "Updating the package and dependecies"
rm ${REPODIR}/dist/*.whl
rm ${REPODIR}/dist/*.tar.gz

poetry update
poetry show
poetry build
poetry install --all-extras --no-root


#1. Cleanup the code (run black, black config in pyproject file)
echo "Running black"
BLACK_LOG=${DEPLOY_DIR}/black_log.log
rm -f ${BLACK_LOG} #clean start
touch ${BLACK_LOG}
cd $REPODIR
poetry run black . 2>&1 | tee -a ${BLACK_LOG}


#3. Run tests
echo "Running testing framework on the tests"
TEST_LOG=${DEPLOY_DIR}/pytest_tests_log.log
rm -f ${TEST_LOG} #clean start
touch ${TEST_LOG}
cd ${REPODIR}/tests
poetry run pytest . --mpl --mpl-generate-summary=html 2>&1 | tee -a ${TEST_LOG}

#4. Run notebook example as tests
echo "Running doc notebooks as test"
NB_LOG=${DEPLOY_DIR}/pytest_on_doc_notebooks_log.log
rm -f ${NB_LOG} #clean start
touch ${NB_LOG}
cd ${REPODIR}/docs/examples
poetry run pytest . --nbval-lax 2>&1 | tee -a ${NB_LOG}
cd ${REPODIR}/docs/topics
poetry run pytest . --nbval-lax 2>&1 | tee -a ${NB_LOG}


#5. Build documentation
echo "Building documentation"
DOCS_LOG=${DEPLOY_DIR}/build_doc_log.log
rm -f ${DOCS_LOG} #clean start
touch ${DOCS_LOG}
cd ${REPODIR}
cd docs
poetry run ./build_doc 2>&1 | tee -a ${DOCS_LOG}



#6. Create a big log file

BIG_LOG=${DEPLOY_DIR}/dev_pipeline_full_log.log
rm -f ${BIG_LOG} # clean start
touch ${BIG_LOG}

cat ${BLACK_LOG} ${TEST_LOG} ${NB_LOG} ${DOCS_LOG} >> ${BIG_LOG}


echo "open logs in geany"
geany ${BIG_LOG} &

cd ${REPODIR}


echo " -----------------------------------"
echo " ---- DONE RUNNING ----"
echo " -----------------------------------"
