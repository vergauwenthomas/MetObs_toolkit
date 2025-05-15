#!/usr/bin/env bash

# Run a full develop pipeline on the toolkit (update + build + doc + tests + code cleanup)

echo " -----------------------------------"
echo " ---- Running DEV pipeline ----"
echo " -----------------------------------"


DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DEP_PIPE_LOG=${DEPLOY_DIR}/dev_pipeline_log
REPODIR=${DEPLOY_DIR}/..

cd ${DEPLOY_DIR}

rm -f ${DEP_PIPE_LOG} #clean start
touch ${DEP_PIPE_LOG}

#remove previous logs
LOGDIR=${DEPLOY_DIR}/logs
rm ${LOGDIR}/*_log


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
cd $REPODIR
poetry run black . 2>&1 | tee -a ${DEP_PIPE_LOG}


#3. Run tests
echo "Running testing framework"
cd ${REPODIR}/toolkit_tests
poetry run pytest --mpl --mpl-generate-summary=html 2>&1 | tee -a ${DEP_PIPE_LOG}


#4. Run notebook example as tests
echo "Running testing framework"
cd ${REPODIR}/docs/examples
poetry run pytest --nbval-lax 2>&1 | tee -a ${DEP_PIPE_LOG}
cd ${REPODIR}/docs/topics
poetry run pytest --nbval-lax 2>&1 | tee -a ${DEP_PIPE_LOG}

#5. Build documentation
#./build_documentation.sh 2>&1 | tee -a ${DEP_PIPE_LOG}
echo "Building documentation"
cd ${REPODIR}
cd docs
poetry run ./build_doc 2>&1 | tee -a ${DEP_PIPE_LOG}



echo "open logs in geany"
geany ${DEP_PIPE_LOG} &



echo " -----------------------------------"
echo " ---- DONE RUNNING ----"
echo " -----------------------------------"
