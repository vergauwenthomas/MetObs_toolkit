#!/usr/bin/env bash

# Run a full develop pipeline on the toolkit (update + build + doc + tests + code cleanup)

echo " -----------------------------------"
echo " ---- Running DEV pipeline ----"
echo " -----------------------------------"


DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPODIR=${DEPLOY_DIR}/..

cd ${DEPLOY_DIR}

#0. Set up GEE authentication environment variable
echo "========================================="
echo "Step 0: Setting up GEE authentication..."
echo "========================================="

# Check for GEE credentials and set environment variable
GEE_CRED_PATH="${HOME}/.config/earthengine/credentials"
if [ -f "$GEE_CRED_PATH" ]; then
    export EARTHENGINE_CREDENTIALS_PATH="$GEE_CRED_PATH"
    echo "✓ GEE credentials found at: $GEE_CRED_PATH"
    echo "  EARTHENGINE_CREDENTIALS_PATH environment variable set"
else
    echo "⚠ WARNING: GEE credentials not found at expected location"
    echo "  Expected: $GEE_CRED_PATH"
    echo "  Tests requiring GEE may fail"
    echo ""
    echo "  To set up GEE credentials, run:"
    echo "    python -c 'import metobs_toolkit; metobs_toolkit.connect_to_gee()'"
    echo ""
fi
echo ""


#1. Update and build package
echo "========================================="
echo "Step 1: Updating package and dependencies..."
echo "========================================="
rm ${REPODIR}/dist/*.whl
rm ${REPODIR}/dist/*.tar.gz

poetry update
poetry show
poetry build
poetry install --all-extras --no-root
echo ""


#1a. Test GEE authentication in Poetry environment (optional)
if [ "$TEST_GEE_AUTH" = "1" ]; then
    echo "========================================="
    echo "Step 1a: Testing GEE authentication in Poetry environment..."
    echo "========================================="
    cd ${DEPLOY_DIR}
    poetry run python ${DEPLOY_DIR}/test_gee_auth.py
    if [ $? -ne 0 ]; then
        echo ""
        echo "⚠ GEE authentication test failed in Poetry environment!"
        echo "  Continue anyway? (tests may fail)"
        read -p "Press Enter to continue or Ctrl+C to abort..."
    fi
    echo ""
fi


#2. Cleanup the code (run black, black config in pyproject file)
echo "========================================="
echo "Step 2: Running black code formatter..."
echo "========================================="
BLACK_LOG=${DEPLOY_DIR}/black_log.log
rm -f ${BLACK_LOG} #clean start
touch ${BLACK_LOG}
cd $REPODIR
poetry run black . 2>&1 | tee -a ${BLACK_LOG}
echo ""


#3. Run tests
echo "========================================="
echo "Step 3: Running test suite..."
echo "========================================="
TEST_LOG=${DEPLOY_DIR}/pytest_tests_log.log
rm -f ${TEST_LOG} #clean start
touch ${TEST_LOG}
cd ${REPODIR}/tests
poetry run pytest . --mpl --mpl-generate-summary=html 2>&1 | tee -a ${TEST_LOG}
echo ""

#4. Run notebook example as tests
echo "========================================="
echo "Step 4: Running notebook examples as tests..."
echo "========================================="
NB_LOG=${DEPLOY_DIR}/pytest_on_doc_notebooks_log.log
rm -f ${NB_LOG} #clean start
touch ${NB_LOG}
cd ${REPODIR}/docs/examples
poetry run pytest . --nbval-lax 2>&1 | tee -a ${NB_LOG}
cd ${REPODIR}/docs/topics
poetry run pytest . --nbval-lax 2>&1 | tee -a ${NB_LOG}
echo ""


#5. Build documentation
echo "========================================="
echo "Step 5: Building documentation..."
echo "========================================="
DOCS_LOG=${DEPLOY_DIR}/build_doc_log.log
rm -f ${DOCS_LOG} #clean start
touch ${DOCS_LOG}
cd ${REPODIR}
cd docs
poetry run ./build_doc 2>&1 | tee -a ${DOCS_LOG}
echo ""


#6. Create a big log file
echo "========================================="
echo "Step 6: Creating combined log file..."
echo "========================================="
BIG_LOG=${DEPLOY_DIR}/dev_pipeline_full_log.log
rm -f ${BIG_LOG} # clean start
touch ${BIG_LOG}

cat ${BLACK_LOG} ${TEST_LOG} ${NB_LOG} ${DOCS_LOG} >> ${BIG_LOG}


echo "Opening logs in default text viewer..."
xdg-open ${BIG_LOG} &

cd ${REPODIR}


echo " -----------------------------------"
echo " ---- DONE RUNNING ----"
echo " -----------------------------------"
