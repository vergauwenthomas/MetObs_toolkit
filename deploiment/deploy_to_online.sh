#!/usr/bin/env bash

DEPLOY_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )



#This script will publish the current build. BE AWARE THAT EXECUTION THIS SCRIPT CAN NOT BE UNDONE!!!!

#if an error occurs that the file already exists, than update the version in the pyproject.

#1 ------  update versions --------
# update version using poetry --> this will update the pyproject.toml
poetry version prerelease
# extract version from pyproject to variable
cd ..
version=$(awk '/^version/{print $NF}' pyproject.toml)

# update the __version__ in the init file
cd metobs_toolkit
sed -i "s/.*__version__=.*/__version__=${version}/g" __init__.py
echo update __init__ version to ${version}


#2 ------- build package ----------


cd ${DEPLOY_DIR}
cd .. #Navigate to workdir


WORKDIR=$(pwd)
DOCDIR=${WORKDIR}/docs
DISTDIR=${WORKDIR}/dist
TESTDIR=${WORKDIR}/tests
EXAMPLESDIR=${WORKDIR}/examples

#1. Update documentation
cd ${DOCDIR}
source build_doc
cd ${DEPLOY_DIR}

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





#3 ------- Upload packgage ---------
#TESTPIPY credentials

USERNAME=__token__
PASSWORD_testpypi=pypi-AgENdGVzdC5weXBpLm9yZwIkYjk3Njk1NTUtYTY0Yy00NjU2LTgwNDItMmFjOGU5NDlkYjI2AAIWWzEsWyJtZXRvYnMtdG9vbGtpdCJdXQACLFsyLFsiMWYzN2IyNDQtYzhlZi00ZjU0LThmMzAtNjIyMmM3MDlkNDNhIl1dAAAGIHgaTgkg0zT_Luf6AaCez1Yj1VVBQaZMR8VqRj0WrfHO
PASSWORD_pypi=pypi-AgEIcHlwaS5vcmcCJDI4N2U1ZWI0LTE4MjctNGQ3Ni04YzMxLTMzMWUzY2UyMjVhNwACFlsxLFsibWV0b2JzLXRvb2xraXQiXV0AAixbMixbImU1YjBmNDMxLTgzM2ItNGQ0Yy05Yzk3LTJjYzc1MTY1MWYwNyJdXQAABiDLqpyqFXsYkX19gfAZnm5QFPTO--yn07Y-J9TxLWglEA


#To testpypi
poetry publish --repository testpypi -u ${USERNAME} -p ${PASSWORD_testpypi}


#To pypi
poetry config pypi-token.pypi ${PASSWORD_pypi}
poetry publish -u ${USERNAME} -p ${PASSWORD_pypi}


