#!/usr/bin/env bash

# update version using poetry --> this will update the pyproject.toml
poetry version prerelease


# extract version from pyproject to variable
cd ..
version=$(awk '/^version/{print $NF}' pyproject.toml)


# update the __version__ in the init file
cd metobs_toolkit
sed -i "s/.*__version__=.*/__version__=${version}/g" __init__.py
echo update __init__ version to ${version}

# navigate back
cd ../deploiment
