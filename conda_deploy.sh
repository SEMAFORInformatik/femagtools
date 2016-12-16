#!/bin/bash

set -e
cd /tmp

echo "Create skeleton.."
conda skeleton pypi femagtools

echo "Build package.."
conda build femagtools

echo "Converting conda package.."
conda convert -f --platform all $HOME/miniconda*/conda-bld/linux-64/femagtools-*.tar.bz2 -o builds/

echo "Deploying..."
anaconda -t $ANACONDA_TOKEN upload builds/linux-64/femagtools-*.tar.bz2
anaconda -t $ANACONDA_TOKEN upload builds/win-64/femagtools-*.tar.bz2
anaconda -t $ANACONDA_TOKEN upload builds/osx-64/femagtools-*.tar.bz2

echo "Successfully deployed to Anaconda.org."

cd -
exit 0
