#!/bin/bash

set -e
cd /tmp

# Upload the linux-64 to anaconda automaticly
echo "Set config"
conda config --set anaconda_upload yes

echo "Create skeleton.."
conda skeleton pypi femagtools

echo "Build package.."
dest_dir=`conda build femagtools --output`
echo "File located under $dest_dir"

echo "Converting conda package.."
conda convert -f --platform all $dest_dir -o builds

echo "Deploying..."
anaconda -t $ANACONDA_TOKEN upload builds/*/femagtools-*.tar.bz2
anaconda -t $ANACONDA_TOKEN upload $dest_dir

echo "Successfully deployed to Anaconda.org."

cd -
exit 0
