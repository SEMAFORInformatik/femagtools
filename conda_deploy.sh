#!/bin/bash

set -e
cd /tmp

echo "Create skeleton.."
conda skeleton pypi femagtools
export VERSION=`date +%Y.%m.%d`
echo "Build package.."
dest_file=$(conda build femagtools --output)
dest_dir=$(dirname $dest_file)
conda build femagtools
echo "File located under $dest_dir"

echo "Converting conda package.."
conda convert -f --platform all $dest_dir/femagtools-*.tar.bz2 -o builds

echo "Deploying..."
anaconda -t $ANACONDA_TOKEN upload builds/**/femagtools-*.tar.bz2
anaconda -t $ANACONDA_TOKEN upload $dest_dir

echo "Successfully deployed to Anaconda.org."

cd -
exit 0
