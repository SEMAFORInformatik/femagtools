# Only need to change these two variables
PKG_NAME=femagtools
USER=semafor

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
for v in 2.7 3.5 3.6
do
  conda-build --python $v .
done
conda convert --platform all
for os in linux-64 linux-32 osx-64 win-32 linux-aarch64 linux-armv6l linux-armv7l linux-ppc64le
do
    anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $CONDA_BLD_PATH/$os/*
done
