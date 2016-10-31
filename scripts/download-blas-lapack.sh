#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "usage: <build directory> [FLAGS]";
  exit 1;
fi

CLONEDIR=$1/openblas
OPENBLAS_LOC=https://github.com/xianyi/OpenBLAS.git
git clone ${OPENBLAS_LOC} ${CLONEDIR}

pushd ${CLONEDIR}
make -j "${@:2}"
make PREFIX=$1/lapack install
popd

