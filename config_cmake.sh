# Auto-Generated on Tue Jun  2 11:37:39 2026
# using library build at /home/clair/repos
# on machine Zelazny1
# settings: --nthread=2

# Setup build and install paths
ROOT_PATH=$(pwd)
BUILD_DIR=$ROOT_PATH/build_release
INSTALL_DIR=$ROOT_PATH/install_release

# Create fresh build directory
rm -rf $BUILD_DIR
mkdir $BUILD_DIR && cd $BUILD_DIR

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DOFT_BUILD_TESTS:BOOL=FALSE \
  -DOFT_PY_KERNEL:STRING=python3 \
  -DOFT_BUILD_EXAMPLES:BOOL=FALSE \
  -DOFT_BUILD_PYTHON:BOOL=TRUE \
  -DOFT_BUILD_DOCS:BOOL=FALSE \
  -DOFT_USE_OpenMP:BOOL=TRUE \
  -DOFT_PACKAGE_BUILD:BOOL=FALSE \
  -DOFT_PACKAGE_NIGHTLY:BOOL=TRUE \
  -DOFT_COVERAGE:BOOL=FALSE \
  -DOFT_DEBUG_STACK:BOOL=FALSE \
  -DOFT_PROFILING:BOOL=FALSE \
  -DOFT_THINCURR_LEGACY:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH=gcc \
  -DCMAKE_CXX_COMPILER:FILEPATH=g++ \
  -DCMAKE_Fortran_COMPILER:FILEPATH=gfortran \
  -DCMAKE_Fortran_FLAGS:STRING="-fallow-argument-mismatch" \
  -DOFT_USE_MPI:BOOL=FALSE \
  -DOFT_METIS_ROOT:PATH=/home/clair/repos/metis-5_1_0 \
  -DHDF5_ROOT:PATH=/home/clair/repos/hdf5-1_14_6 \
  -DBLAS_ROOT:PATH=/home/clair/repos/OpenBLAS-0_3_30 \
  -DLAPACK_ROOT:PATH=/home/clair/repos/OpenBLAS-0_3_30 \
  -DBLA_VENDOR:STRING=OpenBLAS \
  -DLIBXML2_ROOT:PATH=/home/clair/repos/libxml2-v2_15_2 \
  /home/clair/repos/OpenFUSIONToolkit/src
