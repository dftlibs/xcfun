#    # load Intel compilers
#    set +x
#    source /theoryfs2/common/software/intel2018/bin/compilervars.sh intel64
#    set -x
# error with Intel 2018.3, 2019.4
#
#    ALLOPTS="-gnu-prefix=${HOST}- ${OPTS}"

# configure
${BUILD_PREFIX}/bin/cmake \
        -H${SRC_DIR} \
        -Bbuild \
        -DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DENABLE_FC_SUPPORT=ON \
        -DCMAKE_Fortran_COMPILER=${FORTRAN} \
        -DPYTHON_EXECUTABLE="${PYTHON}" \
        -DENABLE_PYTHON_INTERFACE=ON \
        -DCMAKE_PREFIX_PATH="${PREFIX}"

# build
cd build
make -j${CPU_COUNT}

# test
# The Python interface is tested using pytest directly
ctest -E "python-interface" -j${CPU_COUNT} --output-on-failure --verbose

# install
make install
