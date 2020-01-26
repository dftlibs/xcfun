#!/usr/bin/env bash

set -eu -o pipefail

example=$1
lang="${example/_*/}"
var=${lang}_COMPILER
compiler=${var}
echo "-- Running example $example"
(
    cd examples/"$example"
    cmake -H. -Bbuild -DXCFun_DIR="$HOME"/Software/xcfun/share/cmake/XCFun -DCMAKE_"${lang}"_COMPILER="${!compiler}"
    cmake --build build
    cd build
    ctest
)
echo "-- Done with example $example"
