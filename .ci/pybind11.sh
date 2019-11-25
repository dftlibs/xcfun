#!/usr/bin/env bash

set -eu -o pipefail

pybind11_VERSION="2.4.3"
echo "-- Installing pybind11 $pybind11_VERSION"
cd "$HOME"/Downloads
mkdir -p pybind11
curl -Ls https://github.com/pybind/pybind11/archive/v${pybind11_VERSION}.tar.gz | tar -xz -C pybind11 --strip-components=1
cd pybind11
cmake -H. -Bbuild_pybind11 -DCMAKE_INSTALL_PREFIX="$HOME"/Deps/pybind11 -DPYTHON_EXECUTABLE=$(which python) -DPYBIND11_TEST=OFF #&> /dev/null
cmake --build build_pybind11 -- install #&> /dev/null
cd "$TRAVIS_BUILD_DIR"
rm -rf "$HOME"/Downloads/pybind11
echo "-- Done with pybind11 $pybind11_VERSION"
