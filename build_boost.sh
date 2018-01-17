#!/bin/sh -e

# Where the boost include and lib files get sent to
INST_DIR=../boost_build

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --prefix=${INST_DIR} --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 --prefix=${INST_DIR} headers
./b2 --prefix=${INST_DIR} install
cd ../..
