#!/bin/sh -e

# Where the boost include and lib files get sent to
INST_DIR=build

# Git method
# Ensure the boost submodule is present
#git submodule sync
#git submodule update --recursive --init deps/boost

# Using BCP
#bcp algorithm build chrono exception filesystem locale program_options stacktrace system timer deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --prefix=${INST_DIR} --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 headers
./b2 install
cd ../..
