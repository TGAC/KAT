#!/bin/sh -e

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --without-libraries=python,wave,math,graph_parallel,mpi,iostreams
./b2 link=static
cd ../..
