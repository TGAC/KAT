#!/bin/sh -e

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --without-libraries=python
./b2 link=static
cd ../..
