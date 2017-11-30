#!/bin/sh -e

# Make sure we are running from the correct directory
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh
./b2
cd ../..

# Create configure script
autoreconf --force --install --verbose "$srcdir"
