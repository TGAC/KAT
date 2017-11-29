#!/bin/sh -e

# Make sure we are running from the correct directory
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Ensure all required boost submodules are present
cd deps/boost

# Build boost
./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system
./b2

# Build configure script
cd ../..
autoreconf --force --install --verbose "$srcdir"
