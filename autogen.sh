#!/bin/sh -e

# Make sure we are running from the correct directory
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

# Ensure the boost submodule is present
git submodule update --init deps/boost

# Ensure all required boost submodules are present
cd deps/boost
git submodule update --init tools/build/ tools/auto_index/ tools/bcp tools/boostdep tools/inspect tools/litre
git submodule update --init libs/config libs/program_options libs/chrono libs/timer libs/filesystem libs/system libs/stacktrace

# Build boost
./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system,stacktrace
./b2

# Build configure script
cd ../..
autoreconf --force --install --verbose "$srcdir"
