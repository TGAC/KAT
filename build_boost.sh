#!/bin/sh -e

# Using BCP
#bcp algorithm build chrono exception filesystem locale program_options stacktrace system timer deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --prefix=build --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 --ignore-site-config headers
./b2 --ignore-site-config install
cd ../..
