#!/bin/sh -e

# Make sure we are running from the correct directory
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

# Create configure script
autoreconf --force --install --verbose "$srcdir"
