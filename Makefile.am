AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
EXTRA_DIST = README COPYING AUTHORS NEWS

AM_LDFLAGS = -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -I$(top_srcdir)/src/inc
AM_CXXFLAGS = -g -O3


bin_PROGRAMS = $(top_builddir)/bin/kat

__top_builddir__bin_kat_SOURCES = \
                            src/sect/sect.hpp src/sect/sect_main.cc \
                            src/comp/comp_main.cc \
                            src/kat.cc

__top_builddir__bin_kat_LDADD = lib/libjellyfish-1.1.a lib/libz.a


dist_noinst_SCRIPTS = autogen.sh

