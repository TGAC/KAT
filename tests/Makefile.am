AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @TEST_CPPFLAGS@ @AM_CPPFLAGS@
AM_CXXFLAGS=-Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
LIBS=@TEST_LIBS@

TESTS = check_comp check_sect
check_PROGRAMS = check_comp check_sect

check_comp_SOURCES = check_comp.cc
check_sect_SOURCES = check_sect.cc
