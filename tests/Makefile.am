AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = @JELLYFISH_CFLAGS@ -I$(top_srcdir)/src @AM_CXXFLAGS@
AM_LDFLAGS = @JELLYFISH_LIBS@ -all-static -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -I$(top_srcdir)/src/inc

TESTS = check_comp check_sect
check_PROGRAMS = check_comp check_sect

check_comp_SOURCES = check_comp.cc
check_sect_SOURCES = check_sect.cc
