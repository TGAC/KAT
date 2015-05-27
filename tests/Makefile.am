AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
CPPFLAGS=-DCPLUSPLUS -I$(top_srcdir)/deps/seqan-library-2.0.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @AM_CPPFLAGS@
CXXFLAGS=-g -O2 -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -std=c++11 @AM_CXXFLAGS@
LDFLAGS=@AM_LDFLAGS@
LIBS=$(top_srcdir)/deps/jellyfish-2.2.0/libjellyfish-2.0.la @BOOST_UNIT_TEST_FRAMEWORK_LIB@ @AM_LIBS@

TESTS_ENVIRONMENT=BOOST_TEST_REPORT_LEVEL=detailed BOOST_TEST_LOG_LEVEL=all
TESTS = check_kat
check_PROGRAMS = check_kat

check_kat_SOURCES = check_comp.cc \
		    check_sect.cc \
		    check_kat.cc
