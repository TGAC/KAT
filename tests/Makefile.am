AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
CPPFLAGS=-DCPLUSPLUS -I$(top_srcdir)/deps/seqan-library-2.0.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @AM_CPPFLAGS@
CXXFLAGS=-g -O0 -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -std=c++11 @AM_CXXFLAGS@
LDFLAGS=@AM_LDFLAGS@
LIBS=$(top_srcdir)/deps/jellyfish-2.2.0/libjellyfish-2.0.la @BOOST_UNIT_TEST_FRAMEWORK_LIB@ @AM_LIBS@

TESTS_ENVIRONMENT=BOOST_TEST_REPORT_LEVEL=detailed BOOST_TEST_LOG_LEVEL=all

TESTS = check_jellyfish check_comp check_sect

check_PROGRAMS = check_jellyfish check_comp check_sect

check_jellyfish_SOURCES =   ../src/jellyfish_helper.cc \
			    check_jellyfish.cc

check_comp_SOURCES =	../src/jellyfish_helper.cc \
			../src/comp.cc \
			check_comp.cc
	
check_sect_SOURCES =	../src/jellyfish_helper.cc \
			../src/sect.cc \
			check_sect.cc
		    
