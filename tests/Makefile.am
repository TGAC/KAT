AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf temp

# additional include paths necessary to compile the C++ library
CXXFLAGS=-g -O0 -fwrapv -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-missing-field-initializers -Wno-sign-compare -std=c++11 @AM_CXXFLAGS@
CPPFLAGS =  -isystem $(srcdir)/../deps/seqan-library-2.0.0/include \
	    -isystem $(srcdir)/../deps/jellyfish-2.2.0/include \
	    -isystem $(srcdir)/../include \
	    -DDATADIR=\"$(srcdir)/data\" \
	    @AM_CPPFLAGS@

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_SH_LOG_FLAGS =

TESTS = check_unit_tests test_hist.sh test_gcp.sh test_sect.sh test_comp.sh

check_PROGRAMS = check_unit_tests

noinst_HEADERS =    gtest/gtest.h \
		    gtest/src/gtest-all.cc \
		    gtest/src/gtest_main.cc
	

check_unit_tests_SOURCES =	check_jellyfish.cc \
				check_spectra_helper.cc \
				check_compcounters.cc \
				check_main.cc

check_unit_tests_LDFLAGS = @AM_LDFLAGS@
	
check_unit_tests_LDADD =    ../deps/jellyfish-2.2.0/.libs/libkat_jellyfish-2.0.la \
			    libgtest.la \
			    ../lib/libkat-2.1.la \
			    @AM_LIBS@
	
include gtest.mk
