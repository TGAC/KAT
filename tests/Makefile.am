AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
CXXFLAGS=-g -O0 -fwrapv -Wall -Wextra -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-missing-field-initializers -Wno-sign-compare -std=c++11 @AM_CXXFLAGS@

CPPFLAGS=-isystem $(srcdir)/../deps/seqan-library-2.0.0/include -isystem $(srcdir)/../deps/jellyfish-2.2.0/include -isystem $(srcdir)/../include @AM_CPPFLAGS@

LDFLAGS=@AM_LDFLAGS@


TESTS = check_unit_tests

check_PROGRAMS = check_unit_tests

noinst_HEADERS =    gtest/gtest.h \
		    gtest/src/gtest-all.cc \
		    gtest/src/gtest_main.cc
	

check_unit_tests_SOURCES =	check_jellyfish.cc \
				check_spectra_helper.cc \
				check_main.cc
	
check_unit_tests_LDADD =    ../deps/jellyfish-2.2.0/.libs/libkat_jellyfish-2.0.la \
			    libgtest.la \
			    ../libkat-2.1.la \
			    @AM_LIBS@
	
include gtest.mk
