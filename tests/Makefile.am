AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
CXXFLAGS=-g -O0 -fwrapv -Wall -Wextra -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-missing-field-initializers -Wno-sign-compare -std=c++11 @AM_CXXFLAGS@

CPPFLAGS=-DCPLUSPLUS -isystem $(srcdir)/../deps/seqan-library-2.0.0/include -isystem $(srcdir)/../deps/jellyfish-2.2.0/include -I$(srcdir)/../src/inc -I$(srcdir)/../src @AM_CPPFLAGS@

LDFLAGS=@AM_LDFLAGS@


TESTS = check_unit_tests

check_PROGRAMS = check_unit_tests

noinst_HEADERS =    gtest/gtest.h \
		    gtest/src/gtest-all.cc \
		    gtest/src/gtest_main.cc
	

check_unit_tests_SOURCES =	../src/inc/gnuplot/gnuplot_i.cc \
				../src/inc/matrix/matrix_metadata_extractor.cc \
				../src/input_handler.cc \
				../src/jellyfish_helper.cc \
				../src/plot_density.cc \
				../src/plot_spectra_cn.cc \
				../src/plot_spectra_mx.cc \
				../src/plot_spectra_hist.cc \
				../src/plot.cc \
				../src/comp.cc \
				../src/sect.cc \
				check_jellyfish.cc \
				check_spectra_helper.cc \
				check_sect.cc \
				check_main.cc
	
check_unit_tests_LDADD =    $(builddir)/../deps/jellyfish-2.2.0/.libs/libkat_jellyfish-2.0.a \
			    libgtest.la \
			    @AM_LIBS@
	
include gtest.mk
