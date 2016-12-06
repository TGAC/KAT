AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

EXTRA_DIST = \
	data/EcoliK12.fasta \
	data/ecoli.header.jf27 \
	data/kat.hist \
	data/sect_length_test.fa \
	data/sect_test.fa \
	data/ecoli_r1.1K.fastq \
	data/ecoli_r2.1K.fastq \
	data/unknown.dat \
	test_comp.sh \
	test_gcp.sh \
	test_hist.sh \
	test_sect.sh
	
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf temp

# additional include paths necessary to compile the C++ library
AM_CXXFLAGS=-g -O0 -fwrapv -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -Wno-missing-field-initializers -Wno-sign-compare -std=c++11 @AM_CXXFLAGS@
AM_CPPFLAGS =  -isystem $(top_srcdir)/deps/seqan-library-2.0.0/include \
	    -isystem $(top_srcdir)/deps/jellyfish-2.2.0/include \
	    -isystem $(top_srcdir)/lib/include \
	    -DDATADIR=\"$(srcdir)/data\" \
	    @AM_CPPFLAGS@

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_SH_LOG_FLAGS =

TESTS = check_unit_tests test_hist.sh test_gcp.sh test_sect.sh test_comp.sh

check_PROGRAMS = check_unit_tests

noinst_HEADERS = \
	gtest/gtest.h \
	gtest/src/gtest-all.cc \
	gtest/src/gtest_main.cc
	

check_unit_tests_SOURCES = \
	check_jellyfish.cc \
	check_spectra_helper.cc \
	check_compcounters.cc \
	check_main.cc

check_unit_tests_LDFLAGS = \
	-L../lib \
	-L../deps/jellyfish-2.2.0 \
	@AM_LDFLAGS@ \
	@LDFLAGS@
	
check_unit_tests_LDADD = \
	libgtest.la \
	-lkat \
	-lkat_jellyfish \
	@AM_LIBS@
	
include gtest.mk
