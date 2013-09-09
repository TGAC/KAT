AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

AM_LDFLAGS = -all-static -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -I$(top_srcdir)/src/inc

bin_PROGRAMS = kat


kat_SOURCES = \
                inc/gnuplot/gnuplot_i.cc \
		inc/gnuplot/gnuplot_i.hpp \
		inc/jellyfish/fstream_default.hpp \
		inc/jellyfish/jellyfish_helper.hpp \
		inc/matrix/sparse_matrix.hpp \
		inc/matrix/threaded_sparse_matrix.hpp \
		inc/matrix/matrix_metadata_extractor.hpp \
                inc/matrix/matrix_metadata_extractor.cc \
                sect/sect_main.cc \
		sect/sect.hpp \
		sect/sect_args.hpp \
		sect/sect_main.hpp \
                comp/comp_main.cc \
		comp/comp.hpp \
		comp/comp_args.hpp \
		comp/comp_main.hpp \
                gcp/gcp_main.cc \
		gcp/gcp.hpp \
		gcp/gcp_args.hpp \
		gcp/gcp_main.hpp \
		histo/histo_main.cc \
		histo/histogram.hpp \
		histo/histo_args.hpp \
		histo/histo_main.hpp \
                plot/plot_main.cc \
		plot/plot_args.hpp \
		plot/plot_main.hpp \
                plot/flame/flame_plot_main.cc \
		plot/flame/flame_plot_args.hpp \
		plot/flame/flame_plot_main.hpp \
                plot/asm/asm_plot_main.cc \
		plot/asm/asm_plot_args.hpp \
		plot/asm/asm_plot_main.hpp \
                plot/sect/sect_plot_main.cc \
		plot/sect/sect_plot_args.hpp \
		plot/sect/sect_plot_main.hpp \
                kat.cc \
		kat_args.hpp

