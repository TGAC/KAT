AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @AM_CPPFLAGS@
AM_CXXFLAGS=-Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
LIBS=@MAIN_LIBS@

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
                inc/str_utils.hpp \
                inc/file_utils.hpp \
                inc/common_args.hpp \
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
		hist/hist_main.cc \
		hist/histogram.hpp \
		hist/hist_args.hpp \
		hist/hist_main.hpp \
                plot/plot_main.cc \
		plot/plot_args.hpp \
		plot/plot_main.hpp \
                plot/common_plot_args.hpp \
                plot/density/density_plot_main.cc \
		plot/density/density_plot_args.hpp \
		plot/density/density_plot_main.hpp \
                plot/spectra-cn/spectra_cn_plot_main.cc \
		plot/spectra-cn/spectra_cn_plot_args.hpp \
		plot/spectra-cn/spectra_cn_plot_main.hpp \
                plot/profile/profile_plot_main.cc \
                plot/profile/profile_plot_main.hpp \
		plot/profile/profile_plot_args.hpp \
		plot/spectra-hist/spectra_hist_plot_main.cc \
		plot/spectra-hist/spectra_hist_plot_args.hpp \
		plot/spectra-hist/spectra_hist_plot_main.hpp \
		plot/spectra-mx/spectra_mx_plot_main.cc \
		plot/spectra-mx/spectra_mx_plot_args.hpp \
		plot/spectra-mx/spectra_mx_plot_main.hpp \
                kat.cc \
		kat_args.hpp

