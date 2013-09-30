HEADERS += \
    src/comp/comp_main.hpp \
    src/comp/comp_args.hpp \
    src/comp/comp.hpp \
    src/inc/jellyfish/jellyfish_helper.hpp \
    src/inc/jellyfish/fstream_default.hpp \
    src/sect/sect_main.hpp \
    src/sect/sect_args.hpp \
    src/sect/sect.hpp \
    src/plot/plot_main.hpp \
    src/plot/plot_args.hpp \
    src/inc/gnuplot/gnuplot_i.hpp \
    src/plot/contamination/contamination_plot_main.hpp \
    src/plot/contamination/contamination_plot_args.hpp \
    src/plot/profile/profile_plot_main.hpp \
    src/plot/profile/profile_plot_args.hpp \
    src/kat_args.hpp \
    templates/template_main.hpp \
    templates/template_args.hpp \
    src/inc/matrix/threaded_sparse_matrix.hpp \
    src/inc/matrix/sparse_matrix.hpp \
    src/gcp/gcp_main.hpp \
    src/gcp/gcp_args.hpp \
    src/gcp/gcp.hpp \
    src/inc/matrix/matrix_metadata_extractor.hpp \
    src/plot/spectra/spectra_plot_main.hpp \
    src/plot/spectra/spectra_plot_args.hpp \
    src/inc/str_utils.hpp \
    src/inc/common_args.hpp \
    src/plot/common_plot_args.hpp \
    src/plot/density/density_plot_main.hpp \
    src/plot/density/density_plot_args.hpp \
    src/plot/spectra-cn/spectra_cn_plot_main.hpp \
    src/plot/spectra-cn/spectra_cn_plot_args.hpp \
    src/plot/spectra-hist/spectra_hist_plot_main.hpp \
    src/plot/spectra-hist/spectra_hist_plot_args.hpp \
    src/hist/histogram.hpp \
    src/hist/hist_main.hpp \
    src/hist/hist_args.hpp \
    src/plot/spectra-mx/spectra_mx_plot_main.hpp \
    src/plot/spectra-mx/spectra_mx_plot_args.hpp \
    src/inc/file_utils.hpp

SOURCES += \
    src/comp/comp_main.cc \
    src/sect/sect_main.cc \
    src/plot/plot_main.cc \
    src/inc/gnuplot/gnuplot_i.cc \
    src/plot/contamination/contamination_plot_main.cc \
    src/plot/profile/profile_plot_main.cc \
    src/kat.cc \
    templates/template_main.cc \
    src/gcp/gcp_main.cc \
    src/inc/matrix/matrix_metadata_extractor.cc \
    tests/check_comp.cc \
    src/plot/spectra/spectra_plot_main.cc \
    src/plot/density/density_plot_main.cc \
    src/plot/spectra-cn/spectra_cn_plot_main.cc \
    src/plot/spectra-hist/spectra_hist_plot_main.cc \
    src/hist/hist_main.cc \
    src/plot/spectra-mx/spectra_mx_plot_main.cc

OTHER_FILES += \
    src/Makefile.am \
    Makefile.am \
    tests/Makefile.am \
    configure.ac
