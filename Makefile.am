
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		INSTALL \
		ChangeLog \
		tests/data/EcoliK12.fasta \
		tests/data/ecoli.header.jf27 \
		tests/data/kat.hist \
		tests/data/sect_length_test.fa \
		tests/data/sect_test.fa


# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh clean.sh

# Scripts to install
dist_bin_SCRIPTS =	scripts/kat_distanalysis.py \
			scripts/kat_plot_misc.py \
			scripts/kat_plot_colormaps.py \
			scripts/kat_plot_density.py \
			scripts/kat_plot_profile.py \
			scripts/kat_plot_spectra-cn.py \
			scripts/kat_plot_spectra-hist.py \
			scripts/kat_plot_spectra-mx.py

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 src tests

AUTOMAKE_OPTIONS = subdir-objects
