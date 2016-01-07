
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Blank man target (used in doc)
man: 
	cd doc && $(MAKE) $(AM_MAKEFLAGS) man
.PHONY: man
	
# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		ChangeLog \
		tests/data/EcoliK12.fasta \
		tests/data/ecoli.header.jf27 \
		tests/data/kat.hist \
		tests/data/sect_length_test.fa \
		tests/data/sect_test.fa \
		tests/data/ecoli_r1.1K.fastq \
		tests/data/ecoli_r2.1K.fastq \
		tests/test_comp.sh \
		tests/test_gcp.sh \
		tests/test_hist.sh \
		tests/test_sect.sh
	

# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh antigen.sh

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
make_dirs=deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 lib src tests

if test @MAKE_DOCS@ = 1; then
    make_dirs +=  doc

SUBDIRS = $(make_dirs)

AUTOMAKE_OPTIONS = subdir-objects