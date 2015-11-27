
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

SEQAN_ARCHIVE = deps/seqan-library-2.0.0.tar.gz

dist_pkgdata_DATA = $(SEQAN_ARCHIVE)

install-data-hook:
	cd $(DESTDIR)/deps; \
	cat $(SEQAN_ARCHIVE) | gunzip | $(TAR) -H gnu --list > uninstall_manifest.txt; \
	cat $(SEQAN_ARCHIVE) | gunzip | $(TAR) -H gnu --no-same-owner --extract; \
	rm --force $(SEQAN_ARCHIVE); \
	cat uninstall_manifest.txt | sed --expression='s/^\|$$/"/g' | xargs chmod a=rX,u+w

uninstall-local:
	cd $(DESTDIR)/deps; \
	cat uninstall_manifest.txt | sed --expression='s/ /\\ /g' | xargs rm --force; \
	rm --force uninstall_manifest.txt

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		INSTALL \
		ChangeLog \
		deps/seqan-library-2.0.0


# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh clean.sh

# Scripts to install
dist_bin_SCRIPTS =	scripts/dist_analysis.py \
			scripts/colormaps.py \
			scripts/density.py \
			scripts/findpeaks.py \
			scripts/header.py \
			scripts/profile.py \
			scripts/spectra-cn.py \
			scripts/spectra-hist.py \
			scripts/spectra-mx.py

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 src tests

AUTOMAKE_OPTIONS = subdir-objects
