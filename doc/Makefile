# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = build
MKDIR_P	      = mkdir -p
INSTALL	      = /usr/bin/install -c -m 644

# Set autoconf variables
prefix	      = /tgac/software/testing/kat/dev/x86_64
PACKAGE_TARNAME = kat
top_srcdir    = ..
srcdir	      = .
docdir	      = ${datarootdir}/doc/${PACKAGE_TARNAME}
mandir	      = ${datarootdir}/man
htmldir	      = ${docdir}
pdfdir        = ${docdir}
infodir       = ${datarootdir}/info
top_builddir  = ..
datarootdir   = ${prefix}/share
datadir       = ${datarootdir}

# User-friendly check for sphinx-build
#if test $(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?) = 1; then
#$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
#fi

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) -c source $(top_srcdir)/doc/source
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) -c source $(top_srcdir)/doc/source

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext


DPS=$(top_srcdir)/doc/source
DPI=$(DPS)/images
DISTFILES = $(DPS)/index.rst \
	    $(DPS)/installation.rst \
	    $(DPS)/kmer.rst \
	    $(DPS)/using.rst \
	    $(DPS)/walkthrough.rst \
	    $(DPS)/faq.rst \
	    $(DPI)/ccoli_comp.png \
	    $(DPI)/ccoli_gcp.png \
	    $(DPI)/ccoli_hist.png \
	    $(DPI)/contaminant_all.png \
	    $(DPI)/contaminant_MP.png \
	    $(DPI)/contaminant_ope1.png \
	    $(DPI)/contaminant_ope2.png \
	    $(DPI)/contaminant_PE.png \
	    $(DPI)/dist_analysis.png \
	    $(DPI)/gc_bias_a.png \
	    $(DPI)/gc_bias_b.png \
	    $(DPI)/gc_bias_c.png \
	    $(DPI)/gc_bias_d.png \
	    $(DPI)/heterozygous_perfect_1.png \
	    $(DPI)/heterozygous_perfect_2.png \
	    $(DPI)/heterozygous_real.png \
	    $(DPI)/kat_logo.png \
	    $(DPI)/kmer_spectra1.png \
	    $(DPI)/kmer_spectra_breakdown.png \
	    $(DPI)/pe_v_asm_clean.png \
	    $(DPI)/pe_v_asm_wrong.png \
	    $(DPI)/pe_v_mp_after_density.png \
	    $(DPI)/pe_v_mp_after_shared.png \
	    $(DPI)/pe_v_mp_before_density.png \
	    $(DPI)/pe_v_mp_before_shared.png \
	    $(DPI)/pe_v_pe_1_density.png \
	    $(DPI)/pe_v_pe_1_shared.png \
	    $(DPI)/pe_v_pe_2_density.png \
	    $(DPI)/pe_v_pe_2_shared_1.png \
	    $(DPI)/pe_v_pe_2_shared_2.png \
	    $(DPI)/pe_v_rna_bar.png \
	    $(DPI)/pe_v_rna_shared.png \
	    $(DPI)/profile.png \
	    $(DPI)/real_r1_v_r2.png \
	    $(DPI)/real_r1_v_r2_shared.png \
	    $(DPI)/simulated_r1_v_r2.png

	

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"

clean:
	rm -rf $(BUILDDIR)/*

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) $(BUILDDIR)/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/KAT.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/KAT.qhc"

devhelp:
	$(SPHINXBUILD) -b devhelp $(ALLSPHINXOPTS) $(BUILDDIR)/devhelp
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/KAT"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/KAT"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdfja:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through platex and dvipdfmx..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	$(SPHINXBUILD) -b text $(ALLSPHINXOPTS) $(BUILDDIR)/text
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) -b man $(ALLSPHINXOPTS) $(BUILDDIR)/man
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) -b texinfo $(ALLSPHINXOPTS) $(BUILDDIR)/texinfo
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

xml:
	$(SPHINXBUILD) -b xml $(ALLSPHINXOPTS) $(BUILDDIR)/xml
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml:
	$(SPHINXBUILD) -b pseudoxml $(ALLSPHINXOPTS) $(BUILDDIR)/pseudoxml
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."

# Automake targets that we'll ignore
EMPTY_AUTOMAKE_TARGETS = check installcheck tags ctags dvi install-exec install-dvi install-ps install-info
.PHONY: $(EMPTY_AUTOMAKE_TARGETS)
$(EMPTY_AUTOMAKE_TARGETS):

# Automake targets that we'll use (that aren't already described)
.PHONY: all man pdf mostly-clean distclean maintainer-clean install install-data install-man install-html install-pdf installdirs uninstall distdirs

# Hook in conventional make targets for building docs
all: html man
        
# Redirect pdf to latexpdf
pdf: latexpdf

# Redirect all cleaning targets to clean
mostly-clean: clean
distclean: clean
maintainer-clean: clean

# Install and uninstall overrides
install: install-data

install-data: install-man install-pdf install-html

# Fail on error
install-man:
	test -s build && test -s build/man && $(INSTALL) build/man/kat.1 $(DESTDIR)$(mandir)/man1/kat.1

# Skip any errors here... we don't require the PDF to be generated
install-pdf:
	-test -s build && test -s build/latex && test -s build/latex/kat.pdf && $(INSTALL) build/latex/kat.pdf $(DESTDIR)$(pdfdir)/kat.pdf;

install-html:
	test -s build && test -s build/html && $(MKDIR_P) $(DESTDIR)$(htmldir) && cp -rf build/html $(DESTDIR)$(htmldir) && chmod 755 -R $(DESTDIR)$(htmldir)

# Make sure all installation directories (e.g. $(bindir))
# actually exist by making them if necessary.
installdirs: mkinstalldirs
	$(srcdir)/mkinstalldirs \
		$(DESTDIR)$(pdfdir) \
		$(DESTDIR)$(htmldir) \
		$(DESTDIR)$(mandir)/man1

uninstall:
	rm -Rf "$(DESTDIR)$(mandir)/man1/kat.1" "$(DESTDIR)$(pdfdir)/kat.pdf" "$(DESTDIR)$(htmldir)"

	
distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
