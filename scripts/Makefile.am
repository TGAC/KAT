AUTOMAKE_OPTIONS = subdir-objects

if MAKE_PYKAT_NOINSTALL

install-exec-hook:
	mkdir -p $(DESTDIR)$(datadir)/scripts/kat/plot
	cp -r kat/*.py $(DESTDIR)$(datadir)/scripts/kat
	cp -r kat/plot/*.py $(DESTDIR)$(datadir)/scripts/kat/plot

uninstall-hook:
	cat install_files.txt | xargs rm -rf
	rm -rf $(DESTDIR)$(datadir)/scripts

endif

if MAKE_PYKAT_INSTALL

install-exec-hook:
	$(PYTHON_INT_PATH) setup.py install --record install_files.txt
	mkdir -p $(DESTDIR)$(datadir)/scripts/kat/plot
	cp -r kat/*.py $(DESTDIR)$(datadir)/scripts/kat
	cp -r kat/plot/*.py $(DESTDIR)$(datadir)/scripts/kat/plot

uninstall-hook:
	cat install_files.txt | xargs rm -rf
	rm -rf $(DESTDIR)$(datadir)/scripts

endif


if MAKE_PYKAT_INSTALL_CUSTOM

install-exec-hook:
	mkdir -p $(DESTDIR)$(prefix)/lib/python$(PYTHON_VERSION)/site-packages 
	export PYTHONPATH=$(DESTDIR)$(prefix)/lib/python$(PYTHON_VERSION)/site-packages && \
	$(PYTHON_INT_PATH) setup.py install --prefix=$(DESTDIR)$(prefix) --record install_files.txt
	mkdir -p $(DESTDIR)$(datadir)/scripts/kat/plot
	cp -r kat/*.py $(DESTDIR)$(datadir)/scripts/kat
	cp -r kat/plot/*.py $(DESTDIR)$(datadir)/scripts/kat/plot

uninstall-hook:
	cat install_files.txt | xargs rm -rf
	rm -rf $(DESTDIR)$(datadir)/scripts

endif

 


EXTRA_DIST = \
	kat/distanalysis.py \
	kat/plot/__init__.py \
	kat/plot/colormaps.py \
	kat/plot/profile.py \
	kat/plot/spectra_cn.py \
	kat/plot/spectra_hist.py \
	kat/plot/spectra_mx.py \
	kat/plot/cold.py \
	kat/peak.py \
	kat/spectra.py \
	kat/__init__.py

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf dist build kat.egg-info
