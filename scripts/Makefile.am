AUTOMAKE_OPTIONS = subdir-objects


if MAKE_PYKAT_INSTALL

install-exec-hook:
	export PYTHONPATH=$(prefix)/lib/python$(PYTHON_VERSION)/site-packages && \
	mkdir -p $(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	python3 setup.py install --prefix=$(prefix) --record install_files.txt

uninstall-hook:
	cat install_files.txt | xargs rm -rf

else

install-exec-hook:
	mkdir -p $(prefix)/lib/python$(PYTHON_VERSION)/local/ && \
	cp -r kat $(prefix)/lib/python$(PYTHON_VERSION)/local

uninstall-hook:
	rm -rf $(prefix)/lib/python$(PYTHON_VERSION)/local

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
