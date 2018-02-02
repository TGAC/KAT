AUTOMAKE_OPTIONS = subdir-objects

install-exec-hook:
	export PYTHONPATH=$(prefix)/lib/python$(PYTHON_VERSION)/site-packages && \
	mkdir -p $(prefix)/lib/python$(PYTHON_VERSION)/site-packages/ && \
	python3 setup.py install --prefix=$(prefix) --record install_files.txt

uninstall-hook:
	cat install_files.txt | xargs rm -rf

EXTRA_DIST = \
	kat/kat_distanalysis.py \
	kat/kat_plot_colormaps.py \
	kat/kat_plot_misc.py \
	kat/kat_plot_profile.py \
	kat/kat_plot_spectra_cn.py \
	kat/kat_plot_spectra_hist.py \
	kat/kat_plot_spectra_mx.py \
	kat/kat_plot_cold.py

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf dist build kat.egg-info
