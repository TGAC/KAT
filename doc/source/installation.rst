.. _installation:

Installation
============

From brew
~~~~~~~~~

If you have brew installed on your system you should be able to install a recent version of KAT by simply typing:

``brew install brewsci/bio/kat``

Many thanks to @sjackman for this one!

From bioconda
~~~~~~~~~~~~~

If you use bioconda you can install KAT using :

``bioconda install kat``


From source
~~~~~~~~~~~

If you wish to install KAT from source, because you don't have brew installed, or wish to ensure you have the latest version, first ensure these dependencies are installed and configured on your system:

  - **GCC** V4.8+
  - **make**
  - **autoconf** V2.53+
  - **automake** V1.11+
  - **libtool** V2.4.2+
  - **pthreads** (probably already installed)
  - **zlib**
  - **Python** V3.5+ with the *tabulate*, *scipy*, *numpy* and *matplotlib* packages and C API installed.  This is optional but highly recommended, without python KAT functionality is limited: no plots, no distribution analysis, and no documentation.
  - **Sphinx-doc** V1.3+ (Optional: only required for building the documentation.

NOTE ON INSTALLING PYTHON: Many system python installations do not come with the C API immediately available, which prevents KAT from embedding python code.  We typically would recommend installing anaconda3 as this would include the latest version of python, all required python packages as well as the C API.  If you are running a debian system and the C libraries are not available by default and you wish to use the system python installation the you can install them using: ``sudo apt-get install python-dev``. Also if you are using a python installation outside your system directory, please make sure you have your PATH and LD_LIBRARY_PATH (or LD_RUN_PATH) environment variables set appropriately.

Then proceed with the following steps:

  - Clone the git repository (For ssh: ``git clone git@github.com:TGAC/KAT.git``; or for https: ``git clone https://github.com/TGAC/KAT.git``), into a directory on your machine.
  - Change directory into KAT project: ``cd KAT``
  - Build boost (this may take some time): ``./build_boost.sh``
  - Setup the KAT configuration scripts by typing: ``./autogen.sh``.
  - Generate makefiles and confirm dependencies: ``./configure``. The configure script can take several options as arguments.  One commonly modified option is ``--prefix``, which will install KAT to a custom directory.  By default this is ``/usr/local``, so the KAT executable would be found at ``/usr/local/bin`` by default. Python functionality can be disabled using ``--disable-pykat``.  Type ``./configure --help`` for full list of options.  Please check the output to ensure the configuration is setup as you'd expect.
  - Compile software: ``make``.  You can leverage extra cores duing the compilation process using the ``-j <#cores>`` option.  Also you can see all command lines used to build the software by setting ``V=1``.
  - Run tests (optional) ``make check``.  (The ``-j`` and ``V=1`` options described above are also supported here.)
  - Install: ``make install``.  If you've not provided a specific installation directory, you will likely need to prefix this command with ``sudo`` in order to provide the permissions required to install to ``/usr/local``.

If sphinx is installed and detected on your system then html documentation and man
pages are automatically built during the build process.  If it is not detected then this step is skipped.  Should you wish to create a PDF version of the manual you can do so by entering the ``doc`` directory and typing ``make pdf``, this is not executed by default.

NOTE: if KAT is failing at the ``./autogen.sh`` step you will likely need to install autotools.  The following command should do this on MacOS: ``brew install autoconf automake libtool``.  On a debian system this can be done with: ``sudo apt-get install autoconf automake libtool``.

*Python scripts*

KAT will install some python scripts to your ``<prefix>/bin`` directory.  If you selected a custom location for prefix and wish to access these scripts directly, then it may be necessary to modify your $PYTHONPATH environment variable. Ensure that ``<prefix>/lib/python<version>/site-packages``, is on your PYTHONPATH, where <version> represents the python version to used when installing KAT e.g. ``/home/me/kat/lib/python3.6/site-packages``.  Alternatively, you could install the kat python package into a python environment by changing into the ``scripts`` directory and typing ``python setup.py install``.
