.. _installation:

Installation
============

To get a quick summary of how to install KAT please consult the README.md file
in the root directory of the project.  This section of the documentation discusses
some more details.

KAT is primarily a C++ application with a few python scripts.  We use the 
GNU build system "Autotools" (autoconf + automake) to assist with package management and to make the 
software portable across UNIX type operating systems.  Installation of KAT
therefore follows a similar incantation to other autotools based projects::

  ./configure; make; sudo make install;

Should you wish to run tests then you can do this by typing: ``make check``.

If you cloned the software directly from the git repository you must first run 
``./autogen.sh`` to create the configure and make 
files for your project.  If you downloaded a source code distribution tarball those
scripts are already present so you can skip this step.  We also supply a script called
``antigen.sh`` which cleans the KAT configuration to a similar state as if the
repository was cloned.  After running ``antigen.sh`` (may require sudo permissions 
if you installed content or run ``make dist``), ``autogen.sh`` must be
run in order to create the ``configure`` and ``Makefile`` files in order to 
build the problem.

The configure script contains all the usual options you might expect in particular
it has the ``--prefix``, which will install KAT to a custom directory.  By default 
this is "/usr/local", so the KAT executable would be found at "/usr/local/bin" by 
default.  In addition, the make file should support all the usual targets and options
such as ``-j <core>``, to increase the number of cores used to compile the program.

If sphinx is installed and detected on your system then html documentation and man 
pages are automatically built during the build process.  If it is not detected then this step
is skipped.  Should you wish to create a PDF version of the manual you can do so
by typing ```make pdf```, this is not executed by default.


External Dependencies
---------------------

KAT depends on some external software:
 * Make
 * libtool
 * zlib
 * C++11 compiler (e.g. GCC V4.8+)
 * Boost
 * Optional but recommended - Python3 with matplotlib, scipy, numpy and sphinx (anaconda3)
 * Optional - Gnuplot
 * Repo building - automake / autoconf

Please make sure all the mandatory programs and libraries are correctly configured and installed 
on your system prior to building KAT.  The rest of this section describes how to 
install these dependencies in more detail. However, for a quick example of how to install KAT dependencies on a typical clean system
take a look at our travis CI script in the root directory: ``.travis.yml``.


Boost
~~~~~

If you don't wish to install the full suite of boost libraries KAT only uses the following:

 - system
 - filesystem
 - program_options
 - chrono
 - timer. 

KAT offers the use the ability to use boost from non-standard locations without setting
environment variables via the following options when running the configure script::

  - "--with-boost=<path_to_boost>"  for specifying a custom boost installation directory
  
However, if you choose to do this please ensure that the boost libraries are available 
on the LD_LIBRARY_PATH at runtime. 

Zlib
~~~~

Most users will already have this installed and configured on their system, however
if you don't then please install it.  Should you choose to install it in a custom directory
then you can use this option in the configure script::

  - "--with-zlib=<path_to_zlib>"  for specifying a custom zlib installation directory

Again please be sure to have the LD_LIBRARY_PATH set if you choose to do this.


Plotting
~~~~~~~~

To enable plotting functionality we require either python3, with numpy, scipy and
matplotlib installed, or gnuplot.  The python installation must come with the python
shared library, on debian systems you can install this with "sudo apt-get install python3-dev".
Although, if you don't already have python3 installed
on your system we recommend installing anaconda3 as this contains everything you
need.  The type of plotting engine used will be determined when running the configure
script, which will select the first engine detected in the following order: python,
gnuplot, none.  The python plotting method is the preferred
method and will produce nicer results.  There is currently no way to select the plotting directory from
a custom location, so the plotting system needs to be properly installed and configured
on your system: i.e. python3 or gnuplot must be available on the PATH.


Documentation
~~~~~~~~~~~~~

Should you wish to build the documentation, you will need python3 with the sphinx
package installed and the sphinx-build executable on the PATH.  If you have already
installed anaconda3 then you will have this already.


Building from git repo
~~~~~~~~~~~~~~~~~~~~~~

If you have cloned the repository then you will also need a few additional dependencies installed
to generate the configuration script.  These are::
 
 * autoconf V2.53+
 * automake V1.11+



Internal Dependencies
---------------------

KAT contains jellyfish and SeqAn in the source tree.  The user does
not need to do anything special to handle these dependencies as they are automatically
built and managed inside KAT.  However, it is important to note that KAT
will create the jellyfish executable in it's bin directory after installation, which
may conflict with your own jellyfish executable if it was already installed on your
PATH.  If you do not want KAT to potentially override or conflict with an 
existing jellyfish installation you might want to consider installing KAT
to a custom location.  You can do this with the ``--prefix`` option when 
running the configure script.  We might revisit this in the future to remove
this potential issue.


Compilation and Installation
----------------------------

First change into the KAT root directory and run ``./configure``, providing
any options you feel are appropriate.  By default the installation directory is ``/usr/local``, 
so the KAT executable would be found at ``/usr/local/bin`` by default.  If you
want to change this use the ``--prefix`` option as previously described.  For a full
list of configuration options type ``./configure --help``.

Next compile the software.  This can be done by typing ``make``.  The compiles
all internal dependencies and KAT itself.

To check the code compiled correct and is operating as expected you can optionally
type  ``make check`` to runs some tests.  This includes unit tests for jellyfish 
which are embedded in the KAT source tree.  To run only KAT
unit tests go into the ``tests`` subdirectory and run ``make check`` there.

Should you have sphinx installed and wish to create a PDF copy of the manual, you
can do so by typing ``make pdf``.

Finally to install the compiled code to the specified (or default) installation
directory type ``make install``.
