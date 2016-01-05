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


External Dependencies
---------------------

KAT depends on some external software, in particular boost.  If you don' wish to
install the full suite of boost libraries KAT only uses the following:

 - system
 - filesystem
 - program_options
 - chrono
 - timer. 

You will also need "make", "libtool", "zlib"
and a C++11 capable compiler such as "GCC V4.8+" to 
compile the code.  

Please make sure all these programs and libraries are correctly configured and installed 
on your system prior to building KAT.  Consult the each program's installation
guide separately for instructions on how to do this.  For boost and zlib KAT
offers you the ability to build these from non-standard locations without setting
environment variables by using the following options when running the configure script.

  - "--with-boost=<path_to_boost>"  for specifying a custom boost installation directory
  - "--with-zlib=<path_to_zlib>"  for specifying a custom zlib installation directory

However, please ensure that all libraries are available on the LD_LIBRARY_PATH at runtime. 

To enable plotting functionality we require either python3, with numpy, scipy and
matplotlib installed, or gnuplot.  The python plotting method is the preferred
method and will produce nicer results.  If you don't already have python3 installed
on your system we recommend installing anaconda3 as this contains everything you
need.  The type of plotting engine used will be determined when running the configure
script, which will select the first engine detected in the following order: python,
gnuplot, none.  There is currently no way to select the plotting directory from
a custom location, so the plotting system needs to be properly installed and configured
on your system: i.e. python3 or gnuplot must be available on the PATH.

For an example of how to install KAT dependencies on a typical system
take a look at our travis CI script in the root directory: ``.travis.yml``

If you have cloned the repository then you will also need a few additional dependencies installed
to generate the configuration script.  These are:
 
   - autoconf V2.53+
   - automake V1.11+


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

Finally to install the compiled code to the specified (or default) installation
directory type ``make install``.
