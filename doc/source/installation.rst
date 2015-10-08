.. _installation:

Installation
============

KAT is primarily a C++ application with a few python scripts.  We use the 
GNU build system ``Autotools`` to assist with package management and to make the 
software portable across UNIX type operating systems.  Installation of KAT
therefore follows a similar incantation to other autotools based projects::

  ```./configure; make; sudo make install;```

However, there are a few caveats.  If you cloned the software directly from the 
git repository you must first run "./autogen.sh" to create the configure and make 
files for your project.  If you downloaded a source code distribution tarball those
scripts are already present so you can skip this step.

External Dependencies
---------------------

KAT depends on some external software, specifically boost and gnuplot, although you
will also need a C++11 capable compiler to compile the code.  Please 
make sure these programs are correctly configured and installed 
on your system prior to building KAT.  Consult the each program's installation
guide separately for instructions on how to do this.  Should you install these dependencies
into non-standard locations you can direct KAT to them by using the following
options when running the configure script.

  - ```--with-boost``` - for specifiying a custom boost installation directory
  - ```--with-zlib``` - for specifying a custom zlib installation directory

We prefer static linking, which is forced in the configuation script.  Therefore
you do not need to modify your LD_LIBRARY_PATH to recognise boost at runtime. 

In addition, should you wish to run the dist_analysis.py script you will need python with the
numpy and scipy libraries available.


Internal Dependencies
---------------------

KAT contains jellyfish and SeqAn in the source tree.  The user does
not need to do anything specical to handle these dependencies as they are automatically
built and managed inside KAT.  However, it is important to note that KAT
will create the jellyfish executable in it's bin directory after installation, which
may conflict with your own jellyfish executable if it was already installed on your
PATH.  If you do not want KAT to potentially override or conflict with an 
existing jellyfish installation you might want to consider installing KAT
to a custom location.  You can do this with the ```--prefix``` option when 
running the configure script.  We might revisit this in the future to remove
this potential issue.


Compilation and Installation
----------------------------

First change into the KAT root directory and run ```./configure```, providing
any options you feel are appropriate.  By default the installation directory is "/usr/local", 
so the KAT executable would be found at "/usr/local/bin" by default.  If you
want to change this use the ``--prefix`` option as previously described.  For a full
list of configuration options type ```./configure --help```.

Next compile the software.  This can be done by typing ```make```.  The compiles
all internal dependencies and KAT itself.

To check the code compiled correct and is operating as expected you can optionally
type  ```make check``` to runs some tests.  This includes unit tests for jellyfish 
which are embedded in the KAT source tree.  To run only KAT
unit tests go into the ``tests`` subdirectory and run ```make check``` there.

Finally to install the compiled code to the specified (or default) installation
directory type ```make install```.
