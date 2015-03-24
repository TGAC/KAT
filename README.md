
#KAT - The K-mer Analysis Toolkit

KAT is a suite of tools that analyse jellyfish kmer hashes.  The following tools are currently available in KAT:

   - **sect**:  SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in a fasta file using K-mers from a jellyfish hash.
   - **comp**:  K-mer comparison tool.  Creates a matrix of shared K-mers between two jellyfish hashes.
   - **gcp:**   K-mer GC Processor.  Creates a matrix of the number of K-mers found given a GC count and a K-mer count.
   - **hist**:  Create an histogram of k-mer occurrences from a jellyfish hash.  Adds metadata in output for easy plotting.
   - **plot**:  Plotting tool.  Contains several plotting tools to visualise K-mer and compare distributions. Requires gnuplot.  The following plot tools are available:

     - **density**:      Creates a density plot from a matrix created with the "comp" tool.  Typically this is used to compare two K-mer hashes produced by different NGS reads.
     - **profile**:      Creates a K-mer coverage plot for a single sequence.  Takes in fasta coverage output coverage from the "sect" tool
     - **spectra-cn**:   Creates a stacked histogram using a matrix created with the "comp" tool.  Typically this is used to compare a jellyfish hash produced from a read set to a jellyfish hash produced from an assembly. The plot shows the amount of distinct K-mers absent, as well as the copy number variation present within the assembly.
     - **spectra-hist**: Creates a K-mer spectra plot for a set of K-mer histograms produced either by jellyfish-histo or kat-histo.
     - **spectra-mx**:   Creates a K-mer spectra plot for a set of K-mer histograms that are derived from selected rows or columns in a matrix produced by the "comp".

In addition, KAT contains a python script for analysing the mathematical distributions present in the K-mer spectra in order to determine how much content is present in each peak.


##Installation:

Generic installation description can be found in the INSTALL file. Short summary: 

  - Acquire the source code.  Either download and decompress the distributable ("tar -xvf kat-<version>.tar.gz"), or clone the git repository (For ssh: ```git clone git@github.com:TGAC/KAT.git```; or for https: ```git clone https://github.com/TGAC/KAT.git```), into a directory on your machine.
  - "cd" into root directory of the installation
  - Ensure these tools are correctly installed and available on your system:
      - gcc tool chain
      - make
      - jellyfish = V1.1.10 or V1.1.11 - http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz **IMPORTANT NOTE**: Please use jellyfish V1.1, we currently do not support jellyfish 2.   We will update KAT to support newer versions of jellyfish in due course.
      - gnuplot (required for plotting at runtime, must be available on the path to use this functionality) - http://www.gnuplot.info
  - If you cloned the git repository you must first run "./autogen.sh" to create the configure and make files for your project.  Do not worry if this fails due to missing dependencies at this stage.  If you downloaded a source code distribution tarball then you can skip this step.
  - For a typical installation on a machine where you have root access type ```./configure; make; sudo make install;```

The configure script can take several options as arguments.  One commonly modified option is ```--prefix```, which will install KAT to a custom directory.  By default this is "/usr/local", so the KAT executable would be found at "/usr/local/bin" by default.  In addition, some options specific to managing KAT dependencies located in non-standard locations are:

  - ```--with-jellyfish``` - for specifying a custom jellyfish directory
  - ```--with-boost``` - for specifying a custom boost directory (boost is only required for unit testing)
  - ```--with-doxygen``` - for specifying a custom doxygen directory (doxygen is only required for generating code documention.

Type ```./configure --help``` for full details.

The Makefile for KAT can take several goals.  Full details of common make goals can be found in the INSTALL file.  Typically, the following options can optionally used by KAT:

  - ```make check``` - runs unit tests.  *NOTE*: Requires boost unit test framework to be installed and available.
  - ```make dist``` - packages the installation into a tarballed distributable.
  - ```make distcheck``` - runs some sanity tests to ensure the tarballed distributable is likely to work.

KAT also come with a python script called "dist_analysis.py", which allows the user to determine the amount of content under each peak in the K-mer spectra.  This script will automatically be moved into your selected install directory after running "make install".  However, before running this script you will need to install python and the python numpy and scipy libraries.  In future versions this script will be properly integrated into KAT as another subtool written in C++.


##Operating Instructions:

After KAT has been installed, the following tools should be available:

 - **kat** - a single executable binary file that contains a number of subtools.
 - **kat_comp_reads.sh** - a bash script demonstrating a simple pipeline to compare the K-mers in two read files
 - **dist_analysis.py** - a python script for determining the amount of content in each peak in the K-mer spectra

Running ```kat --help``` will bring up a list of available tools within kat.  To get help on any of these subtools simple type: ```kat <tool> --help```.  For example: ```kat sect --help``` will show details on how to use the sequence coverage estimator tool.

Specifically, jellyfish must be available for dynamic linking at runtime.  In addition, in order to use the plotting tools it is necessary for "gnuplot" to be available in the PATH.


##Extending KAT:

Developers can extend KAT by adding additional tools, whilst leveraging some of the shared resources that KAT and Jellyfish have made available.  In order to add an additional tool to KAT, developers will need a reasonable working knowledge of C++ programming and have GNU auto tools available on their system.  The process for adding a new subtool is as follows:

1. Create a new directory with the tools name in the "src" directory
2. Copy the template_args.hpp file into this directory and rename to whatever you wish.  Modify the template file so that it contains details of how to use your tool.  Comments have been added to the template to indicate places where you will have to add your custom code. The args template file makes use of getopt.h so developers familiar with this library should have no issues here.  For those unfamiliar with this library, please read the getopt documentation: http://www.gnu.org/software/libc/manual/html_node/Getopt.html
3. Copy the template_main.cc and _main.hpp files into the new directory and write whatever code is necessary for your tool.
4. Add an include and extend the validMode method in "src/kat.cc" so that your tool is recognised.  Also add your tool to the longDescription method.
5. Update the "src/kat_args.hpp" to extend the KAT help messages.
6. Update the Makefile.am file to include your _main.cc file.
7. Run ```aclocal; autoconf; automake``` to generate the actual configure script and initial Makefiles.
8. Run ```./configure```, with any appropriate options, to make the final Makefiles.
9. Run ```make``` to compile the new version of KAT with your tools included.  The KAT binary will be available in the "./bin" directory.
10. Run ```sudo make install``` to install the software.

See INSTALL file for more details on configuring steps 8-10.

There are some shared resources available which might aid the generation of a subtool.  It is worth browsing the ./src/inc directory to see what is available.  There are libraries for:

- Easing generation of gnuplot commands.  Code was taken and modified from: http://ndevilla.free.fr/gnuplot/
- "jellyfish_helper.hpp" provides some convienient functionality for loading an managing jellyfish hashes from a simple file path.
- Sparse Matrix implementation.  In order to avoid loading heavy dependencies such as boost a simple sparse matrix implementation has been added to store matricies in a relatively memory efficient way.  The code was originally taken from: http://www.cplusplus.com/forum/general/8352/ and modified for use in KAT.  If more functionality is required than is available here, either extend this class or use a dedicated matrix library.
- string and file utils.  Some shortcuts to commonly used string and file operations that would otherwise only be available by adding another library as a dependency to this project.

If you think your subtool is useful and want it available in the official KAT release then please contact daniel.mapleson@tgac.ac.uk or bernardo.clavijo@tgac.ac.uk for discussions on how to harmonise the code.  The job will be easier if you maintain a branch from a clone or fork of the KAT repository on github.


##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

Bernardo Clavijo
Daniel Mapleson
Darren Heavens
Sarah Ayling
Mario Caccamo

See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
