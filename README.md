
#KAT - The K-mer Analysis Toolkit

KAT is a suite of tools that analyse jellyfish hashes or sequence files (fasta or fastq) using kmer counts.  The following tools are currently available in KAT:

   - **sect**:   SEquence Coverage estimator Tool.  Estimates the coverage of each sequence in a file using K-mers from another sequence file.
   - **comp**:   K-mer comparison tool.  Creates a matrix of shared K-mers between two (or three) sequence files or hashes.
   - **gcp:**    K-mer GC Processor.  Creates a matrix of the number of K-mers found given a GC count and a K-mer count.
   - **hist**:   Create an histogram of k-mer occurrences from a sequence file.  Adds metadata in output for easy plotting.
   - **filter**: Filtering tools.  Contains tools for filtering k-mer hashes and FastQ/A files:
     - **kmer**:         Produces a k-mer hash containing only k-mers within specified coverage and GC tolerances.
     - **seq**:          Filters a sequence file based on whether or not the sequences contain k-mers within a provided hash.
   - **plot**:   Plotting tools.  Contains several plotting tools to visualise K-mer and compare distributions. Requires gnuplot.  The following plot tools are available:
     - **density**:      Creates a density plot from a matrix created with the "comp" tool.  Typically this is used to compare two K-mer hashes produced by different NGS reads.
     - **profile**:      Creates a K-mer coverage plot for a single sequence.  Takes in fasta coverage output coverage from the "sect" tool
     - **spectra-cn**:   Creates a stacked histogram using a matrix created with the "comp" tool.  Typically this is used to compare a jellyfish hash produced from a read set to a jellyfish hash produced from an assembly. The plot shows the amount of distinct K-mers absent, as well as the copy number variation present within the assembly.
     - **spectra-hist**: Creates a K-mer spectra plot for a set of K-mer histograms produced either by jellyfish-histo or kat-histo.
     - **spectra-mx**:   Creates a K-mer spectra plot for a set of K-mer histograms that are derived from selected rows or columns in a matrix produced by the "comp".

In addition, KAT contains a python script for analysing the mathematical distributions present in the K-mer spectra in order to determine how much content is present in each peak.

This README contain some brief details of how to install and use KAT.  For more
extensive documentation please visit: https://readthedocs.org/projects/kat/


##Installation:

Generic installation description can be found in the INSTALL file. There are two ways to install KAT from source, either by cloning the git repository, or by downloading a distributable package, the later method is generally recommended as it reduces the number of installation steps and dependencies required to be on your system.

Installing from distributable:
 - Confirm dependencies are installed and configured:
  - GCC V4.8+
  - make
  - libtool V2.4.2+
  - Boost (system,filesystem,program_options,chrono,timer) V1.53+
  - Plotting engine:
   - Option 1 (preferred) python3 with matplotlib.  We recommend installing anaconda as this has all the required packages preinstalled.
   - Option 2 gnuplot 
 - Download tarball from https://github.com/TGAC/KAT/releases
 - Decompress and untar: ```tar -xvf kat-<version>.tar.gz```
 - Change into directory: ```cd kat-x.x.x```
 - Generate makefiles and confirm dependencies: ```./configure```
 - Compile software: ```make```
 - Run tests (optional) ```make check```
 - Install: ```sudo make install```


Installing from cloned repository
  - Clone the git repository (For ssh: ```git clone git@github.com:TGAC/KAT.git```; or for https: ```git clone https://github.com/TGAC/KAT.git```), into a directory on your machine.
  - "cd" into root directory of the installation
  - Ensure these tools are correctly installed and available on your system:
      - autoconf V2.53+
      - automake V1.11+
  - Create configuration script by typing: ```./autogen.sh```.
  - Follow all steps described in "Installing from a distributable" (except for the download and decompress tarball steps).```

The configure script can take several options as arguments.  One commonly modified option is ```--prefix```, which will install KAT to a custom directory.  By default this is "/usr/local", so the KAT executable would be found at "/usr/local/bin" by default.  In addition, some options specific to managing KAT dependencies located in non-standard locations are:

  - ```--with-boost``` - for specifying a custom boost directory

Type ```./configure --help``` for full details.

KAT can also make plots.  To enable plotting functionality we require either python3, with numpy, scipy and matplotlib installed, or gnuplot.  The python plotting method is the preferred method and will produce nicer results.  
If you don't already have python3 installed on your system we recommend installing anaconda3 as this contains everything you need.  

The type of plotting engine used will be determined when running the configure script, which will select the first engine detected in the following order: python, gnuplot, none.  
There is currently no way to select the plotting directory from a custom location, so the plotting system needs to be properly installed and configured on your system: i.e. python3 or gnuplot must be available on the PATH.


##Operating Instructions:

After KAT has been installed, the following tools should be available:

 - **kat** - a single executable binary file that contains a number of subtools.
 
Running ```kat --help``` will bring up a list of available tools within kat.  To get help on any of these subtools simple type: ```kat <tool> --help```.  For example: ```kat sect --help``` will show details on how to use the sequence coverage estimator tool.

KAT supports file globbing for input, this is particularly useful when trying to count and analyse kmers for paired end files.  For example,
assuming you had two files: LIB_R1.fastq, LIB_R2.fastq in the current directory then ```kat hist -C -m27 LIB_R?.fastq```, will consume any 
files matching the pattern LIB_R?.fastq as input, i.e. LIB_R1.fastq, LIB_R2.fastq.  The same result could be achieved listing the files at
the command line: ```kat hist -C -m27 LIB_R1.fastq LIB_R2.fastq```

Note, the KAT comp subtool takes 2 or three groups of inputs as positional arguments therefore we need to distinguish between the file groups.
This is achieved by surrounding any glob patterns or file lists in single quotes.  For example, assuming we have LIB1_R1.fastq, LIB1_R2.fastq, 
LIB2_R1.fastq, LIB2_R2.fastq in the current directory, and we want to compare LIB1 against LIB2, instead of catting the files together, we might 
run either: ```kat comp -C -D 'LIB1_R?.fastq' 'LIB2_R?.fastq'```; or ```kat comp -C -D 'LIB1_R1.fastq LIB1_R2.fastq' 'LIB2_R1.fastq LIB2_R2.fastq'```.
Both commands do the same thing.



##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

* Bernardo Clavijo
* Daniel Mapleson
* Darren Heavens
* Sarah Ayling
* Mario Caccamo
* George Kettleborough
* Gonzalo Garcia
* Jon Wright


See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)

We would also like to thank the authors of Jellyfish: https://github.com/gmarcais/Jellyfish; 
and SeqAn: http://www.seqan.de/.  Both are embedded inside KAT.

