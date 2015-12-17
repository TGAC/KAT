.. _using:

Using KAT
================

KAT is a C++ program containing a number of subtools which can be used in
isolation or as part of a pipeline.  Typing ```kat --help``` will show a
list of the available subtools.  Each subtool has its own help system which you 
can access by typing ```kat <subtool> --help```.  


HIST
----

Creates a histogram with the number of distinct k-mers having a given frequency, 
derived from the input. The input can take the form of one or more FastA or FastQ 
files, or a jellyfish hash.  The last bucket in the histogram behaves as a catchall: 
it tallies all k-mers with a count greater or equal to the low end point of this bucket. 

This tool is very similar to the \"histo\" tool in jellyfish itself.  The primary 
difference being that the output contains metadata that make the histogram easier 
for the user to plot, and that this version is faster because we do not need to 
dump the hash to disk and read it back.

Basic usage:: 

    kat hist [options] (<input>)+


Applications:

 * Assess data quality: estimates of kmers deriving from errors; sequencing bias
 * Determine completeness of sequencing
 * Identify genomic properties: Heterozygosity, homozygosity, karyotype, repeat content.
 * Limited contamination detection



GCP
---

This tool takes in either a single jellyfish hash or one or more FastA or FastQ 
input files and then counts the GC nucleotides for each distinct K-mer in the hash.  
For each GC count and K-mer coverage level, the number of distinct K-mers are counted 
and stored in a matrix.  This matrix can be used in much that same way as a kmer
spectra histogram, although it provides richer output by incorperating GC content
into the picture.  This helps to distinguish legitimate content from contamination, 
which often appear as separate spots at unexpected GC and coverage levels.

Basic usage::

    kat gcp (<input>)+

Applications:

 * Assess data quality: estimates of kmers deriving from errors; sequencing bias
 * Determine completeness of sequencing
 * Identify genomic properties: Heterozygosity, homozygosity, karyotype, repeat content.
 * Contamination detection


Comp
----

Compares jellyfish K-mer count hashes.

The most common use case for this tool is to compare two (or three) K-mer hashes.  
The typical use case for this tool is to compare K-mers from two K-mer hashes 
both representing K-mer counts for reads.  However, it is also common to compare 
K-mers generated from reads to those generated from an assembly. If comparing 
K-mers from reads to K-mers from an assembly, the larger (most likely the read) 
K-mer hash should be provided first, then the assembly second. The third 
optional input acts as a filter, restricting the analysis to the K-mers present 
on that set.  The manual contains more details on specific use cases.

Basic usage::

    kat comp [options] <input_1> <input_2> [<input_3>]

Should the user wish to group multiple files to be concatenated into a single input 
group they may do so by surrounding the input group in single quotes.  The following
example groups the full input read set into the first input and compares against
an assembly::

    kat comp -t 8 -o pe_v_asm_test 'PE1.R1.fq PE1.R2.fq' asm.fa

... or more compactly::

    kat comp -t 8 -o pe_v_asm_test 'PE1.R?.fq' asm.fa


Applications:

 * Determine sequencing bias between left and right read pairs.
 * Compare the kmer spectrum of input reads against an assembly to gauge assembly completeness



SECT
----

Estimates coverage levels across sequences in the provided input sequence file.
This tool will produce a fasta style representation of the input sequence file 
containing K-mer coverage counts mapped across each sequence.  K-mer coverage is 
determined from the provided counts input file, which can be either one jellyfish 
hash, or one or more FastA / FastQ files.  In addition, a space separated table 
file containing the mean coverage score and GC of each sequence is produced.  The 
row order is identical to the original sequence file.

NOTE: K-mers containing any Ns derived from sequences in the sequence file not be 
included.

Basic usage::

    kat sect [options] <sequence_file> (<input>)+


Applications:

 * Analyse kmer coverage across assembled sequences
 * Compare assemblies using kmers, helpful to levels of contamination of a specific organism.
 * Contamination detection - Compare kmer spectrum against assembly providing average coverage and GC values for each contig, which can be 2D binned and plot as a heatmap


Filtering tools
---------------

KAT comes with two filtering tools allowing the user to slice and dice their data
in a rapid and simple way.


K-mer filtering
~~~~~~~~~~~~~~~

This tool allows the user to produce K-mer hashes, within and outside user defined 
GC and k-mer coverage bounds. This is useful for isolating k-mers that could be 
attributable to contamination, or for contamination removal.  Normally, the user 
would identify such regions using plots from the GCP tool.

Basic usage::

    kat filter kmer [options] (<input>)+

Applications:

 * Extracting k-mers with defined GC and coverage
 * Contamination extraction (from k-mer hash)


Sequence filtering
~~~~~~~~~~~~~~~~~~

The user loads a k-mer hash and then filters sequences (either in or out) depending 
on whether those sequences contain the k-mer or not.  The user can also apply a 
threshold requiring X% of k-mers to be in the sequence before filtering is applied.

Basic usage::

    kat filter seq [options] <seq_file> (<k-mer_hash>)

Applications:

 * Contamination extraction from read file or assembly file


Plotting tools
--------------

KAT comes with a selection of plotting tools for representing and comparing
K-mer spectra in various ways.  All plotting tools come with the ability to manually
modify axis, titles, limits, size, resolution, etc, although they will all try to pick 
intelligent defaults directly from the data provided.  


Density
~~~~~~~

Creates a scatter plot, where the density or "heat" at each point represents the 
number of distinct K-mers at that point.  Typically this is used to visualise a 
matrix produced by the ```kat comp``` tool to compare frequencies from two K-mer 
hashes produced by different NGS reads, or to visualise the GC vs K-mer matrices 
produced by the ```kat gcp``` tool.

Basic usage::

    kat plot density <matrix_file>


Applications:

 * Visualise GC vs coverage matrices
 * Visualise coverage vs coverage matrices


.. image:: images/ccoli_gcp.png
    :scale: 20%
.. image:: images/ccoli_comp.png
    :scale: 20%


Profile
~~~~~~~

Shows K-mer coverage level across an sequence

Basic usage::

    kat plot profile <sect_counts_file>

Applications:

 * Visualise coverage (and optionally GC) levels across a sequence or set of sequences

.. image:: images/profile.png
    :scale: 30%


Spectra_CN
~~~~~~~~~~

Shows K-mer duplication levels, which correspond to copy number variation within 
an assembly by comparing K-mers found in sequenced reads, to K-mers found in an 
assembly of those reads. Uses matrix output from the ```kat comp``` tool.

Basic usage::

    kat plot spectra-cn <matrix_file>

Applications:

 * Visualise the copy number spectra of WGS data compared against an assembly

.. image:: images/heterozygous_real.png
    :scale: 75%


Spectra_hist
~~~~~~~~~~~~

Visualises the K-mer spectra from ```kat hist``` or ```jellyfish histo``` output.  
This tool is designed to plot line graphs of one or more histograms.  The idea is 
to be able to compare total K-mer counts between different datasets.

Basic usage::

    kat plot spectra-hist <hist_file>

Applications:

 * Basic K-mer spectra visualisation

.. image:: images/ccoli_hist.png
    :scale: 20%



Spectra_mx
~~~~~~~~~~

Produces K-mer spectras from rows or columns in a matrix generated by ```kat comp```.  
This tool is designed to plot line graphs for one or more histograms, each histogram 
being represented by a single row or column in the matrix.

This tool also has a special mode for showing shared and exclusive content between 
two different samples. This mode takes the first row and column of the matrix representing 
content which is found exclusively in each sample.  Two more lines are plotting, 
one which has each following row summed, and the other that has each following column 
summed.  These two plots represent the shared content for each sample.  This mode 
can be activated using the ```--intersection``` flag.

Alternatively, you can select specific rows and columns from the matrix using a 
comma separated list identified with the ```--list``` option.  Each element in the 
list should start with either a 'c' or a 'r' indicating whether or not the column 
or row is requested.  Then the element should contain a number indicating which 
column or row to select.  For example: ```--list c0,r1``` will select column 0 and 
row 1. Note: spaces are not tolerated in this list.

Basic usage::

    kat plot spectra-mx <matrix_file>


Applications:

 * Visualising shared and exclusive content between two datasets
 * RNAseq to WGS comparison
 * Visualising k-mer spectra of arbitrary columns and rows from a matrix

.. image:: images/pe_v_pe_1_shared.png
    :scale: 50%

