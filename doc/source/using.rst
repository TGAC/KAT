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

Usage: kat hist [options] (<input>)+


Applications::

 * Assess data quality: estimates of kmers deriving from errors; sequencing bias
 * Determine completeness of sequencing
 * Identify genomic properties: Heterozygosity, homozygosity, karyotype, repeat content.
 * Limited contamination detection


Runtime:



GCP
---

This tool takes in either a single jellyfish hash or one or more FastA or FastQ 
input files and then counts the GC nucleotides for each distinct K-mer in the hash.  
For each GC count and K-mer coverage level, the number of distinct K-mers are counted 
and stored in a matrix.  This matrix can be used in much that same way as a kmer
spectra histogram, although it provides richer output by incorperating GC content
into the picture.  This helps to distinguish legitimate content from contamination, 
which often appear as separate spots at unexpected GC and coverage levels.

Usage: kat gcp (<input>)+

Applications::

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

Usage: kat comp [options] <input_1> <input_2> [<input_3>]

TODO: Describe how to do globbing with comp here using quotes to define input groups.


Applications::

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

Usage: kat sect [options] <sequence_file> (<input>)+


Applications::

 * Analyse kmer coverage across assembled sequences
 * Compare assemblies using kmers, helpful to levels of contamination of a specific organism.
 * Contamination detection - Compare kmer spectrum against assembly providing average coverage and GC values for each contig, which can be 2D binned and plot as a heatmap