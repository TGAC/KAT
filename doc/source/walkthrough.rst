.. _walkthrough:

KAT Walkthrough
===============

This part of the manual describes particular use cases where KAT can be useful:


Comparing R1 v R2 in an Illumina PE dataset
-------------------------------------------



Detecting GC bias
-----------------



Checking library consistency
----------------------------

PE vs PE
~~~~~~~~


PE vs LMP
~~~~~~~~~



Contamination detection and extraction
--------------------------------------

Breaking WGS data into k-mers provides a nice way of identifying contamination, or
otherwise unexpected content, in your reads or assemblies.  This section will walk
you through how you might be able to identify and extract contamination in your 
data.


In reads
~~~~~~~~

Contamination detection in your WGS datasets are reliant on the contamination having
differing levels of coverage and/or GC value from your target species.  The process
in KAT is simple::

    kat gcp [options] (WGS_file)+

Running this tool will produce a matrix containing distinct k-mer counts at varying 
frequency and GC value.  It will also produce a plot, such as the one shown here:

.. figure:: images/contaminant_MP.png
    :scale: 50%
    :alt: Contamination detection 1
    :align: center
    :figclass: align-center

    KAT GCP output run through the density plotting tool.  Error k-mers shown at
    very low level with wide GC spread, genuine content between 10-100X with GC 
    spread from 5-25, unexpected content shown at approx 200X with GC 15-25.

The high coverage hot-spot is already suspicious but it becomes even more so after
consider other WGS libraries of the same sample:

.. image:: images/contaminant_ope1.png
    :scale: 33%
.. image:: images/contaminant_ope2.png
    :scale: 33%
.. image:: images/contaminant_PE.png
    :scale: 33%

No other library contains the such a hotspot at GC 15-25.  After merging all libraries
into one the contaminant becomes obvious as the coverage has not altered, meaning
that k-mers in that region where not also found in the other libraries:

.. figure:: images/contaminant_all.png
    :scale: 50%
    :alt: Contamination detection 2
    :align: center
    :figclass: align-center

    After merging all WGS libraries of the same sample, the original suspicious
    region has not altered in coverage.

We can then use the filtering tools in KAT to extract k-mers inside, or outside
defined coverage and GC limits.  In this case we could take the original MP library
and run the following command::

    kat filter kmer --low_count=100 --high_count=250 --low_gc=13 --high_gc=25 <path_to_MP_lib>

This produces a k-mer hash containing only those k-mer found in the defined region.
We can then get the reads (or assembled contigs) associated with these k-mers by
running the following command::

    kat filter seq --threshold=0.5 <path_to_seq_file_to_filter> <filtered_k-mer hash>

BLASTing some of those filtered sequences might then identify the contaminant.


In assemblies
~~~~~~~~~~~~~

Detecting contaminants in assemblies involves a similar process to that described 
in the previous section.  It involves marking contigs in an assembly their average 
k-mer coverage and GC%.  

To obtain the average coverage and GC% scores for each contig use the following
command::

    kat sect [options] <assembly> (WGS_data)+

By extracting the median coverage and gc% columns from the stats file it is possible
to create a scatter plot which can be used in a similar way to that described in
the previous section.

A second use case assumes you already know the contaminant genome and have
access to the reference assembly of that contaminant.  In this case you can 
directly inspect your assembly for signs of the contaminant using the following command::

    kat sect [options] <assembly> <contaminant_genome>

This counts k-mers in the contaminant genome and applies them to sequences in your
assembly.  By reverse sorting the stats file produced by the "%_non_zero_corrected" column
you can identify contigs belonging to the contaminant.  Normally, assuming the 
contaminant is the exact same species as that found in your assembly you expect
to see very high percentage scores (> 90%).  Moderate scores (20-80%) might indicate
either some shared content or chimeric sequences and should be investigated more
thoroughly.



Finding repetitive genomic regions
----------------------------------

Sometimes it's useful to identify regions in a dataset that are repetitive.  This
can easily be done with the following command::

    kat sect -E -F [options] <genome_file> <genome_file>

This counts k-mers in the provided assembly and then marksup the assembly with
those k-mer counts at each position.  Regions that have a count of 1 are extracted
into a new FastA file containing non-repetitive content.  Regions that have a count
of 2-20 (maximum threshold can be adjusted) are extracted to FastA file containing
the repetitive content.


Checking RNAseq consistency with genomic content
------------------------------------------------



Assembly analysis using k-mer spectra
-------------------------------------
