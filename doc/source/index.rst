.. KAT documentation master file, created by
   sphinx-quickstart on Thu Oct  8 16:26:54 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to KAT's documentation!
===============================

.. image:: images/kat_logo.png
    :scale: 30%

KAT provides a suite of tools that, through the use of k-mer counts, help the user address or identify issues such as determining sequencing completeness for assembly, assessing sequencing bias, identifying contaminants, validating genomic assemblies and filtering content.  KAT is geared primarily to work with high-coverage genomic reads from Illumina devices, although can work with any fasta or fastq sequence file.  

At it’s core KAT exploits the concept of k-mer spectra (histograms plotting number of distinct k-mers at each frequency).  By studying properties of the k-mer spectra it’s possible to discover important information about the data quality (level of errors, sequencing biases, completeness of sequencing coverage and potential contamination) and genomic complexity (size, karyotype, levels of heterozygosity and repeat content).  Further information can be gleaned through pairwise comparison of spectra, making KAT useful for WGS library comparisons and assembly validation.

The K-mer counting itself, a critical element for all KAT tools, is accomplished through an integrated and modified version of Jellyfish2's counting method http://www.genome.umd.edu/jellyfish.html.  We selected Jellyfish for this task because it supports large K values and is one of the fastest k-mer counting programs currently available.



.. toctree::
    :maxdepth: 2
    :numbered:
    
    installation
    using
    kmer
    walkthrough
    faq

.. _system:

System requirements
===================

KAT supports Unix, linux or Mac systems.  Windows, with something like cygwin, may work but hasn't been
tested.  A minimum of 8GB RAM, which will enable you to process small - medium sized datasets.  
Large datasets will require more RAM (potentially a lot more), the actual amount of
memory required depends on the size of the genome's to be processed, the k-mer size
selected and the size of your datasets.


.. _citing:

Citing
======

The KAT paper is currently in submission.  In the meantime, if you use our software
and wish to cite us please use our bioRxiv preprint:

Daniel Mapleson, Gonzalo Garcia Accinelli, George Kettleborough, Jonathan Wright, and Bernardo J. Clavijo. 
**KAT: A K-mer Analysis Toolkit to quality control NGS datasets and genome assemblies.**
Bioinformatics, 2016. `doi: 10.1093/bioinformatics/btw663 <http://bioinformatics.oxfordjournals.org/content/early/2016/10/20/bioinformatics.btw663.abstract>`_



.. _issues:

Issues
======

Should you discover any issues with KAT, or wish to request a new feature please raise a `ticket here <https://github.com/TGAC/KAT/issues>`_.
Alternatively, contact Daniel Mapleson at: daniel.mapleson@earlham.ac.uk; or Bernardo Clavijo at: bernardo.clavijo@earlham.ac.uk.  
However, please consult the `Frequently Asked Questions <faq>`_ page first in case your 
question is already answered there.



.. _availability:

Availability and License
========================

Open source code available on github: `https://github.com/TGAC/KAT.git <https://github.com/TGAC/KAT>`_

This documentation is hosted publicablly on read the docs: `https://kat.readthedocs.org/en/latest/ <https://kat.readthedocs.org/en/latest/>`_

KAT is available under `GNU GLP V3 <http://www.gnu.org/licenses/gpl.txt>`_


Acknowledgements
================

We owe a big acknowledgment to all TGAC staff that has been bored eternally
with k-mers, you have all been incredible patient and supportive with us.

Thanks to Mario Caccamo, Sarah Ayling, Federica Di Palma and David Swarbreck for 
all the support, feedback and encouragement.

Thanks to Richard Leggett, Daniel Zerbino and Zamin Iqbal for all the interesting 
discussions, comments and input.

Thanks to Dan Sargent for the use of his P.micrantha datasets for tests, and
their inclusion as figures on this document.

Thanks to all the KAT early adopters users who have provided invaluable feedback 
on the tool in its early stages: Paul Bailey, Jose De Vega, Rocio Enriquez-Gasca,
Marco Ferrarini and Dharanya Sampath.  And more recently, those from further afield
who have contributed on github.

A big thanks to the author of jellyfish, Guillaume Marcais. Jellyfish is fantastic
piece of software and is critical to enabling KAT to do what it does in an efficient
and timely fashion.

Last but not least a very special thanks to the Lab guys on their white coats, trying to make
sense of all our comments, giving us better data each day and trying to get into
our heads all the complex explanations for the biases and extra variability we
were finding day after day.


Credits
=======

 - Daniel Mapleson (The software architect and developer)
 - Gonzalo Garcia (KAT superuser and primary tester)
 - George Kettleborough (For the recent python plotting functionality)
 - Jon Wright (KAT superuser and documentation writer)
 - Bernardo Clavijo (KAT's godfather, evangelist and all-round k-mer guru)
 
