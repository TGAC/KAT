.. KAT documentation master file, created by
   sphinx-quickstart on Thu Oct  8 16:26:54 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to KAT's documentation!
===============================

.. image:: images/kat_logo.png
    :scale: 30%

KAT is a suite of tools that generating, analysing and comparing k-mer spectra 
produced from sequence files (fasta or fastq).

.. toctree::
    :maxdepth: 2
    :numbered:
    
    installation
    using
    kmer
    walkthrough

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

We are currently planning to submit multiple publications around different aspects of KAT.
In the meantime if you use KAT in your work, please reference our github page `https://github.com/TGAC/KAT <https://github.com/TGAC/KAT>`_



.. _issues:

Issues
======

Should you discover any issues with spectre, or wish to request a new feature please raise a `ticket here <https://github.com/TGAC/KAT/issues>`_.
Alternatively, contact Daniel Mapleson at: daniel.mapleson@tgac.ac.uk


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

Last but not least a very special thanks to the Lab guys on their white coats, trying to make
sense of all our comments, giving us better data each day and trying to get into
our heads all the complex explanations for the biases and extra variability we
were finding day after day.


Credits
=======

 - Bernardo Clavijo (KAT's godfather, evangelist and all-round k-mer guru)
 - Daniel Mapleson (The software architect and developer)
 - George Kettleborough (For the recent python plotting functionality)
 - Gonzalo Garcia and Jon Wright (KAT superusers and providers of invaluable feedback)
