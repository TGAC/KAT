.. _kmer:

K-mer spectra
=============

K-mer spectra is a representation of a dataset showing how many short
fixed length words (y-axis) appear a certain number of times (x-axis). The k-mer 
spectra is composed of distributions representing groups of motifs at different 
frequencies in the sample, plus biases. Given not too many biases, this becomes 
shape of the distributions provides a useful set of properties describing the 
biological sample and the sequencing process and the amount of useful data in the
dataset.

A typical nice 31-mer spectrum of S.cerevisae S288C WGS dataset is shown in the
following figure:

.. image:: images/kmer_spectra1.png
    :scale: 75%
    

This is composed of an error component containing a huge amount of
rare motifs, and a several other components as distributions with different modes
according to how many times a motif appear on the genome. The decomposition
showing this distributions can be seen here:

.. image:: images/kmer_spectra_breakdown.png
    :scale: 50%