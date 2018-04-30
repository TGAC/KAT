
.. _faq:

Frequently Asked Questions
==========================

Can KAT handle compressed sequence files?
-----------------------------------------

Yes, as of V2.4.0, KAT has native support for gzip decompression, so just treat
gzipped files as regular uncompressed fastq or fasta files.

If you wish to decompress other files such as bzip (or if you are using a pre V2.4.0 KAT), then
this is supported via named pipes.  Anonymous named pipes (process substitution)
is also supported.

For example, say we wanted to run ``kat hist`` using
bz2 paired end dataset, we can use a named pipe to do this as follows::

    mkfifo pe_dataset.fq && bzip2 -d -c pe_dataset_?.fq.bz2 > pe_dataset.fq &
    kat hist -o pe_dataset.hist pe_dataset.fq

Where ``pe_dataset_?.fq.bz2``, represents ``pe_dataset_1.fq.bz2`` and ``pe_dataset_2.fq.bz2``.

For those unfamiliar with named pipes, the first line will create an empty file
in your working directory called pe_dataset.fq and then specifies that anything
consuming from the named pipe will take data that has been gunzipped first.  To be
clear this means you do not have to decompress the gzipped files to disk, this happens
on the fly as consumed by KAT.

Alternatively, using process substitution we could write the previous example more
concisely in a single line like this::

    kat hist -o oe_dataset.hist <(bzip2 -d -c pe_dataset_?.fq.bz2)

As a more complex example, the KAT comp tool can be driven in spectra-cn mode using
both compressed paired end reads and a compressed assembly as follows::

    kat comp -o oe_spectra_cn <(bzip2 -d -c pe_dataset_?.fq.bz2) <(bzip2 -d -c asm.fa.bz2)

Thanks to John Davey and Torsten Seeman for suggesting this.


Why is jellyfish bundled with KAT?
----------------------------------

We require a stable interface to the k-mer hash arrays produced by jellyfish hence,
we are reliant on a particular version of jellyfish to guarantee that KAT works
correctly.  Instead of potentially requiring the user to install multiple jellyfish instances
on their machine, we bundle our own version, with all jellyfish binaries prefixed
with `kat_` in order to avoid any naming clashes with official jellyfish releases.
We have also made several modifications to jellyfish which make it more suitable
to processing via KAT.


Where's the distributable tarball for post-V2.4.0 releases?
-----------------------------------------------------------

As of V2.4.0, we no longer support installation via tarball.  We did this in order
to ensure boost is built alongside KAT, and this just didn't fit well into the
```make dist``` mechanism.  Please follow the new installation instructions and download
KAT via ```git clone```.


Should I dump jellyfish hashes to disk?
---------------------------------------

Most KAT tools have an option ```-d``` to write out the jellyfish hash to disk.  While it may seem
logical to do this if you need to reprocess the data, in normal usage do not recommend it.
Reading in K-mers straight from fastq or fasta, in most circumstances, will be faster than loading a 
jellyfish hash directly, particularly when using 4 or more threads (and especially, if the 
input files are gzipped).  Also storing jellyfish hashes can consume a large amount of storage space.

