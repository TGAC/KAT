
.. _faq:

Frequently Asked Questions
==========================

Can KAT handle compressed sequence files?
-----------------------------------------

Yes, via named pipes.  Anonymous named pipes (process substitution) is also supported.
For example, say we wanted to run ``kat hist`` using
gzipped paired end dataset, we can use a named pipe to do this as follows::

    mkfifo pe_dataset.fq && gunzip -c pe_dataset_?.fq.gz > pe_dataset.fq &
    kat hist -o pe_dataset.hist pe_dataset.fq

Where ``pe_dataset_?.fq.gz``, represents ``pe_dataset_1.fq.gz`` and ``pe_dataset_2.fq.gz``.

For those unfamiliar with named pipes, the first line will create an empty file
in your working directory called pe_dataset.fq and then specifies that anything 
consuming from the named pipe will take data that has been gunzipped first.  To be
clear this means you do not have to decompress the gzipped files to disk, this happens
on the fly as consumed by KAT.

Alternatively, using process substitution we could write the previous example more 
concisely in a single line like this::

    kat hist -o oe_dataset.hist <(gunzip -c pe_dataset_?.fq.gz)

As a more complex example, the KAT comp tool can be driven in spectra-cn mode using
both compressed paired end reads and a compressed assembly as follows::

    kat comp -o oe_spectra_cn <(gunzip -c pe_dataset_?.fq.gz) <(gunzip -c asm.fa.gz)

Thanks to John Davey and Torsten Seeman for suggesting this.


Why is jellyfish bundled with KAT?
----------------------------------

We require a stable interface to the k-mer hash arrays produced by jellyfish hence,
we are reliant on a particular version of jellyfish to guarantee that KAT works
correctly.  Instead of potentially requiring the user to install multiple jellyfish instances
on their machine, we bundle our own version, with all jellyfish binaries prefixed 
with `kat_` in order to avoid any naming clashes with official jellyfish releases.


I downloaded a release from github but it doesn't contain the configure script.  What gives?
--------------------------------------------------------------------------------------------

Github offers a feature which allows you to download source code bundles for all 
releases.  However, these bundles do not contain the distributable form of KAT, i.e.
they are not produced by calling ``make dist``.  Although they come from the same 
place you can distinguish between the github bundles and the proper distributable form
of KAT because it will have the name: ``kat-x.x.x.tar.gz``.
