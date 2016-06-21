#! /bin/sh

. ./compat.sh

$KAT comp -m13 -o temp/spectra-cn_test ${data}/ecoli_r1.1K.fastq ${data}/EcoliK12.fasta
$KAT comp -m13 -v -n -o temp/density_test ${data}/ecoli_r1.1K.fastq ${data}/ecoli_r2.1K.fastq
$KAT comp -m13 -o temp/glob_test ${data}'/ecoli_r?.1K.fastq' ${data}/EcoliK12.fasta

