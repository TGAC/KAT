#! /bin/sh

. ./compat.sh

$KAT gcp -m17 -o temp/gcp_test data/ecoli_r?.1K.fastq
