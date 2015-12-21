#! /bin/sh

. ./compat.sh

$KAT sect -o temp/sect_length ${data}/sect_length_test.fa ${data}/ecoli.header.jf27
$KAT sect -o temp/sect_test ${data}/sect_length_test.fa ${data}/ecoli.header.jf27
