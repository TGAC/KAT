#!/bin/bash

# Assumes jellyfish, gnuplot and kat are on the path

# Arguments
# 1 - threads
# 2 - output prefix
# 3 - read1 fastq
# 4 - read2 fastq
# 5 - z_cap for plotting

jellyfish count -c16 -t $1 -s 10000000000 -C -m31 -o $3.jf31 $3;
jellyfish count -c16 -t $1 -s 10000000000 -C -m31 -o $4.jf31 $4;

kat comp -t $1 -o $2 $3.jf31_0 $4.jf31_0;

kat plot flame -o $2_main.mx.png -x $3 -y $4 -z $5 $2_main.mx;
