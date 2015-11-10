#!/usr/bin/env python

def readheader(input_file):
    header = {}
    for line in input_file:
        if line[0:2] == "# ":
            s = line[2:-1].split(":")
            n = s[0]
            v = ":".join(s[1:])
            header[n] = v
        elif line[:-1] == "###":
            break
        else:
            break
    return header
