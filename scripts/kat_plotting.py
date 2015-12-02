#!/usr/bin/env python3

import numpy as np

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

def findpeaks(a):
    a = np.squeeze(np.asarray(a))
    ad = np.sign(np.diff(a))
    # remove zeros to find end of plateaus
    ad[ad == 0] = 1
    return np.where(np.diff(ad) == -2)[0] + 1

