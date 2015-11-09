#!/usr/bin/env python

import numpy as np

def findpeaks(a):
    ad = np.sign(np.diff(a))
    # remove zeros to find end of plateaus
    ad[ad == 0] = 1
    return np.where(np.diff(ad) == -2)[0]

