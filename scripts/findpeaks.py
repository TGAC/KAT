#!/usr/bin/env python3

import numpy as np

def findpeaks(a):
    a = np.squeeze(np.asarray(a))
    ad = np.sign(np.diff(a))
    # remove zeros to find end of plateaus
    ad[ad == 0] = 1
    return np.where(np.diff(ad) == -2)[0] + 1

