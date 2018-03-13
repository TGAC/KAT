import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap

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

def correct_filename(filename):
    split = filename.split('.')
    if len(split) > 1:
        ext = split[-1]
    else:
        ext = ''

    types = list(plt.gcf().canvas.get_supported_filetypes().keys())
    if ext in types:
        return filename
    elif "png" in types:
        return filename + ".png"
    elif "pdf" in types:
        return filename + ".pdf"
    else:
        return filename + "." + types[0]

def wrap(name):
    return "\n".join(textwrap.wrap(name, 60))

