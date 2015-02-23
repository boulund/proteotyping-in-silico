#!/usr/bin/env python
# Fredrik Boulund 2015
# Create summary statistics of proteotyping in-silico mutation study

from __future__ import division
from sys import argv, exit

# if using SEABORN
# import seaborn as sns
#sns.set_context("poster", font_scale=2.0)
#sns.set_context({"legend.frameon":True})

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame, Panel
import cPickle

pd.options.display.mpl_style="default"

if len(argv) < 2:
    print "usage: script.py file.results [...]"
    exit()

def extract(filename):
    """Extract juicy bits from filename"""
    name = filename.split("/")[-1].split(".fasta.results")[0]
    taxname = name.split("_",1)[0]
    if "p0." in name:
        mutrate = float("0."+ name.split("p0.")[1].split(".",1)[0])
    else:
        mutrate = 0.0
    return taxname, mutrate

d = {}
for filename in argv[1:]:
    taxname, mutrate = extract(filename)
    if not d.has_key(taxname):
        d[taxname] = {}
    with open(filename) as f:
        [f.readline() for i in xrange(0,3)]
        line = f.readline()
        sensitivity, M, _ = line.split(None, 2)
        for line in f:
            if line.startswith(" Total:"):
                N = line.split()[1]
                break
    sensitivity = float(sensitivity)
    M = int(M)
    N = int(N)
    FNR = (N-M)/M
    if not d[taxname].has_key("FNR"):
        d[taxname]["FNR"] = [FNR]
        d[taxname]["Sensitivity"] = [sensitivity]
        d[taxname]["Total"] = [N]
        d[taxname]["Discriminative"] = [M]
    else:
        d[taxname]["FNR"].append(FNR)
        d[taxname]["Sensitivity"].append(sensitivity)
        d[taxname]["Total"].append(N)
        d[taxname]["Discriminative"].append(M)

dfs = {}
for taxname in d.keys():
    dfs[taxname] = DataFrame(d[taxname], 
        columns=d[taxname].keys(), 
        index=[0.0, 0.01, 0.02, 0.03, 0.05, 0.10])
wp = Panel(dfs)
wp.to_pickle("panel.pkl")

with open("dict.pkl", 'wb') as f:
    cPickle.dump(d, f)

plot_spec = wp.minor_xs("Sensitivity").transpose().plot(kind="bar")
plt.legend(loc=3)
plt.title("Sensitivity")
plt.xlabel("Mutation rate")
plt.ylabel("Sensitivity")
plt.ylim((80, 100))
plt.gcf().subplots_adjust(bottom=.25)

fnr = wp.minor_xs("FNR")
fnr.plot(kind="bar")
plt.legend(loc=0)
plt.title("False Negative Rate")
plt.xlabel("Mutation rate")
plt.ylabel("FNR")

print "Total"
print wp.minor_xs("Total")
print "Discriminative"
print wp.minor_xs("Discriminative")
print "Sensitivity"
print wp.minor_xs("Sensitivity")
print "False Negative Rate"
print wp.minor_xs("FNR")
