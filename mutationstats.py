#!/usr/bin/env python2.7
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
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from pandas import DataFrame, Panel
import cPickle

pp = PdfPages("plots.pdf")
pd.options.display.mpl_style="default"
pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)
pd.set_option("display.width", 1000)

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
mutrates = set()
for filename in argv[1:]:
    taxname, mutrate = extract(filename)
    mutrates.add(mutrate)
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
        index=list(mutrates))
wp = Panel(dfs)

with open("dict.pkl", 'wb') as f:
    cPickle.dump(d, f)

ax_spec = wp.minor_xs("Sensitivity").transpose().plot(kind="bar")
plt.legend(loc=3)
plt.title("Sensitivity")
plt.xlabel("Species")
plt.ylabel("Sensitivity")
plt.ylim((80, 100))
#plt.gcf().subplots_adjust(bottom=.25)
plt.setp(ax_spec.get_xticklabels(), rotation="horizontal", fontsize=20)
pp.savefig()

ax_fnr = wp.minor_xs("FNR").transpose().plot(kind="bar")
plt.legend(loc=0)
plt.title("False Negative Rate")
plt.xlabel("Species")
plt.ylabel("FNR")
plt.setp(ax_fnr.get_xticklabels(), rotation="horizontal", fontsize=20)
pp.savefig()
pp.close()

print "Total"
print wp.minor_xs("Total")
print "Discriminative"
print wp.minor_xs("Discriminative")
print "Sensitivity"
print wp.minor_xs("Sensitivity")
print "False Negative Rate"
print wp.minor_xs("FNR")
