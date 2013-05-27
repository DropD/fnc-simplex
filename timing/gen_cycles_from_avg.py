#!/usr/bin/env python

from math import *
import subprocess
import os
import pylab
from fncplot import fncplot
import cPickle

outdir = "tmp"
problemdir = "../problems/gen/"
prog = "../bin/main"
# problem_sizes should reflect the available *_nn_*.dlp problems
problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600, 1000]


with open(os.path.join(outdir, 'cycles_avg')) as cycpi:
    avgs = cPickle.load(cycpi)

pylab.figure()
for key in avgs:
    pylab.plot(problem_sizes, avgs[key], label=key)
fncplot.title(r'Average Runtime', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('Cycles')
pylab.yscale('log')
pylab.grid(True)
pylab.legend(loc='lower right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
