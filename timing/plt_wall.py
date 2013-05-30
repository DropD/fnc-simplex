#!/usr/bin/env python

from math import *
import os
import pylab
from fncplot import fncplot
import cPickle

# algorithms filtered by intersection of keep and ignore, keep can be left undefined
#~ keep = [ ]
keep = [ 'baseline', 'soplex', 'gurobi', 'glpk' ]
ignore = [ ]

outdir = "tmp"
prog = "../bin/main"

if len(sys.argv) > 1:
    outdir = sys.argv[1]

with open(os.path.join(outdir, 'problem_sizes')) as cycpi:
    problem_sizes = cPickle.load(cycpi)

with open(os.path.join(outdir, 'wall_avg')) as wallpi:
    avgs = cPickle.load(wallpi)
print('Available algorithms: ' + str([key for key in avgs]) )

try:
  keep
except NameError:
  keep = avgs


pylab.figure()

for key in sorted(avgs.keys()):
    if key in keep and key not in ignore:
        print('Plotting ' + key)
        pylab.plot(problem_sizes, avgs[key], label=key)

fncplot.title(r'Average Runtime', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('Walltime')
pylab.yscale('log')
pylab.grid(True)
pylab.legend(loc='lower right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
