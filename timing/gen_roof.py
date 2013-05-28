#!/usr/bin/env python

from math import *
import os
import pylab
from fncplot import fncplot
import cPickle

# algorithms filtered by intersection of keep and ignore, keep can be left undefined
keep = [ 'avx', 'baseline', 'block_swap', 'block_avx', 'block_2', 'soplex', 'gurobi', 'glpk' ]
ignore = [ ]

outdir = "tmp"
prog = "../bin/main"

with open(os.path.join(outdir, 'problem_sizes')) as cycpi:
    problem_sizes = cPickle.load(cycpi)

with open(os.path.join(outdir, 'fpc_avg')) as fpcpi:
    avgs = cPickle.load(fpcpi)
with open(os.path.join(outdir, 'ci_avg')) as cipi:
    avgs2 = cPickle.load(cipi)
print('Available algorithms: ' + str([key for key in avgs]) )

try:
  keep
except NameError:
  keep = avgs


sizepi = [0.01, 3]
pi = [2, 2]
sizebeta = [2.**(-6), 0.5, 2**1]
beta = [2.**(-4), 2, 2**3]


fig = pylab.figure()
pylab.plot(sizepi, pi)
pylab.plot(sizebeta, beta)

for key in avgs:
    if key in keep and key not in ignore:
        print('Plotting ' + key)
        pylab.plot(avgs2[key], avgs[key], label=key)

fncplot.title(r'Roofline plot', fontstyle='italic')
fncplot.xlabel('Operational intensity')
fncplot.ylabel('flop/cycle')
ax = fig.add_subplot(111)
ax.xaxis.set_major_locator(pylab.LinearLocator())
ax.yaxis.set_major_locator(pylab.LinearLocator())
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
pylab.xlim([2.**(-3), 2.5])
pylab.ylim([2.**(-3), 2.5])

pylab.grid(True)
pylab.legend(loc='lower right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
