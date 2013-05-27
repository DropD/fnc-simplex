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

with open(os.path.join(outdir, 'fpc_avg')) as fpcpi:
    avgs = cPickle.load(fpcpi)
with open(os.path.join(outdir, 'ci_avg')) as cipi:
    avgs2 = cPickle.load(cipi)

sizepi = [0.01, 3]
pi = [2, 2]
sizebeta = [2.**(-6), 0.5, 2**1]
beta = [2.**(-4), 2, 2**3]


fig = pylab.figure()
pylab.plot(sizepi, pi)
pylab.plot(sizebeta, beta)
#~ keep = [ 'avx', 'baseline', 'block_swap', 'block_avx', 'block_2', 'soplex', 'gurobi', 'glpk', 'ssa' ]
for key in avgs:
  if key not in [ 'avx', 'array', 'block_swap', 'block_avx', 'block_2', 'sse', 'block-sse', 'nta' ]:
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
