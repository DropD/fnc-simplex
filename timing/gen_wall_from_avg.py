#!/usr/bin/env python

from math import *
import subprocess
import os
import pylab
from fncplot import fncplot
import cPickle

outdir = "tmp"
#~ problemdir = "../problems/gen/"
problemdir = "../problems/gen_verylight/"
prog = "../bin/main"
# problem_sizes should reflect the available *_nn_*.dlp problems
problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600, 1000]

with open(os.path.join(outdir, 'wall_avg')) as wallpi:
    avgs = cPickle.load(wallpi)

pylab.figure()
#~ keep = [ 'avx', 'baseline', 'block_swap', 'block_avx', 'block_2', 'soplex', 'gurobi', 'glpk', 'ssa' ]
for key in avgs:
  if key not in [ 'array', 'block_2', 'sse', 'block-sse', 'nta' ]:    pylab.plot(problem_sizes, avgs[key], label=key)
fncplot.title(r'Average Runtime', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('Walltime')
pylab.yscale('log')
pylab.grid(True)
pylab.legend(loc='lower right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
