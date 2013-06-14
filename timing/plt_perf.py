#!/usr/bin/env python

from math import *
import os
import sys
import pylab
from fncplot import fncplot
import cPickle

# algorithms filtered by intersection of keep and ignore, keep can be left undefined
#~ keep = [ 'baseline', 'avx', 'block2_swap', 'block2_swap_avx', 'block4_swap_avx', 'block8_swap_avx', 'nta' ]
#~ keep = [ 'baseline', 'avx', 'block2_swap', 'block2_swap_avx', 'block4_swap_avx', 'block8_swap_avx', 'nta', 'soplex', 'gurobi' ]
#~ keep = [ 'baseline', 'soplex', 'gurobi', 'glpk' ]
#~ keep = [ 'baseline', 'block16x16_swap', 'block16x1_swap', 'block16x8_swap', 'block2x8_swap', 'block1x4_swap', 'block2x4_swap' ]
ignore = [ ]
#~ ignore = [ 'block2_swap', 'avx' ]
#~ ignore = [ 'block2_swap', 'avx' ]

outdir = "tmp"
#~ outdir = "tmp_donj_upper_0_O3novec"
prog = "../bin/main"

if len(sys.argv) > 1:
    outdir = sys.argv[1]

with open(os.path.join(outdir, 'problem_sizes')) as cycpi:
    problem_sizes = cPickle.load(cycpi)

with open(os.path.join(outdir, 'fpc_avg')) as fpcpi:
    avgs = cPickle.load(fpcpi)
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

fncplot.title(r'Average performance', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('flop/cycle')
pylab.ylim([0,2])
#pylab.xscale('log')
pylab.grid(True)
pylab.legend(loc='upper right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
