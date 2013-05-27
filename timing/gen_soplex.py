#!/usr/bin/env python

from math import *
import subprocess
import os
import time
import pylab
from fncplot import fncplot

problemdir = "../problems/gen/"
prog = "../bin/main"
prog2 = "../../soplex-1.7.1/bin/soplex"
# problem_sizes should reflect the available *_nn_*.dlp problems
#~ problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600, 1000]
problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300]


def run(prog):
    start = time.time()
    subprocess.call(prog, shell=True)
    return (time.time() - start)


files = [f for f in os.listdir(problemdir) if os.path.splitext(f)[1] == ".dlp"]
files.sort()
files2 = [f for f in os.listdir(problemdir) if os.path.splitext(f)[1] == ".lp"]
files2.sort()
avgs = {}
avgs2 = {}

for k in problem_sizes:
    avg = {};
    avg2= {};
    tok = "%04d" % (k)
    problems = [ f for f in files if tok in f ]
    for p in problems:
        print("_______________________\n"+os.path.join(problemdir, p))
        t = run(prog + " " + os.path.join(problemdir, p))
        t2 = run(prog2 + " " + os.path.join(problemdir, p))
        for line in data:
            key = line[0]
            if not avg.get(key):
                avg[key] = 0
            if not avg2.get(key):
                avg2[key] = 0
            avg[key] += float(line[2])  # line[2] == fpc
            avg2[key] += float(line[3])  # line[3] == ci
    for key in avg:
        val = avg[key] / len(problems)
        val2 = avg2[key] / len(problems)
        if not avgs.get(key):
            avgs[key] = []
        if not avgs2.get(key):
            avgs2[key] = []
        avgs[key].append(val);
        avgs2[key].append(val2);

sizepi = [0.01, 3]
pi = [2, 2]
sizebeta = [2.**(-6), 0.5, 2**1]
beta = [2.**(-4), 2, 2**3]


fig = pylab.figure()
pylab.plot(sizepi, pi)
pylab.plot(sizebeta, beta)
for key in avgs:
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
