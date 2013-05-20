#!/usr/bin/env python

from math import *
import subprocess
import os
import pylab
from fncplot import fncplot

problemdir = "../../problems/gen/"
prog = "../bin/main"
problem_sizes = [10, 20, 50, 100, 500, 1000]


def run(prog):
    subprocess.call(prog, shell=True)
    f = open("rdtsc", 'r')
    c = float(f.readline())
    fpc = float(f.readline())
    f.close()
    return fpc


files = [f for f in os.listdir(problemdir) if os.path.splitext(f)[1] == ".dlp"]
files.sort()
avgs = []

for k in problem_sizes:
    avg = 0;
    tok = "%04d" % (k)
    problems = [ f for f in files if tok in f ]
    for p in problems:
        print("_______________________\n"+problemdir + p)
        fpc = run(prog + " " + problemdir + p)
        avg += fpc
    avg /= len(problems)
    avgs.append(avg)



pylab.figure()
pylab.plot(problem_sizes, avgs, label="baseline")
fncplot.title(r'Averaged performance', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('flop/cycle')
pylab.ylim([0,2])
pylab.grid(True)
pylab.legend(loc='center right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
