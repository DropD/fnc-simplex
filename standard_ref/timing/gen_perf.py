#!/usr/bin/env python

from math import *
import subprocess
import os
import pylab
from fncplot import fncplot

problemdir = "../../problems/gen/"
prog = "../bin/main"
#~ problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600, 1000]
problem_sizes = [30]


def run(prog):
    subprocess.call(prog, shell=True)
    with open('rdtsc', 'r')as f:
        lines = f.read().split('\n')
    data = []
    for line in lines:
        if line != "":
            data.append( line.split(',') )
    return data


files = [f for f in os.listdir(problemdir) if os.path.splitext(f)[1] == ".dlp"]
files.sort()
avgs = {}

for k in problem_sizes:
    avg = {};
    tok = "%04d" % (k)
    problems = [ f for f in files if tok in f ]
    for p in problems:
        print("_______________________\n"+problemdir + p)
        data = run(prog + " " + problemdir + p)
        for line in data:
            if not avg.get(line[0]):
                avg[line[0]] = 0
            avg[line[0]] += float(line[2])
    for key in avg:
        val = avg[key] / len(problems)
        if not avgs.get(key):
            avgs[key] = []
        avgs[key].append(val);

pylab.figure()
for key in avgs:
    pylab.plot(problem_sizes, avgs[key], label=key)
fncplot.title(r'Average performance', fontstyle='italic')
fncplot.xlabel('Number of variables $n$')
fncplot.ylabel('flop/cycle')
pylab.ylim([0,2])
pylab.grid(True)
pylab.legend(loc='center right')
#~ pylab.savefig('baseline_performance.png')

pylab.show()
