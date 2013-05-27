#!/usr/bin/env python

from math import *
import subprocess
import os
import pylab
from fncplot import fncplot
import cPickle

problemdir = "../problems/gen/"
prog = "../bin/main"
# problem_sizes should reflect the available *_nn_*.dlp problems
problem_sizes = [10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600, 1000]


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
avgs0 = {}
avgs1 = {}
avgs2 = {}
avgs3 = {}

for k in problem_sizes:
    avg0 = {};
    avg1 = {};
    avg2= {};
    avg3= {};
    tok = "%04d" % (k)
    problems = [ f for f in files if tok in f ]
    for p in problems:
        print("_______________________\n"+os.path.join(problemdir, p))
        data = run(prog + " " + os.path.join(problemdir, p))
        for line in data:
            key = line[0]
            if not avg0.get(key):
                avg0[key] = 0
            if not avg1.get(key):
                avg1[key] = 0
            if not avg2.get(key):
                avg2[key] = 0
            if not avg3.get(key):
                avg3[key] = 0
            avg0[key] += float(line[1])  # line[1] == cycles
            avg1[key] += float(line[2])  # line[2] == fpc
            avg2[key] += float(line[3])  # line[3] == ci
            avg3[key] += float(line[4])  # line[4] == walltime
    for key in avg:
        val0 = avg0[key] / len(problems)
        val1 = avg1[key] / len(problems)
        val2 = avg2[key] / len(problems)
        val3 = avg3[key] / len(problems)
        if not avgs0.get(key):
            avgs0[key] = []
        if not avgs1.get(key):
            avgs1[key] = []
        if not avgs2.get(key):
            avgs2[key] = []
        if not avgs3.get(key):
            avgs3[key] = []
        avgs0[key].append(val0);
        avgs1[key].append(val1);
        avgs2[key].append(val2);
        avgs3[key].append(val3);

cPickle.dump(avgs0, 'cycles_pickle')
cPickle.dump(avgs1, 'fpc_pickle')
cPickle.dump(avgs2, 'ci_pickle')
cPickle.dump(avgs3, 'wall_pickle')
