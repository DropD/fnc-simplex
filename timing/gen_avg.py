#!/usr/bin/env python

from math import *
import subprocess
import os
import sys
import cPickle
import re

problemdir = "../problems/gen_light/"
#~ problemdir = "../problems/gen_std"
#~ problemdir = "../problems/gen_heavy/"
#~ problemdir = "../problems/gen_upper/"

if len(sys.argv) > 1:
    problemdir = sys.argv[1]

outdir = "tmp"
prog = "../bin/main"

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
problem_sizes = [ re.search(r"_[0-9]+_", f).group() for f in files ]
problem_sizes = list(set(problem_sizes))  # make unique
problem_sizes.sort()
#~ print(problem_sizes)

avgs0 = {}
avgs1 = {}
avgs2 = {}
avgs3 = {}

for token in problem_sizes:
    avg0 = {};
    avg1 = {};
    avg2= {};
    avg3= {};
    problems = [ f for f in files if token in f ]
    for p in problems:
        print("_______________________\n"+os.path.join(problemdir, p)+"\n")
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
    for key in avg0:
        print key
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


problem_sizes = [ s.replace('_', '') for s in problem_sizes ]
with open(os.path.join(outdir, 'problem_sizes'), 'w') as cycpi:
    cPickle.dump(problem_sizes, cycpi)

with open(os.path.join(outdir, 'cycles_avg'), 'w') as cycpi:
    cPickle.dump(avgs0, cycpi)
with open(os.path.join(outdir, 'fpc_avg'), 'w') as fpcpi:
    cPickle.dump(avgs1, fpcpi)
with open(os.path.join(outdir, 'ci_avg'), 'w') as cipi:
    cPickle.dump(avgs2, cipi)
with open(os.path.join(outdir, 'wall_avg'), 'w') as wapi:
    cPickle.dump(avgs3, wapi)
