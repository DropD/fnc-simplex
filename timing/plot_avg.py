import os
from matplotlib import pyplot as plt
from fncplot import fncplot as fpl
import cPickle

avgdir = "tmp"

def load_avg(avgfile):
    with open(os.path.abspath(os.path.join(avgdir, avgfile))) as avgpi:
        avgs = cPickle.load(avgpi)
    return avgs

def plt_perf():
    fpl.title(r'Average Performance', fontstyle='italic')
    fpl.xlabel('Number of variables $n$')
    fpl.ylabel('flop/cycle')
    plt.ylim([0,2])
    plt.grid(True)
    plt.legend(loc='upper right')

def plt_cycles():
    fpl.title(r'Average Runtime', fontstyle='italic')
    fpl.xlabel('Number of variables $n$')
    fpl.ylabel('Cycles')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='lower right')

def plt_wall():
    fpl.title(r'Average Runtime', fontstyle='italic')
    fpl.xlabel('Number of variables $n$')
    fpl.ylabel('Walltime')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc='lower right')

def plt_roof():
    sizepi = [0.01, 3]
    pi = [2, 2]
    sizebeta = [2.**(-6), 0.5, 2**1]
    beta = [2.**(-4), 2, 2**3]

    plt.plot(sizepi, pi)
    plt.plot(sizebeta, beta)

    fpl.title(r'Roofline plot', fontstyle='italic')
    fpl.xlabel('Operational intensity')
    fpl.ylabel('flop/cycle')
    ax = plt.gcf().get_axes()[0]
    ax.xaxis.set_major_locator(plt.LinearLocator())
    ax.yaxis.set_major_locator(plt.LinearLocator())
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    plt.xlim([2.**(-3), 2.5])
    plt.ylim([2.**(-3), 2.5])
    plt.grid(True)
    plt.legend(loc='lower right')

plt_what = {
            'perf' : plt_perf,
            'cycles' : plt_cycles,
            'wall' : plt_wall,
            'roof' : plt_roof
           }

def plot_avg(avgs, what=['perf'], who=['baseline']): 
    n = avgs['problem_sizes']
    for w in what:
        avg = avgs[w]
        plt.figure()
        for impl in who:
            if w == 'roof':
                plt.plot(avg[impl], avgs['perf'][impl], label=impl)
            else:
                plt.plot(n, avg[impl], label=impl)
        plt_what[w]()
        plt.show()
