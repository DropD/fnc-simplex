import ethpy as ep
import numpy as np
from ethpy import fncplot as fpl
import matplotlib.pyplot as plt

data = []
data.append(ep.read('0008.dat'))
data.append(ep.read('0016.dat'))
data.append(ep.read('0032.dat'))
data.append(ep.read('0064.dat'))
data.append(ep.read('0128.dat'))
data.append(ep.read('0256.dat'))
data.append(ep.read('0512.dat'))
data.append(ep.read('1024.dat'))

x = [8, 16, 32, 64, 128, 256, 512, 1024]

r = np.array([np.average(i, 0) for i in data])
s = np.array([np.std(i, 0) for i in data])

plt.figure()
plt.semilogy(x, r[:,0], label='baseline')
plt.semilogy(x, r[:,1], label='glpk')
fpl.xlabel('problem size')
fpl.ylabel('cycles')
fpl.title('Baseline')
plt.legend()

plt.savefig('baseline.pdf')
