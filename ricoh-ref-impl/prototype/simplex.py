import numpy as np
from numpy import *

def fstr(flt):
    return '{:> .3f}'.format(flt)

def vstr(vec):
    return '[{}]'.format(', '.join([fstr(i) for i in vec]))

def print_stf(c, A, b):
    print 'min '+' '.join(['{0:> .3f} x{1}'.format(c[i], i) for i in range(len(c))])
    print 's.t.'
    for i in range(len(b)):
        print ' + '.join(['{0:> .3f} x{1}'.format(A[i,j], j) for j in range(len(c))]) + ' = {0: .3f}'.format(b[i])
    print ', '.join(['x{0}'.format(i) for i in range(len(b))]) + ' >= 0'

def make_std_form(c, A, b):
    print 'Warning: transformation not implemented. input must be in standard form:'
    print "min c' x s.t."
    print "A x == b"
    print "  x >= 0"
    print '\nIs this correct?'
    print_stf(c, A, b)
    return c, A, b

def SimpleX(c, A, b):
    eps = 1e-10
    # std form
    c, A, b = make_std_form(c, A, b)

    # dimensions
    n = len(c)
    m = len(b)
    assert A.shape[0] == m and A.shape[1] == n, 'Wrong Matrix shape for A'

    # initial basis
    B = [i for i in range(n-m, n)]
    N = [i for i in range(n-m)]
    Ai = A[:,B]
    x = zeros(n)
    x[B] = b

    # step
    while True:
        p = dot(c[B], Ai)
        c_bar = c - dot(p, A)
        if not any(c_bar < -eps):
            return x
        j = [i for i in range(len(c_bar)) if c_bar[i] < 0][0]
        u = dot(Ai, A[:,j])
        if all(u <= 0):
            print '\nobjective is -inf'
            return array([])
        xbu = [u[i] > eps and x[B][i] / u[i] or inf for i in range(m)]
        ll = argmin(xbu)
        d = xbu[ll]
        l = B[ll]

        x[B] = x[B] - d * u

        B.remove(l)
        N.append(l)

        N.remove(j)
        B.append(j)

        B.sort()
        N.sort()

        x[j] = d

        Aa = dot(Ai, A[:,B])
        Q = identity(m)
        for i in range(m):
            Q[i, ll] = - (u[i] / u[ll])
        Q[ll, ll] = 1 / u[ll]
        Ai = dot(Q, Ai)

if __name__ == '__main__':
    # test system, solution is 164.699
    # min c' x; Ax == b; x >= 0
    test_c = np.array([ -3.0, -9.0, 0.0, 0.0,   0.0])
    test_A = np.array([
                       [ 3.0,  1.0, 1.0, 0.0, 0.0], 
                       [ 2.0,  3.0, 0.0, 1.0, 0.0], 
                       [-2.0,  3.0, 0.0, 0.0, 1.0],
                      ])
    test_b = np.array([ 15.0, 18.0, 6.0])

    x = SimpleX(test_c, test_A, test_b)
    print '\ndone.'
    print 'x = {0}'.format(vstr(x))
