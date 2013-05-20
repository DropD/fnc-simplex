import sys
import numpy as np

class TiergartenLP(object):
    def __init__(self, n, interest):
        self.n = n
        self.interest = interest
        self.required = np.random.uniform(0, 100, n)
        self.c = np.random.uniform(1, 1.5, n)
        self.b = np.zeros(n)
        self.A = np.zeros((n, n))

        intrpow = np.array([interest**(n - i - 1) for i in range(n)])
        for i in range(n):
            self.A[i,:i+1] = intrpow[n-i-1:]

        self.b = np.dot(self.A, self.required)

class PositiveLP(object):
    def __init__(self, n):
        self.n = n
        self.c = np.random.uniform(0, 1e6, n)
        self.b = np.random.uniform(0, 1e6, n)
        self.A = np.random.uniform(0, 1e6, (n,n))

class MixedLP(object):
    def __init__(self, n):
        self.n = n
        self.c = np.random.uniform(0, 1e6, n)
        self.b = np.random.uniform(0, 1e6, n)
        self.A = np.random.uniform(-1e6, 1e6, (n,n))

if __name__ == '__main__':
    import os.path as path
    from cplexlp import cplexlp
    from donjlp import donjlp
    from cpplexlp import cpplexlp

    usage = 'usage: python gen_lp.py <int: size> <path: output_file_name>'

    if len(sys.argv) != 3:
        print usage
        sys.exit()

    n_s, f = sys.argv[1:]
    n = int(n_s)

    f = path.splitext(f)[0]
    fclp = f+'.lp'
    fdlp = f+'.dlp'
    fcpplp = f+'.cpplp'

    interest = np.random.uniform(1.0, 1.1, 1)[0]

    #p = TiergartenLP(n, interest)
    p = PositiveLP(n)
    #p = MixedLP(n)

    with open(fclp, 'w') as out:
        out.write(str(cplexlp(searchList = [{'problem' : p}])))

    with open(fdlp, 'w') as out:
        out.write(str(donjlp(searchList = [{'problem' : p}])))

    with open(fcpplp, 'w') as out:
        out.write(str(cpplexlp(searchList = [{'problem' : p}])))
