from timerH import timerH
from timerC import timerC
import sys
from os import path

class Func(object):
    def __init__(self, fname, hname):
        self.name = fname
        self.header = path.split(hname)[1]
        self.signature = None
        self.args = None
        with open(hname, 'r') as header:
            for line in header:
                line = line.strip()
                fstart = line.find(fname)
                if(fstart == -1):
                    pass
                else:
                    self.rtype = line[:fstart].strip()
                    sigstart = line.find('(') + 1
                    sigend   = line.find(')')
                    self.signature = line[sigstart:sigend]
                    break
            if not self.signature:
                raise Exception('function not found in header!')

        sigspl = self.signature.split()
        argl = [i[:-1] for i in sigspl if i[-1] == ',']
        argl.append(sigspl[-1])
        self.args = ', '.join(argl)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: python gen-timer.py <function> <header> <prefix>'
        sys.exit()

    fname, hname = sys.argv[1:3]
    prefix = sys.argv[3]
    f = Func(fname, hname)

    thn = 'include/time_{0}.hpp'.format(fname)
    tcn = 'src/time_{0}.cpp'.format(fname)

    with open(path.join(prefix, thn), 'w') as th:
        th.write(str(timerH(searchList = [{'func' : f}])))

    with open(path.join(prefix, tcn), 'w') as tc:
        tc.write(str(timerC(searchList = [{'func' : f}])))
