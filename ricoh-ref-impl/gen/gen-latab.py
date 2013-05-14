from table import table
import sys
from os import path
import ethpy as ep

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: python gen-latab.py <datafile> <outfile>'
        sys.exit()

    dfile, ofile = sys.argv[1:3]

    data = [i[0] for i in ep.read(dfile)]

    with open(ofile, 'w') as tex:
        tex.write(str(table(searchList = [{'data' : data}])))
