from block import block
from block_avx import block_avx

if __name__ == '__main__':
    blen = [1,2,4,8,16]
    for i in blen:
        for j in blen:
            of = 'code/block{0}x{1}_swap.hpp'.format(i, j)
            ofa = 'code/block{0}x{1}_swap_avx.hpp'.format(i, j)
            with open(of, 'w') as out:
                out.write(str(block(searchList=[{'N':i, 'M':j}])))
            if j%4 == 0:
                with open(ofa, 'w') as out:
                    out.write(str(block_avx(searchList=[{'N':i, 'M':j}])))
