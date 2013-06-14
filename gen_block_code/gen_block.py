from block import block
from block_swap import block_swap
from block_avx import block_avx
from block_swap_avx import block_swap_avx
import os

if __name__ == '__main__':
    blen = [1,2,4,8,16]
    dirnames = [os.path.join('code',i) for i in ['block', 'block_swap', 'block_avx', 'block_swap_avx']]
    for dn in dirnames:
        if not os.path.exists(dn):
            os.mkdir(dn)

    for i in blen:
        for j in blen:
            b = 'code/block/block{0}x{1}.hpp'.format(i, j)
            bs = 'code/block_swap/block{0}x{1}_swap.hpp'.format(i, j)
            with open(b, 'w') as out:
                out.write(str(block(searchList=[{'N':i, 'M':j}])))
            with open(bs, 'w') as out:
                out.write(str(block_swap(searchList=[{'N':i, 'M':j}])))
            if j%4 == 0:
                ba = 'code/block_avx/block{0}x{1}_avx.hpp'.format(i, j)
                bsa = 'code/block_swap_avx/block{0}x{1}_swap_avx.hpp'.format(i, j)
                with open(ba, 'w') as out:
                    out.write(str(block_avx(searchList=[{'N':i, 'M':j}])))
                with open(bsa, 'w') as out:
                    out.write(str(block_swap_avx(searchList=[{'N':i, 'M':j}])))

    hb = 'code/block.hpp'
    hbs = 'code/block_swap.hpp'
    hba = 'code/block_avx.hpp'
    hbsa = 'code/block_swap_avx.hpp'

    with open(hb, 'w') as hbf, open(hbs, 'w') as hbsf:
        for i in blen:
            for j in blen:
                b = 'block/block{0}x{1}.hpp'.format(i, j)
                bs = 'block_swap/block{0}x{1}_swap.hpp'.format(i, j)
                hbf.write('#include "{0}"\n'.format(b))
                hbsf.write('#include "{0}"\n'.format(bs))
 
    with open(hba, 'w') as hbaf, open(hbsa, 'w') as hbsaf:
        for i in blen:
            for j in blen:
                if j%4 == 0:
                    ba = 'block_avx/block{0}x{1}_avx.hpp'.format(i, j)
                    bsa = 'block_swap_avx/block{0}x{1}_swap_avx.hpp'.format(i, j)
                    hbaf.write('#include "{0}"\n'.format(ba))
                    hbsaf.write('#include "{0}"\n'.format(bsa))
