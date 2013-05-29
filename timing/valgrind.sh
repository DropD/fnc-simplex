
#~ LP='../problems/gen_light/lpb_00050_1.dlp'
#~ LP='../problems/gen_light/lpb_00200_1.dlp'
LP='../problems/gen_light/lpb_01000_1.dlp'
#~ LP='../problems/gen_upper/lpb_04000_1.dlp'

COLLECT='Simplex*::solve*'

curdir=`pwd`
cd ../src
make debug
cd "$curdir"

valgrind \
   --tool=callgrind --callgrind-out-file=callgrind.out \
   --toggle-collect="$COLLECT" \
   --dump-before="$COLLECT" \
   --dump-instr=yes \
   --cache-sim=yes \
   ../bin/main-dbg $LP
#~ callgrind_annotate --inclusive=yes --auto=yes --tree=both callgrind.out
callgrind_annotate --inclusive=yes --tree=both callgrind.out
rm callgrind.out

valgrind \
   --tool=cachegrind --cachegrind-out-file=cachegrind.out \
   --cache-sim=yes \
   ../bin/main-dbg $LP
cg_annotate --auto=yes cachegrind.out
rm cachegrind.out
