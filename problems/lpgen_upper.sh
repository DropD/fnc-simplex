
LPGEN='python ../lp_gen/build/gen_lp.py'
LPDIR='gen_upper'

mkdir -p $LPDIR

run_gen() {
    echo "generate $1 lps of size $2"
    for i in $(seq 1 $1); do
        $LPGEN 10 $LPDIR/lpb_`printf '%04d' $2`_"$i" ;
    done
}

run_gen 3 500
run_gen 3 600
run_gen 3 1000
run_gen 3 2000
run_gen 3 5000
run_gen 3 10000
