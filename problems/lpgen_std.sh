
LPGEN='python ../lp_gen/build/gen_lp.py'
LPDIR='gen_std'

mkdir -p $LPDIR

run_gen() {
    echo "generate $1 lps of size $2"
    for i in $(seq 1 $1); do
        $LPGEN $2 $LPDIR/lpb_`printf '%05d' $2`_"$i" ;
    done
}

run_gen 20 10
run_gen 20 20
run_gen 20 30
run_gen 20 50
run_gen 10 80
run_gen 10 100
run_gen 10 150
run_gen 10 200
run_gen 5 300
run_gen 5 400
run_gen 5 500
run_gen 3 600
run_gen 3 800
run_gen 3 1000
