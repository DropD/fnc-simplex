
LPGEN='python ../lp_gen/build/gen_lp.py'
LPDIR='gen_heavy'

mkdir -p $LPDIR

run_gen() {
    echo "generate $1 lps of size $2"
    for i in $(seq 1 $1); do
        $LPGEN $2 $LPDIR/lpb_`printf '%05d' $2`_"$i" ;
    done
}

run_gen 100 10
run_gen 100 20
run_gen 50 30
run_gen 50 50
run_gen 30 80
run_gen 20 100
run_gen 20 150
run_gen 20 200
run_gen 10 300
run_gen 10 400
run_gen 5 500
run_gen 5 600
run_gen 3 1000
run_gen 3 2000
