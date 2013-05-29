
LPGEN='python ../lp_gen/build/gen_lp.py'
LPDIR='gen_light'

mkdir -p $LPDIR

run_gen() {
    echo "generate $1 lps of size $2"
    for i in $(seq 1 $1); do
        $LPGEN $2 $LPDIR/lpb_`printf '%05d' $2`_"$i" ;
        echo -n '.'
    done
    echo ''
}

run_gen 5 10
run_gen 5 20
run_gen 5 30
run_gen 5 50
run_gen 5 80
run_gen 3 100
run_gen 3 150
run_gen 3 200
run_gen 3 300
run_gen 3 400
run_gen 2 500
run_gen 2 600
run_gen 1 800
run_gen 1 1000
