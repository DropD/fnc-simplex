
LPGEN='python ../lp_gen/build/gen_lp.py'

#~ echo generate 3 lps of size 4
#~ for i in {1..3}; do $LPGEN 4 gen/lpb_4_$"$i" ; done

echo generate 100 lps of size 10
for i in {1..100}; do $LPGEN 10 gen/lpb_0010_$"$i" ; done

echo generate 100 lps of size 20
for i in {1..100}; do $LPGEN 20 gen/lpb_0020_$"$i" ; done

echo generate 50 lps of size 50
for i in {1..50}; do $LPGEN 50 gen/lpb_0050_$"$i" ; done

echo generate 20 lps of size 100
for i in {1..20}; do $LPGEN 100 gen/lpb_0100_$"$i" ; done

echo generate 5 lps of size 500
for i in {1..5}; do $LPGEN 500 gen/lpb_0500_$"$i" ; done

echo generate 3 lp of size 1000
for i in {1..3}; do $LPGEN 1000 gen/lpb_1000_$"$i" ; done
