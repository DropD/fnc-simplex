
LPGEN='python ../lp_gen/build/gen_lp.py'

#~ echo generate 3 lps of size 4
#~ for i in {1..3}; do $LPGEN 4 gen_verylight/lpb_4_$"$i" ; done

echo generate 10 lps of size 10
for i in {1..10}; do $LPGEN 10 gen_verylight/lpb_0010_$"$i" ; done

echo generate 10 lps of size 20
for i in {1..10}; do $LPGEN 20 gen_verylight/lpb_0020_$"$i" ; done

echo generate 10 lps of size 30
for i in {1..10}; do $LPGEN 30 gen_verylight/lpb_0030_$"$i" ; done

echo generate 10 lps of size 50
for i in {1..10}; do $LPGEN 50 gen_verylight/lpb_0050_$"$i" ; done

echo generate 10 lps of size 80
for i in {1..10}; do $LPGEN 80 gen_verylight/lpb_0080_$"$i" ; done

echo generate 5 lps of size 100
for i in {1..5}; do $LPGEN 100 gen_verylight/lpb_0100_$"$i" ; done

echo generate 5 lps of size 150
for i in {1..5}; do $LPGEN 150 gen_verylight/lpb_0150_$"$i" ; done

echo generate 5 lps of size 200
for i in {1..5}; do $LPGEN 200 gen_verylight/lpb_0200_$"$i" ; done

echo generate 5 lps of size 300
for i in {1..5}; do $LPGEN 300 gen_verylight/lpb_0300_$"$i" ; done

echo generate 5 lps of size 400
for i in {1..5}; do $LPGEN 400 gen_verylight/lpb_0400_$"$i" ; done

echo generate 3 lps of size 500
for i in {1..3}; do $LPGEN 500 gen_verylight/lpb_0500_$"$i" ; done

echo generate 2 lps of size 600
for i in {1..2}; do $LPGEN 600 gen_verylight/lpb_0600_$"$i" ; done

echo generate 1 lps of size 1000
for i in {1..1}; do $LPGEN 1000 gen_verylight/lpb_1000_$"$i" ; done
