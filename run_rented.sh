#!/usr/bin/env bash
# =====================================================================
#  Twenty hours on the rented box (72 vCPU / 36 physical cores, 256 GB)
#
#  Run from an MSYS2 MINGW64 shell.  Each part is independent; if one is
#  cut short its output is still usable, so run them in this order and
#  stop wherever the budget runs out.
#
#  Everything appends, so nothing is ever overwritten.
# =====================================================================
set -u
T=36                     # physical cores, not the 72 logical ones
cd /c/primes || exit 1

# ---------------------------------------------------------------------
# 0.  BUILD  (about a minute)
# ---------------------------------------------------------------------
for f in family_sweep4 mers hunt2; do
  gcc -O3 -march=native -fopenmp -static $f.c -lgmp -o $f.exe || exit 1
done
echo "build ok"

# ---------------------------------------------------------------------
# 1.  CALIBRATE THIS MACHINE FIRST  (about 15 minutes)
#
#     Everything below is sized from the laptop's numbers.  This box has
#     multi-channel memory, so its parallel speedup should be far better
#     than the laptop's 2.3x -- but measure it rather than assume it.
#     Compare the "tested" counts: 36 threads should be roughly four
#     times 8 threads.  If it is much less, cut the target sizes in
#     part 3 by about a quarter for each halving of throughput.
# ---------------------------------------------------------------------
./hunt2.exe -scheme A -digits 8000 -budget 300 -t 8  -B 100000000
./hunt2.exe -scheme A -digits 8000 -budget 300 -t $T -B 100000000
mv hunt2.csv calibration.csv
mv hunt2_primes.txt calibration_primes.txt 2>/dev/null
echo "=== calibration done -- check the two 'tested' numbers before going on ==="

# ---------------------------------------------------------------------
# 2.  THE POWER-OF-TWO EFFECT           (about 10 hours)
#
#     The one positive claim in the paper: offsets with Pr = 1 produced
#     23% more primes than offsets with Pr prime, on a sample of only 15
#     offsets (permutation p = 0.0019).  Raising C from 50000 to 800000
#     brings the count of Pr = 1 offsets to 19 and multiplies the primes
#     found by roughly sixteen, which is what decides the question.
#
#     -nosave skips the multi-gigabyte text dump; the CSV keeps p, sign
#     and c, and every number can be rebuilt from those.
# ---------------------------------------------------------------------
./family_sweep4.exe -N 700 -e 1 -C 800000 -t $T -B 200000000 -nosave -q
mv sweep4.csv      pow2_e1_N700_C800k.csv
mv sweep4_time.csv pow2_e1_N700_C800k_time.csv

python features.py pow2_e1_N700_C800k.csv --cmax 800000 --lo 100 --hi 2003 \
       > pow2_report.txt
echo "=== part 2 done -- see pow2_report.txt ==="

# ---------------------------------------------------------------------
# 3.  THE THREE STRATEGIES AT LARGER SIZES   (about 6 hours)
#
#     The equal-size comparison currently runs from 473 to 1829 digits
#     and says A, B and C are indistinguishable.  This extends it to
#     about 9000 digits.  If the verdict survives an order of magnitude
#     in size, the paper can state it without hedging.
#
#     J is doubled to 120 so each size contributes twice the candidates.
# ---------------------------------------------------------------------
for b in $(seq 3000 50 30000); do
  for s in A B C; do
    ./mers.exe -bits $b -scheme $s -J 120 -t $T -B 10000000 -s $b
  done
done
mv mers.csv mers_3k_30k_bits.csv
echo "=== part 3 done ==="

# ---------------------------------------------------------------------
# 4.  ONE RECORD RUN                     (about 3 hours)
#
#     No -digits here on purpose: hunt2 times a real modular
#     exponentiation on this machine and picks the largest size with a
#     realistic chance inside the budget.  On the laptop that logic gave
#     12000 digits in five hours; here it should choose somewhere around
#     40000 to 50000.
# ---------------------------------------------------------------------
./hunt2.exe -scheme A -budget 10800 -t $T -B 1000000000
mv hunt2.csv        record_run.csv
mv hunt2_primes.txt record_primes.txt

echo
echo "all done.  Files to keep:"
echo "  calibration.csv                 machine speed, both thread counts"
echo "  pow2_e1_N700_C800k.csv          the power-of-two dataset"
echo "  pow2_report.txt                 the feature analysis"
echo "  mers_3k_30k_bits.csv            three-strategy comparison"
echo "  record_run.csv, record_primes.txt"
