#!/usr/bin/env bash
# =====================================================================
#  One ten-hour run, strategy B (every exponent 2, so M is the squared
#  primorial), on the laptop.  Target: about 20500 digits, expected
#  three primes.
#
#  Run from MSYS2 MINGW64, in the folder holding hunt2.c
# =====================================================================
set -u
cd /c/primes || exit 1

gcc -O3 -march=native -fopenmp -static hunt2.c -lgmp -o hunt2.exe || exit 1

# --- 10 minutes first: confirm the machine still calibrates the same ---
#     Look at "one 2000-digit test".  On this laptop it has been 0.060 s.
#     If it comes out much slower (a hot room, another program running),
#     drop -digits below by a thousand for every 10% of slowdown.
./hunt2.exe -scheme B -digits 6000 -budget 300 -t 4 -B 100000000
mv hunt2.csv       precheck.csv
mv hunt2_primes.txt precheck_primes.txt 2>/dev/null
echo "=== precheck done ==="

# --- the real run ---
#
#  -digits 20500  fixed on purpose, so the run is reproducible: the
#                 automatic sizing varies a little with the calibration
#                 and would shift the whole candidate set between runs.
#  -B 1000000000  a deeper sieve than usual.  At this size one test costs
#                 about 11 core-seconds while one reduction costs microsec-
#                 onds, so sieving stays profitable far past the usual
#                 bound; going from 10^8 to 10^9 removes roughly a ninth
#                 of the remaining tests.  It needs about 1.2 GB for a
#                 minute at startup.
#  -t 4           on this machine four threads give about 1.4x at these
#                 sizes, not 4x: single-channel memory is the limit, and
#                 the FFT working set no longer fits in cache.  Two
#                 threads are within a few per cent if the fans bother you.
#
./hunt2.exe -scheme B -digits 20500 -budget 36000 -t 4 -B 1000000000

mv hunt2.csv       B_20500.csv
mv hunt2_primes.txt B_20500_primes.txt

echo
echo "done.  The last lines of the run tell you how many of the survivors"
echo "were actually tested; if candidates remain, the printed -from value"
echo "continues from where this stopped instead of repeating the same ones."
