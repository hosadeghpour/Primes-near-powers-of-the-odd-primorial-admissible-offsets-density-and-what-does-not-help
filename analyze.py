#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze.py -- reads the CSV written by family_sweep2 and produces the
tables the paper needs.  Standard library only: no numpy, no sympy.

    python analyze.py sweep2.csv --cmax 50000

Reports
-------
1. sign balance                     (is + or - more productive?)
2. primes by power of two           (is some k special?)
3. powers of two vs other offsets   (the 2.5-sigma effect, tested fairly)
4. the Fortunate-type function      F(p), and a check of the conjecture
5. observed vs predicted counts     (Bateman-Horn, by family)
"""

import argparse
import csv
import math
import statistics
from collections import Counter

EULER = math.exp(0.5772156649015329)


def sieve(limit):
    f = bytearray([1]) * (limit + 1)
    f[0] = f[1] = 0
    i = 2
    while i * i <= limit:
        if f[i]:
            f[i * i::i] = bytearray(len(f[i * i::i]))
        i += 1
    return [i for i in range(limit + 1) if f[i]]


def odd_part(c):
    while c % 2 == 0:
        c //= 2
    return c


def is_pow2(c):
    return c & (c - 1) == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csvfile")
    ap.add_argument("--cmax", type=int, default=20000,
                    help="the -C value used in the sweep")
    args = ap.parse_args()

    rows = []
    with open(args.csvfile, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append((int(r["p"]), int(r["digits"]), r["sign"], int(r["c"])))
    if not rows:
        raise SystemExit("no rows")

    C = args.cmax
    prime_list = sieve(C)
    prime_set = set(prime_list)
    fams = sorted({p for p, _, _, _ in rows})
    digits = {p: d for p, d, _, _ in rows}
    print(f"{len(rows)} primes, {len(fams)} families, p = {fams[0]}..{fams[-1]}, "
          f"c <= {C}\n")

    # ---- 1. sign balance -------------------------------------------------
    minus = sum(1 for _, _, s, _ in rows if s == "-")
    plus = len(rows) - minus
    sd = math.sqrt(len(rows))
    print("1. SIGN")
    print(f"   minus {minus}   plus {plus}   difference {minus - plus:+d}"
          f"   (1 sigma = {sd:.0f}, so {(minus - plus) / sd:+.1f} sigma)\n")

    # ---- 2. by power of two ---------------------------------------------
    print("2. PURE POWERS OF TWO, BY EXPONENT")
    bk = Counter((int(round(math.log2(c))), s) for _, _, s, c in rows if is_pow2(c))
    kmax = int(math.log2(C))
    print(f"   {'k':>4} {'minus':>7} {'plus':>7} {'total':>7}")
    tot2 = 0
    for k in range(1, kmax + 1):
        a, b = bk.get((k, "-"), 0), bk.get((k, "+"), 0)
        tot2 += a + b
        print(f"   {k:>4} {a:>7} {b:>7} {a + b:>7}")
    print(f"   total {tot2}  ({100 * tot2 / len(rows):.1f}% of all primes found)\n")

    # ---- 3. powers of two vs the rest, on equal footing ------------------
    # An offset 2^j*m is admissible only for families with p < m, so compare
    # only over a family range where both groups are admissible everywhere.
    print("3. POWERS OF TWO VS OTHER OFFSETS (same families, same eligibility)")
    print(f"   {'p range':>16} {'fam':>5} {'pow2':>8} {'other':>8} {'diff':>8} {'sigma':>7}")
    for hi in (300, 601, 1009, 2003, 5003):
        fam = [p for p in fams if p <= hi]
        if len(fam) < 20:
            continue
        famset = set(fam)
        got = Counter(c for p, _, _, c in rows if p in famset)
        pw, oth = [], []
        for c in range(2, C + 1, 2):
            o = odd_part(c)
            if o == 1:
                pw.append(got.get(c, 0))
            elif o in prime_set and o > hi:
                oth.append(got.get(c, 0))
        if len(pw) < 2 or len(oth) < 2:
            continue
        mp, mo = statistics.mean(pw), statistics.mean(oth)
        se = math.sqrt(statistics.pvariance(pw) / len(pw) +
                       statistics.pvariance(oth) / len(oth))
        print(f"   {fam[0]:>7}..{hi:<7} {len(fam):>5} {mp:>8.2f} {mo:>8.2f} "
              f"{mp - mo:>+8.2f} {(mp - mo) / se if se else 0:>7.1f}")
    print("   (mean number of primes contributed by one offset of each kind;\n"
          "    three sigma or more in every row would make this a real effect)\n")

    # ---- 4. the Fortunate-type function ----------------------------------
    print("4. FORTUNATE-TYPE FUNCTION  F(p) = least admissible c giving a prime")
    least = {}
    for p, _, s, c in rows:
        k = (p, s)
        if k not in least or c < least[k]:
            least[k] = c
    bad_odd, over_p2, censored = [], [], []
    for p in fams:
        for s in "-+":
            c = least.get((p, s))
            if c is None:
                censored.append((p, s))
                continue
            o = odd_part(c)
            if o > 1 and o not in prime_set:
                bad_odd.append((p, s, c, o))
            if c >= p * p:
                over_p2.append((p, s, c))
    print(f"   values computed: {len(least)}")
    print(f"   odd part neither 1 nor prime: {len(bad_odd)}  {bad_odd[:4]}")
    print(f"   values reaching p^2         : {len(over_p2)}  {over_p2[:4]}")
    print(f"   families with no prime found: {len(censored)}  {censored[:6]}")
    if least:
        mx = max(least.items(), key=lambda kv: kv[1])
        print(f"   largest F: {mx[1]} at p={mx[0][0]} sign {mx[0][1]} "
              f"(p^2 = {mx[0][0] ** 2})")
    print("   first values:")
    print(f"   {'p':>6} {'F-':>8} {'F+':>8}")
    for p in fams[:16]:
        print(f"   {p:>6} {str(least.get((p, '-'), '-')):>8} "
              f"{str(least.get((p, '+'), '-')):>8}")
    print()

    # ---- 5. observed vs predicted ---------------------------------------
    print("5. OBSERVED / PREDICTED  (Bateman-Horn: density = e^gamma*ln p / ln N)")
    n_adm = {}
    for p in fams:
        n = 0
        for c in range(2, C + 1, 2):
            o = odd_part(c)
            if o == 1 or (o in prime_set and o > p):
                n += 1
        n_adm[p] = n
    got_p = Counter(p for p, _, _, _ in rows)
    ratios = []
    for p in fams:
        lnN = digits[p] * math.log(10)
        exp = 2 * n_adm[p] * EULER * math.log(p) / lnN   # two signs
        if exp > 3:
            ratios.append(got_p[p] / exp)
    if ratios:
        print(f"   families used: {len(ratios)}")
        print(f"   mean {statistics.mean(ratios):.3f}   "
              f"median {statistics.median(ratios):.3f}   "
              f"sd {statistics.pstdev(ratios):.3f}")
        print("   (a mean near 1.00 confirms the heuristic for this family)")


if __name__ == "__main__":
    main()
