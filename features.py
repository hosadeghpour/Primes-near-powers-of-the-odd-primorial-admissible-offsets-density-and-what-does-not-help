#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
features.py -- does any property of the offset predict whether it produces a prime?

    python features.py e1_N700_C50k.csv --cmax 50000 --lo 100 --hi 2003

Reads a sweep CSV with columns p,digits,sign,c and pools every admissible
offset across the chosen family range.  For each offset it computes a battery
of features of the odd part m (the offset is 2^j * m, and by the admissibility
theorem m is 1 or a prime above p), then reports, for each feature, the
correlation with the number of primes that offset produced.

The comparison is only fair over a family range where every offset is equally
admissible, so --lo/--hi restrict to families with lo <= p <= hi, and only
offsets with m > hi are pooled.  Powers of two (m = 1) are reported separately
because they have no odd neighbourhood.

Adding your own feature: write a function of m returning a number, and put it
in the FEATURES dict.  Everything else is automatic.

Standard library only.
"""

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict


def sieve_flags(n):
    f = bytearray([1]) * (n + 1)
    f[0] = f[1] = 0
    i = 2
    while i * i <= n:
        if f[i]:
            f[i * i::i] = bytearray(len(f[i * i::i]))
        i += 1
    return f


def spf_table(n):
    spf = list(range(n + 1))
    i = 2
    while i * i <= n:
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


class Arith:
    """omega, Omega and the primality degree, from one spf table."""

    def __init__(self, limit):
        self.spf = spf_table(limit)
        self.isp = sieve_flags(limit)
        self.limit = limit

    def factor(self, n):
        f = {}
        while n > 1:
            q = self.spf[n]
            while n % q == 0:
                n //= q
                f[q] = f.get(q, 0) + 1
        return f

    def omega(self, n):
        return len(self.factor(n)) if n > 1 else 0

    def Omega(self, n):
        return sum(self.factor(n).values()) if n > 1 else 0

    def degree(self, n):
        """The user's primality degree: omega + (1 - omega/Omega).

        A prime gives 1, p^2 gives 1.5, pq gives 2, p^2*q gives 2.33 ...
        so SMALL means prime-like.  The value determines (omega, Omega)
        uniquely: floor recovers omega (except for squarefree n, where the
        value is the integer omega itself), and Omega follows.  It is a
        re-encoding of the classical pair, not new information.
        """
        if n <= 1:
            return 0.0
        w = self.omega(n)
        O = self.Omega(n)
        return w + (1.0 - w / O)

    def prev_prime(self, n):
        k = n - 1
        while k > 1 and not self.isp[k]:
            k -= 1
        return k

    def next_prime(self, n):
        k = n + 1
        while k <= self.limit and not self.isp[k]:
            k += 1
        return k


NEIGHBOURS = [2, 4, 6, 8, 10, 12, 18]


def build_features(A):
    """Every feature is a function of the odd part m."""
    D = [d for x in NEIGHBOURS for d in (-x, x)]

    def n_prime_nbrs(m):
        return sum(1 for d in D if 1 < m + d <= A.limit and A.isp[m + d])

    def mean_degree(m):
        v = [A.degree(m + d) for d in D if 1 < m + d <= A.limit]
        return statistics.mean(v) if v else 0.0

    def min_degree(m):
        v = [A.degree(m + d) for d in D if 1 < m + d <= A.limit]
        return min(v) if v else 0.0

    def gap_prev(m):
        return m - A.prev_prime(m)

    def gap_next(m):
        return A.next_prime(m) - m

    def gap_total(m):
        return A.next_prime(m) - A.prev_prime(m)

    return {
        "prime_neighbours": n_prime_nbrs,       # how prime-dense around m
        "mean_degree":      mean_degree,        # user's degree, averaged
        "min_degree":       min_degree,
        "gap_prev":         gap_prev,
        "gap_next":         gap_next,
        "gap_total":        gap_total,          # isolation of m
        "twin_up":          lambda m: 1 if m + 2 <= A.limit and A.isp[m + 2] else 0,
        "twin_down":        lambda m: 1 if m > 2 and A.isp[m - 2] else 0,
        "sophie_germain":   lambda m: 1 if 2 * m + 1 <= A.limit and A.isp[2 * m + 1] else 0,
        "omega_m_minus_1":  lambda m: A.omega(m - 1),
        "omega_m_plus_1":   lambda m: A.omega(m + 1),
        "m_mod_3":          lambda m: m % 3,
        "m_mod_4":          lambda m: m % 4,
        "m_mod_8":          lambda m: m % 8,
    }


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return 0.0, 0.0
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sx = sum((a - mx) ** 2 for a in xs)
    sy = sum((b - my) ** 2 for b in ys)
    if sx == 0 or sy == 0:
        return 0.0, 0.0
    r = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / math.sqrt(sx * sy)
    return r, r * math.sqrt(n - 2) / math.sqrt(max(1e-12, 1 - r * r))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csvfile")
    ap.add_argument("--cmax", type=int, default=50000)
    ap.add_argument("--lo", type=int, default=100)
    ap.add_argument("--hi", type=int, default=2003)
    args = ap.parse_args()

    C, LO, HI = args.cmax, args.lo, args.hi
    A = Arith(C + 200)

    rows = [(int(r["p"]), int(r["c"])) for r in csv.DictReader(open(args.csvfile))]
    fams = {p for p, _ in rows if LO <= p <= HI}
    got = Counter(c for p, c in rows if p in fams)
    print("families %d..%d used: %d" % (LO, HI, len(fams)))

    pow2, offs = [], []
    for c in range(2, C + 1, 2):
        m = c
        while m % 2 == 0:
            m //= 2
        if m == 1:
            pow2.append(got.get(c, 0))
        elif A.isp[m] and m > HI:
            offs.append((m, got.get(c, 0)))
    if not offs:
        raise SystemExit("no offsets pooled -- lower --hi or raise --cmax")

    ys = [y for _, y in offs]
    base = statistics.mean(ys)
    print("offsets pooled: %d   mean primes each: %.3f  (sd %.3f)"
          % (len(offs), base, statistics.pstdev(ys)))
    if pow2:
        print("pure powers of two: %d offsets, mean %.3f  (%+.1f%%)"
              % (len(pow2), statistics.mean(pow2),
                 100 * (statistics.mean(pow2) - base) / base))
    print()

    feats = build_features(A)
    print("%-18s %9s %8s %s" % ("feature", "r", "t", "group means"))
    results = []
    for name, fn in feats.items():
        xs = [fn(m) for m, _ in offs]
        r, t = pearson(xs, ys)
        buckets = defaultdict(list)
        for x, y in zip(xs, ys):
            buckets[round(x, 1)].append(y)
        shown = sorted(b for b in buckets if len(buckets[b]) >= 30)[:6]
        desc = "  ".join("%g:%.2f" % (b, statistics.mean(buckets[b])) for b in shown)
        results.append((abs(t), name, r, t, desc))
    results.sort(reverse=True)
    for _, name, r, t, desc in results:
        print("%-18s %+9.4f %8.2f  %s" % (name, r, t, desc))

    print()
    n = len(offs)
    print("With %d offsets, |t| above about 2 is worth a second look, and with "
          "%d features\ntested at once about %.1f is expected from chance alone."
          % (n, len(feats), statistics.NormalDist().inv_cdf(1 - 0.5 / len(feats))))


if __name__ == "__main__":
    main()
