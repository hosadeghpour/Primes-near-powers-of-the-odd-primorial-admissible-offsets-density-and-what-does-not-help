/*
 * indices.c  --  three arithmetic indices for every n up to a limit
 *
 *   omega(n)  = number of distinct prime factors            (OEIS A001221)
 *   Omega(n)  = number of prime factors with multiplicity   (OEIS A001222)
 *   D(n)      = omega(n) + ( 1 - omega(n)/Omega(n) )
 *   W(n)      = product over distinct primes q | n of q/(q-1)
 *
 * D is the combined index: a prime gives 1, p^2 gives 1.5, p*q gives 2,
 * p^2*q gives 2.33, and so on -- SMALL means prime-like.  Two properties are
 * worth knowing before using it.
 *
 *   1. D is squarefree-integral: if n is squarefree then Omega = omega and
 *      D(n) = omega(n) exactly; otherwise omega < D < omega + 1.  So
 *      floor(D) recovers omega, and Omega = omega / (omega + 1 - D) recovers
 *      Omega.  D therefore carries exactly the information in the pair
 *      (omega, Omega) -- no more, no less.  It is a convenient single number,
 *      not a new invariant.
 *
 *   2. n = 1 has omega = Omega = 0 and the quotient is undefined.  This
 *      program reports D(1) = 0 by convention; change ONE line below if you
 *      prefer something else.
 *
 * W is the Mertens boost, and it is the index that matters if the question is
 * "does n sit where primes are easy to find?".  D counts every prime factor
 * alike; W does not, and it is right not to: dividing by 2 is worth a factor
 * of 2, by 3 a factor of 1.5, by 1000003 a factor of 1.000001.  A number
 * divisible by 2*3*5*7 and a number divisible by 101*103*107*109 have the same
 * omega and the same D, but W = 4.375 against W = 1.04.  It is W that makes
 * the primorial construction work: for n = p# one gets W ~ e^gamma * ln p.
 *
 * A warning that came out of measuring it.  W(n) does NOT predict how many
 * primes lie in an interval around n: over 1.7 million values of n with a
 * window of +-50, the correlation is r = -0.015.  The reason is that if q | n
 * then among the offsets k the ones killed by q are exactly those with q | k,
 * which is the same proportion 1/q as for a random n.  Divisibility of n does
 * not create protection, it only decides which offsets receive it.  What W
 * does raise is the density among offsets COPRIME to n -- and the number of
 * such offsets falls by exactly the compensating factor.
 *
 * All arrays are filled by sieving rather than by factoring each n, so the
 * whole run is O(N log log N) and one pass over memory per prime.
 *
 * Build:
 *   gcc -O3 -march=native indices.c -o indices        (no libraries needed)
 *
 * Run:
 *   ./indices -N 1000000                    -> indices.csv, one row per n
 *   ./indices -N 1000000 -o out.csv
 *   ./indices -N 100000000 --summary        -> distribution only, no file
 *   ./indices -N 1000000 --min 2.9 --max 3.0  -> only rows in that D range
 *
 * Memory: 2 bytes per n, so 2 MB at N = 10^6 and 2 GB at N = 10^9.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static double index_D(unsigned w, unsigned O)
{
    if (O == 0) return 0.0;                 /* the n = 1 convention */
    return (double) w + (1.0 - (double) w / (double) O);
}

int main(int argc, char **argv)
{
    long N = 1000000, n;
    const char *outfile = "indices.csv";
    int summary = 0, i;
    double dmin = -1e30, dmax = 1e30;
    uint8_t *w, *O;
    double *W;
    FILE *f = NULL;

    for (i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-N")  && i + 1 < argc) N = atol(argv[++i]);
        else if (!strcmp(argv[i], "-o")  && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "--min") && i + 1 < argc) dmin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max") && i + 1 < argc) dmax = atof(argv[++i]);
        else if (!strcmp(argv[i], "--summary")) summary = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (N < 1) { fprintf(stderr, "N must be at least 1\n"); return 1; }

    w = calloc((size_t) N + 1, 1);
    O = calloc((size_t) N + 1, 1);
    W = malloc(((size_t) N + 1) * sizeof(double));
    if (!w || !O || !W) { fprintf(stderr, "out of memory at N = %ld\n", N); return 1; }
    for (n = 0; n <= N; n++) W[n] = 1.0;

    /* omega: every n with w[n] still 0 when we reach it is prime, and each
       prime marks all of its multiples exactly once                        */
    for (n = 2; n <= N; n++) {
        if (w[n] == 0) {                        /* n is prime */
            long m;
            double boost = (double) n / (double) (n - 1);
            for (m = n; m <= N; m += n) { w[m]++; W[m] *= boost; }
            /* Omega: one pass per prime POWER, so p^k contributes k        */
            {
                long q = n;
                while (q <= N / n || q <= N) {
                    for (m = q; m <= N; m += q) O[m]++;
                    if (q > N / n) break;       /* next power would overflow */
                    q *= n;
                }
            }
        }
    }

    if (!summary) {
        f = fopen(outfile, "w");
        if (!f) { fprintf(stderr, "cannot write %s\n", outfile); return 1; }
        fprintf(f, "n,omega,Omega,D,W\n");
    }

    {
        long counts[64];
        long shown = 0;
        double dsum = 0.0;
        memset(counts, 0, sizeof counts);
        for (n = 1; n <= N; n++) {
            double D = index_D(w[n], O[n]);
            if (w[n] < 63) counts[w[n]]++;
            dsum += D;
            if (!summary && D >= dmin && D <= dmax) {
                fprintf(f, "%ld,%u,%u,%.6f,%.6f\n", n, w[n], O[n], D, W[n]);
                shown++;
            }
        }
        if (f) { fclose(f); printf("wrote %ld rows to %s\n", shown, outfile); }
        printf("\nN = %ld\n", N);
        printf("mean D = %.6f\n", dsum / (double) N);
        printf("\n%8s %14s %10s\n", "omega", "count", "share");
        for (i = 0; i < 20; i++)
            if (counts[i])
                printf("%8d %14ld %9.4f%%\n", i, counts[i],
                       100.0 * (double) counts[i] / (double) N);
        printf("\n(omega = 1 covers the primes and the prime powers;\n"
               " D separates them: primes give exactly 1, p^k gives 2 - 1/k.)\n");
    }
    free(w); free(O); free(W);
    return 0;
}
