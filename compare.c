/*
 * compare.c  --  which exponent strategy is best?
 *
 *      M(e) = product over odd primes q <= p of q^(e_q),   e_q >= 0
 *      N    = | M(e) -/+ 2^k -/+ l | ,   1 <= k <= 20,
 *             l in {0,2,4,6,8,10,12,14,16,18,20,30}
 *
 *  A : e_q = 1 for every q          (the plain primorial)
 *  B : e_q = 2 for every q
 *  C : e_q drawn uniformly from {0,1,2}, redrawn for every p,
 *      seeded by (master_seed, p) so the whole run is reproducible
 *
 * Three comparison modes, because the three of them answer different
 * questions and only reporting all three is honest:
 *
 *   -mode digits   every strategy is grown until N has about the same number
 *                  of digits.  This is the fair test of DENSITY: same cost per
 *                  test, so whoever finds more primes wins.
 *   -mode pmax     every strategy uses the same largest prime p.  B and C then
 *                  produce different sizes; this isolates the effect of the
 *                  exponents on the set of primes dividing M.
 *   -mode time     each strategy gets the same wall-clock budget and climbs to
 *                  larger and larger sizes; the answer is the largest prime it
 *                  managed to find.  This is the test that matters if your
 *                  goal is a record.
 *
 * A correctness point that decides the whole comparison: an offset c is
 * inadmissible only when its odd part shares a prime with the SUPPORT of e,
 * i.e. with the primes that actually divide M.  Under strategy C a prime with
 * e_q = 0 is absent from M, so it is allowed to divide N and it no longer
 * protects the candidate.  Testing C against the wrong admissibility rule
 * would make it look better than it is.
 *
 * Build:  gcc -O3 -march=native -fopenmp -static compare.c -lgmp -o compare.exe
 * Run:    ./compare.exe -mode digits -D 2000 -t 4
 *         ./compare.exe -mode pmax   -p 4000 -t 4
 *         ./compare.exe -mode time   -T 600  -t 4
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static const int L_LIST[] = {0,2,4,6,8,10,12,14,16,18,20,30};
static const int N_L = 12;
static int KMAX = 20;
static long OFF;                       /* offsets live in [-OFF, OFF] */

static double now(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return 0.0;
#endif
}

static uint32_t *sieve_primes(uint32_t limit, size_t *cnt)
{
    uint8_t *f = calloc((size_t) limit + 1, 1);
    uint32_t i, j; size_t n = 0; uint32_t *out;
    for (i = 2; (uint64_t) i * i <= limit; i++)
        if (!f[i]) for (j = i * i; j <= limit; j += i) f[j] = 1;
    for (i = 2; i <= limit; i++) if (!f[i]) n++;
    out = malloc(n * sizeof(uint32_t)); n = 0;
    for (i = 2; i <= limit; i++) if (!f[i]) out[n++] = i;
    free(f); *cnt = n; return out;
}

static uint32_t *spf_table(uint32_t limit)
{
    uint32_t *spf = malloc(((size_t) limit + 1) * sizeof(uint32_t));
    uint32_t i, j;
    for (i = 0; i <= limit; i++) spf[i] = i;
    for (i = 2; (uint64_t) i * i <= limit; i++)
        if (spf[i] == i)
            for (j = i * i; j <= limit; j += i) if (spf[j] == j) spf[j] = i;
    return spf;
}

/* deterministic small PRNG so a run can be reproduced from the seed alone */
static uint64_t rng_state;
static void rng_seed(uint64_t s) { rng_state = s * 6364136223846793005ULL + 1; }
static uint32_t rng_next(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17; return (uint32_t) (rng_state >> 32);
}

typedef struct {
    long   candidates, survivors, primes;
    double seconds;
    size_t digits, max_digits;
    uint32_t pmax;
    long   support;
} result_t;

/* one experiment: build M from the exponent vector, sieve, test */
static result_t run_family(int scheme, uint32_t pmax, uint64_t seed,
                           const uint32_t *odd, size_t nodd,
                           const uint32_t *small, size_t nsmall,
                           const uint32_t *spf, uint64_t Bmax, int verbose)
{
    result_t R; mpz_t M, cand;
    uint8_t *mark, *insup;
    size_t i, j;
    long idx;
    int k, li, s1, s2;
    double t0 = now();

    memset(&R, 0, sizeof R);
    R.pmax = pmax;
    mpz_init(M); mpz_init(cand);
    mpz_set_ui(M, 1);
    insup = calloc((size_t) pmax + 2, 1);
    rng_seed(seed);
    for (i = 0; i < nodd && odd[i] <= pmax; i++) {
        int ex = (scheme == 0) ? 1 : (scheme == 1) ? 2 : (int) (rng_next() % 3);
        if (ex == 0) continue;
        { mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, odd[i], ex);
          mpz_mul(M, M, t); mpz_clear(t); }
        insup[odd[i]] = 1; R.support++;
    }
    R.digits = mpz_sizeinbase(M, 10);

    /* ---- distinct offsets, filtered by admissibility w.r.t. supp(e) ---- */
    mark = calloc(2 * OFF + 1, 1);
    for (k = 1; k <= KMAX; k++)
        for (li = 0; li < N_L; li++)
            for (s1 = -1; s1 <= 1; s1 += 2)
                for (s2 = -1; s2 <= 1; s2 += 2) {
                    long v = (long) s1 * (1L << k) + (long) s2 * L_LIST[li];
                    if (v == 0) continue;
                    mark[v + OFF] = 1;
                }
    for (idx = 0; idx <= 2 * OFF; idx++) {
        unsigned long a;
        if (!mark[idx]) continue;
        a = labs(idx - OFF);
        while ((a & 1) == 0) a >>= 1;
        while (a > 1) {                       /* every prime factor of the */
            uint32_t f = spf[a];              /* odd part must be outside  */
            if (f <= pmax && insup[f]) { mark[idx] = 0; break; }
            while (a % f == 0) a /= f;
        }
        if (mark[idx]) R.candidates++;
    }

    /* ---- sieve the offset axis ---------------------------------------- */
    {   uint64_t q_ceiling = 0xFFFFFFFFFFFFFFFFull;
        if (mpz_fits_ulong_p(M)) {
            unsigned long qv = mpz_get_ui(M);
            q_ceiling = (qv > (unsigned long) OFF + 2) ? qv - OFF : 0;
        }
        for (j = 0; j < nsmall; j++) {
            uint32_t q = small[j];
            long r, w;
            if (q > Bmax || q >= q_ceiling) break;
            if (q <= pmax && insup[q]) continue;   /* cannot divide N      */
            r = (long) mpz_fdiv_ui(M, q);
            w = r - (long) q * ((r + OFF) / (long) q);
            while (w < -OFF) w += q;
            for (; w <= OFF; w += q) if (mark[w + OFF]) mark[w + OFF] = 0;
        }
    }

    /* ---- test what is left -------------------------------------------- */
#pragma omp parallel
    {
        mpz_t c2; mpz_init(c2);
#pragma omp for schedule(dynamic, 4)
        for (idx = 0; idx <= 2 * OFF; idx++) {
            long v;
            if (!mark[idx]) continue;
            v = idx - OFF;
            if (v < 0) mpz_add_ui(c2, M, (unsigned long) (-v));
            else       mpz_sub_ui(c2, M, (unsigned long) v);
            mpz_abs(c2, c2);
#pragma omp atomic
            R.survivors++;
            if (mpz_cmp_ui(c2, 2) < 0) continue;
            if (!mpz_probab_prime_p(c2, 25)) continue;
#pragma omp critical
            {
                R.primes++;
                if (mpz_sizeinbase(c2, 10) > R.max_digits)
                    R.max_digits = mpz_sizeinbase(c2, 10);
            }
        }
        mpz_clear(c2);
    }
    R.seconds = now() - t0;
    if (verbose)
        printf("    %-8s p<=%-7u supp=%-6ld %7zu digits  cand %5ld  surv %5ld"
               "  primes %4ld  %8.2f s\n",
               scheme == 0 ? "A e=1" : scheme == 1 ? "B e=2" : "C rand",
               pmax, R.support, R.digits, R.candidates, R.survivors,
               R.primes, R.seconds);
    free(mark); free(insup); mpz_clear(M); mpz_clear(cand);
    return R;
}

/* find the pmax that brings a strategy to about D digits */
static uint32_t pmax_for_digits(int scheme, uint64_t seed, long D,
                                const uint32_t *odd, size_t nodd)
{
    double acc = 0.0; size_t i;
    rng_seed(seed);
    for (i = 0; i < nodd; i++) {
        int ex = (scheme == 0) ? 1 : (scheme == 1) ? 2 : (int) (rng_next() % 3);
        acc += ex * log10((double) odd[i]);
        if (acc >= (double) D) return odd[i];
    }
    return odd[nodd - 1];
}

int main(int argc, char **argv)
{
    int threads = 4, opt, mode = 0, trials = 3;
    long D = 2000;
    uint32_t pfix = 4000;
    double budget = 600.0;
    uint64_t seed = 12345, Bmax = 20000000ULL;
    size_t nodd, nsmall;
    uint32_t *odd, *small, *spf;
    FILE *csv;

    for (opt = 1; opt < argc; opt++) {
        if      (!strcmp(argv[opt], "-t")     && opt+1 < argc) threads = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-D")     && opt+1 < argc) D = atol(argv[++opt]);
        else if (!strcmp(argv[opt], "-p")     && opt+1 < argc) pfix = strtoul(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-T")     && opt+1 < argc) budget = atof(argv[++opt]);
        else if (!strcmp(argv[opt], "-s")     && opt+1 < argc) seed = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-r")     && opt+1 < argc) trials = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-B")     && opt+1 < argc) Bmax = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-k")     && opt+1 < argc) KMAX = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-mode")  && opt+1 < argc) {
            const char *m = argv[++opt];
            mode = !strcmp(m, "pmax") ? 1 : !strcmp(m, "time") ? 2 : 0;
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    OFF = (1L << KMAX) + 30;
    odd = sieve_primes(4000000, &nodd);
    { size_t k = 0, i; for (i = 0; i < nodd; i++) if (odd[i] > 2) odd[k++] = odd[i]; nodd = k; }
    small = sieve_primes((uint32_t) (Bmax < 4000000000ULL ? Bmax : 4000000000ULL), &nsmall);
    spf = spf_table((uint32_t) OFF + 1);

    csv = fopen("compare.csv", "a");
    fprintf(csv, "mode,scheme,trial,pmax,support,digits,candidates,survivors,"
                 "primes,seconds,max_prime_digits\n");

    if (mode == 0) {
        int sc, tr;
        printf("MODE digits: every strategy grown to about %ld digits\n", D);
        for (sc = 0; sc < 3; sc++)
            for (tr = 0; tr < (sc == 2 ? trials : 1); tr++) {
                uint64_t sd = seed + 1000 * tr;
                uint32_t pm = pmax_for_digits(sc, sd, D, odd, nodd);
                result_t R = run_family(sc, pm, sd, odd, nodd, small, nsmall,
                                        spf, Bmax, 1);
                fprintf(csv, "digits,%c,%d,%u,%ld,%zu,%ld,%ld,%ld,%.3f,%zu\n",
                        "ABC"[sc], tr, R.pmax, R.support, R.digits,
                        R.candidates, R.survivors, R.primes, R.seconds,
                        R.max_digits);
                fflush(csv);
            }
    } else if (mode == 1) {
        int sc, tr;
        printf("MODE pmax: every strategy uses the primes up to %u\n", pfix);
        for (sc = 0; sc < 3; sc++)
            for (tr = 0; tr < (sc == 2 ? trials : 1); tr++) {
                uint64_t sd = seed + 1000 * tr;
                result_t R = run_family(sc, pfix, sd, odd, nodd, small, nsmall,
                                        spf, Bmax, 1);
                fprintf(csv, "pmax,%c,%d,%u,%ld,%zu,%ld,%ld,%ld,%.3f,%zu\n",
                        "ABC"[sc], tr, R.pmax, R.support, R.digits,
                        R.candidates, R.survivors, R.primes, R.seconds,
                        R.max_digits);
                fflush(csv);
            }
    } else {
        int sc;
        printf("MODE time: %.0f s per strategy, climbing in size\n", budget);
        for (sc = 0; sc < 3; sc++) {
            double spent = 0.0; long d = 300; size_t best = 0; long tot = 0;
            while (spent < budget) {
                uint32_t pm = pmax_for_digits(sc, seed, d, odd, nodd);
                result_t R = run_family(sc, pm, seed, odd, nodd, small, nsmall,
                                        spf, Bmax, 1);
                spent += R.seconds; tot += R.primes;
                if (R.max_digits > best) best = R.max_digits;
                fprintf(csv, "time,%c,0,%u,%ld,%zu,%ld,%ld,%ld,%.3f,%zu\n",
                        "ABC"[sc], R.pmax, R.support, R.digits, R.candidates,
                        R.survivors, R.primes, R.seconds, R.max_digits);
                fflush(csv);
                d = (long) (d * 1.4) + 50;
            }
            printf("  -> %c: %ld primes in %.0f s, largest %zu digits\n",
                   "ABC"[sc], tot, spent, best);
        }
    }
    fclose(csv);
    printf("\nresults appended to compare.csv\n");
    free(odd); free(small); free(spf);
    return 0;
}
