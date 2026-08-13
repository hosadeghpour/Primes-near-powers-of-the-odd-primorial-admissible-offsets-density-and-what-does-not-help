/*
 * hunt2.c  --  largest prime of the family within a fixed time budget,
 *              for each of the three exponent strategies
 *
 *      M(e) = product of q^(e_q) over the odd primes q <= p
 *      N    = | M(e) - c | ,  c even and admissible
 *
 *      A : e_q = 1        B : e_q = 2        C : e_q uniform on {0,1,2}
 *
 * Given -budget T seconds, the program chooses the target size itself.  That
 * matters: a hunt that climbs blindly wastes the budget.  Testing a candidate
 * of D digits costs about a*D^2.24, and one prime costs
 *
 *      tests = ln N / (e^gamma * ln B)          (B = sieve bound)
 *
 * so the largest D with a realistic chance inside T satisfies
 * tests(D) * a*D^2.24 = T.  The constant a is measured on the machine at
 * startup rather than assumed, so the answer is right for whatever hardware
 * this runs on.
 *
 * The offset bound C is then set so that the number of admissible offsets is
 * a few times the number of tests expected to be needed -- enough shots to
 * make the expected count about two, without wasting memory.
 *
 * Build:  gcc -O3 -march=native -fopenmp -static hunt2.c -lgmp -o hunt2.exe
 * Run:    ./hunt2.exe -scheme A -budget 3600 -t 4
 *         ./hunt2.exe -scheme C -budget 3600 -t 4 -s 7
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

#define EG 1.7810724179901979

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

static uint64_t rng_state;
static void rng_seed(uint64_t s) { rng_state = s * 6364136223846793005ULL + 1; }
static uint32_t rng_next(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17; return (uint32_t) (rng_state >> 32);
}

static int exponent_of(int scheme) /* call in prime order */
{
    return scheme == 0 ? 1 : scheme == 1 ? 2 : (int) (rng_next() % 3);
}

/* build M for a scheme, stopping when it reaches `target` digits */
static uint32_t build_M(mpz_t M, uint8_t *insup, int scheme, uint64_t seed,
                        double target, const uint32_t *odd, size_t nodd,
                        long *support)
{
    double acc = 0.0; size_t i; uint32_t last = 3;
    mpz_set_ui(M, 1); *support = 0;
    rng_seed(seed);
    for (i = 0; i < nodd; i++) {
        int ex = exponent_of(scheme);
        last = odd[i];
        if (ex) {
            mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, odd[i], ex);
            mpz_mul(M, M, t); mpz_clear(t);
            if (insup) insup[odd[i]] = 1;
            (*support)++;
            acc += ex * log10((double) odd[i]);
        }
        if (acc >= target) break;
    }
    if (acc < target)
        fprintf(stderr, "WARNING: ran out of primes at p=%u, M is only "
                        "%.0f of the %.0f target digits\n", last, acc, target);
    return last;
}

int main(int argc, char **argv)
{
    int scheme = 0, threads = 4, opt;
    double budget = 3600.0, yield_target = 2.5, fixed_digits = 0.0;
    uint64_t seed = 12345, B = 100000000ULL;
    size_t nodd, nsmall, i;
    uint32_t *odd, *small, *spf = NULL, pmax;
    mpz_t M, cand;
    uint8_t *insup, *alive[2];
    long n_c, C, adm = 0, surv = 0, found = 0, support = 0, tested = 0;
    long from = 2, last_j = 0;
    double a_const, D, t0, tstart;
    size_t digits;
    FILE *out;

    for (opt = 1; opt < argc; opt++) {
        if      (!strcmp(argv[opt], "-t")      && opt+1 < argc) threads = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-budget") && opt+1 < argc) budget = atof(argv[++opt]);
        else if (!strcmp(argv[opt], "-s")      && opt+1 < argc) seed = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-yield")  && opt+1 < argc) yield_target = atof(argv[++opt]);
        else if (!strcmp(argv[opt], "-from")   && opt+1 < argc) from = atol(argv[++opt]);
        else if (!strcmp(argv[opt], "-digits") && opt+1 < argc) fixed_digits = atof(argv[++opt]);
        else if (!strcmp(argv[opt], "-B")      && opt+1 < argc) B = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-scheme") && opt+1 < argc) {
            char c = argv[++opt][0];
            scheme = (c=='B'||c=='b') ? 1 : (c=='C'||c=='c') ? 2 : 0;
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    tstart = now();
    odd = sieve_primes(8000000, &nodd);
    { size_t k = 0; for (i = 0; i < nodd; i++) if (odd[i] > 2) odd[k++] = odd[i]; nodd = k; }

    /* ---------- calibrate the machine at a 2000-digit reference ---------- */
    mpz_init(M); mpz_init(cand);
    build_M(M, NULL, scheme, seed, 2000.0, odd, nodd, &support);
{   /* Time one modular exponentiation, which is exactly what a survivor
           costs.  Timing mpz_probab_prime_p on arbitrary neighbours of M is
           wrong: most of them have a small factor and GMP throws them out by
           trial division in microseconds, so the calibration reads zero. */
        double s, el; mpz_t e, r, b;
        mpz_init(e); mpz_init(r); mpz_init_set_ui(b, 3);
        mpz_sub_ui(e, M, 1);
        s = now();
        mpz_powm(r, b, e, M);
        el = now() - s;
        a_const = el / pow(2000.0, 2.24);
        mpz_clear(e); mpz_clear(r); mpz_clear(b); }
    printf("calibration: one 2000-digit test = %.4f s  ->  a = %.3e\n",
           a_const * pow(2000.0, 2.24), a_const);

    /* ---------- choose the target size from the budget ------------------ */
    {   double lo = 500, hi = 500000, mid;
        int it;
        for (it = 0; it < 80; it++) {
            double p_est, tests, cost;
            mid = 0.5 * (lo + hi);
            p_est = mid * log(10.0) / (scheme == 1 ? 2.0 : 1.0);
            if (p_est < 10) p_est = 10;
            tests = mid * log(10.0) / (EG * log((double) B));
            cost  = tests * a_const * pow(mid, 2.24);
            if (cost * yield_target < budget) lo = mid; else hi = mid;
        }
        D = lo;
        if (D > 200000) D = 200000;          /* sanity cap */
        /* an explicit size makes the run reproducible: the calibration varies
           a little from run to run, which otherwise shifts the target and with
           it the whole candidate set */
        if (fixed_digits > 0) D = fixed_digits;
    }
    printf("budget %.0f s, aiming for ~%.1f primes  ->  target size %.0f digits\n",
           budget, yield_target, D);

    /* ---------- build the real M ---------------------------------------- */
    insup = calloc((size_t) odd[nodd - 1] + 2, 1);   /* exactly as big as needed */
    pmax = build_M(M, insup, scheme, seed, D, odd, nodd, &support);
    digits = mpz_sizeinbase(M, 10);
    printf("scheme %c: p <= %u, support %ld primes, N has %zu digits\n",
           "ABC"[scheme], pmax, support, digits);

    /* ---------- choose the offset bound so there are enough shots -------- */
    {   double tests = digits * log(10.0) / (EG * log((double) B));
        long need = (long) (3 * tests) + 100;
        C = 200000;
        while (C < 400000000L) {
            long cnt = 0, j;
            free(spf); spf = spf_table((uint32_t) C + 1);
            for (j = 1; j <= C / 2; j++) {
                uint32_t m = (uint32_t) (2 * j);
                while ((m & 1) == 0) m >>= 1;
                if (m == 1 || spf[m] > pmax) cnt++;
            }
            if (2 * cnt >= need) { adm = 2 * cnt; break; }
            C *= 4;
        }
        printf("need ~%ld tests, offset bound C = %ld gives %ld admissible\n",
               need, C, adm);
    }
    n_c = C / 2;
    alive[0] = malloc(n_c + 1); alive[1] = malloc(n_c + 1);
    { long j; int msq = (scheme == 1);      /* all exponents 2 -> M is square */
      for (j = 1; j <= n_c; j++) {
          uint32_t m = (uint32_t) (2 * j);
          int ok;
          while ((m & 1) == 0) m >>= 1;
          ok = (m == 1 || spf[m] > pmax);
          alive[0][j] = alive[1][j] = ok;
          if (ok && msq) {
              unsigned long c = 2UL * j, r = (unsigned long) (sqrt((double) c) + 0.5);
              if (r * r == c) alive[0][j] = 0;    /* difference of two squares */
          }
      } }

    /* ---------- sieve ---------------------------------------------------- */
    t0 = now();
    small = sieve_primes((uint32_t) (B < 4000000000ULL ? B : 4000000000ULL), &nsmall);
#pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long) nsmall; ii++) {
        uint32_t q = small[ii]; uint64_t r, t, c0; int s;
        if (q == 2) continue;
        if (q <= pmax && insup[q]) continue;
        r = mpz_fdiv_ui(M, q);
        for (s = 0; s < 2; s++) {
            t = s ? (q - r) % q : r;
            if (t == 0) t = q;
            c0 = (t & 1) ? t + q : t;
            for (; c0 <= (uint64_t) C; c0 += 2ULL * q) alive[s][c0 / 2] = 0;
        }
    }
    { long j; for (j = 1; j <= n_c; j++) surv += alive[0][j] + alive[1][j]; }
    printf("sieve: %.1f s -> %ld survivors,  expected primes %.2f\n",
           now() - t0, surv,
           surv * EG * log((double) B) / (digits * log(10.0)));

    /* ---------- test until the budget runs out --------------------------- */
    out = fopen("hunt2_primes.txt", "a");
    t0 = now();
#pragma omp parallel
    {
        mpz_t c2; mpz_init(c2);
#pragma omp for schedule(dynamic, 1)
        for (long j = from / 2; j <= n_c; j++) {
            int s;
            if (now() - tstart > budget) continue;
            if (j > last_j) {
#pragma omp critical (progress)
                if (j > last_j) last_j = j;
            }
            for (s = 0; s < 2; s++) {
                if (!alive[s][j]) continue;
                if (s) mpz_add_ui(c2, M, (unsigned long) (2 * j));
                else   mpz_sub_ui(c2, M, (unsigned long) (2 * j));
                mpz_abs(c2, c2);
                if (mpz_cmp_ui(c2, 2) < 0) continue;
#pragma omp atomic
                tested++;
                if (!mpz_probab_prime_p(c2, 25)) continue;
#pragma omp critical
                {
                    found++;
                    {   /* c = 2^k * Pr is the paper's notation; by the
                           admissibility theorem Pr is 1 or a prime > p */
                        unsigned long cc = 2UL * j, Pr = 2UL * j; int kk = 0;
                        while ((Pr & 1) == 0) { Pr >>= 1; kk++; }
                        printf("  PRIME  scheme %c  c = %c%lu = 2^%d * %lu   "
                               "%zu digits\n", "ABC"[scheme], s ? '+' : '-',
                               cc, kk, Pr, mpz_sizeinbase(c2, 10));
                        gmp_fprintf(out, "# scheme=%c p<=%u sign=%c k=%d Pr=%lu "
                                    "c=%lu %zu digits\n%Zd\n", "ABC"[scheme],
                                    pmax, s ? '+' : '-', kk, Pr, cc,
                                    mpz_sizeinbase(c2, 10), c2); }
                    fflush(stdout); fflush(out);
                }
            }
        }
        mpz_clear(c2);
    }
    fclose(out);
    printf("\nscheme %c: %ld primes of %zu digits in %.0f s\n",
           "ABC"[scheme], found, digits, now() - tstart);
    printf("   tested %ld of %ld survivors (%.0f%%)  ->  expected %.2f primes"
           " for the tests actually done\n", tested, surv,
           100.0 * tested / (surv ? surv : 1),
           tested * EG * log((double) B) / (digits * log(10.0)));
    {   double need = digits * log(10.0) / (EG * log((double) B));
        if ((double) tested < need)
            printf("   NOTE: only %ld tests were done but about %.0f are needed\n"
                   "         per prime at this size - raise the budget, lower the\n"
                   "         target, or use more threads.\n", tested, need);
        else if (tested < surv)
            printf("   (the candidate list was deliberately larger than needed;\n"
                   "    %ld untested candidates remain, use -from %ld to continue)\n",
                   surv - tested, 2 * last_j + 2); }
    { FILE *f = fopen("hunt2.csv", "a");
      fprintf(f, "%c,%u,%ld,%zu,%ld,%ld,%ld,%ld,%.1f\n", "ABC"[scheme], pmax,
              support, digits, adm, surv, tested, found, now() - tstart);
      fclose(f); }
    mpz_clear(M); mpz_clear(cand);
    free(alive[0]); free(alive[1]); free(insup); free(spf); free(odd); free(small);
    return 0;
}
