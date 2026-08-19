/*
 * pow2.c  --  the branch  N = | M(e) - 2^n |,  sieved on the EXPONENT axis
 *
 * The other programs in this collection sieve the offset axis: for a prime q
 * they compute M mod q once and strike one arithmetic progression of offsets.
 * That is right when the offsets are dense.  It is wrong here.  The offsets
 * 2^n are exponentially sparse, and the useful range of n is not a few dozen
 * but all of 1 .. log2(M) -- several thousand values for a number of a few
 * hundred digits.  Capping n at 20, as we did at first, threw away almost the
 * whole family.
 *
 * The right sieve works on n directly.  For a prime q not dividing M,
 *
 *      q | M - 2^n     <=>     2^n = M  (mod q),
 *
 * and the solutions form an arithmetic progression of difference ord_q(2).
 * Rather than compute a discrete logarithm, we simply walk n upward keeping
 * 2^n mod q in a register and compare it with M mod q: one multiplication and
 * one reduction per (q, n) pair, no table, no factoring of q-1, and every
 * solution below n_max is found on the way. The idea of sieving this axis is
 * due to a correspondent on mersenneforum who implemented it in Python.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -static pow2.c -lgmp -o pow2.exe
 * Run:
 *   ./pow2.exe -p 5000 -scheme A -t 4 -B 5000000
 *   ./pow2.exe -digits 3000 -scheme B -t 4 -B 20000000
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

static uint64_t rng_state;
static void rng_seed(uint64_t s) { rng_state = s * 6364136223846793005ULL + 1; }
static uint32_t rng_next(void)
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17; return (uint32_t) (rng_state >> 32);
}

int main(int argc, char **argv)
{
    uint32_t pmax = 5000, *odd, *small;
    uint64_t B = 5000000;
    int scheme = 0, threads = 4, opt;
    long target_digits = 0, nmax, alive_cnt = 0, found = 0, support = 0;
    uint32_t plimit;
    size_t nodd, nsmall, i;
    mpz_t M, cand;
    uint8_t *alive, *insup;
    double t0;
    FILE *out;

    for (opt = 1; opt < argc; opt++) {
        if      (!strcmp(argv[opt], "-p")      && opt+1 < argc) pmax = strtoul(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-digits") && opt+1 < argc) target_digits = atol(argv[++opt]);
        else if (!strcmp(argv[opt], "-t")      && opt+1 < argc) threads = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-B")      && opt+1 < argc) B = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-scheme") && opt+1 < argc) {
            char c = argv[++opt][0];
            scheme = (c=='B'||c=='b') ? 1 : (c=='C'||c=='c') ? 2 : 0;
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    odd = sieve_primes(8000000, &nodd);
    { size_t k = 0; for (i = 0; i < nodd; i++) if (odd[i] > 2) odd[k++] = odd[i]; nodd = k; }

    /* ---- build M ---- */
    mpz_init(M); mpz_init(cand);
    mpz_set_ui(M, 1);
    insup = calloc((size_t) odd[nodd-1] + 2, 1);
    rng_seed(12345);
    plimit = pmax;
    { double acc = 0.0;
      for (i = 0; i < nodd; i++) {
          int ex = (scheme == 0) ? 1 : (scheme == 1) ? 2 : (int) (rng_next() % 3);
          if (target_digits == 0 && odd[i] > plimit) break;
          if (ex) {
              mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, odd[i], ex);
              mpz_mul(M, M, t); mpz_clear(t);
              insup[odd[i]] = 1; support++;
              acc += ex * log10((double) odd[i]);
          }
          pmax = odd[i];
          if (target_digits && acc >= (double) target_digits) break;
      } }

    /* ---- the exponent range is set by M itself, not by a cutoff we choose --- */
    nmax = (long) mpz_sizeinbase(M, 2) - 2;
    printf("scheme %c: p <= %u, support %ld primes, M has %zu digits\n",
           "ABC"[scheme], pmax, support, mpz_sizeinbase(M, 10));
    printf("exponent range n = 1 .. %ld   (that is log2(M), not an arbitrary cap)\n",
           nmax);

    alive = malloc(nmax + 2);
    memset(alive, 1, nmax + 2);
    alive[0] = 0;
    /* if M is a perfect square, even n factors as a difference of squares */
    if (scheme == 1) { long n; for (n = 2; n <= nmax; n += 2) alive[n] = 0; }

    /* ---- sieve the exponent axis ---- */
    t0 = now();
    small = sieve_primes((uint32_t) (B < 4000000000ULL ? B : 4000000000ULL), &nsmall);
#pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long) nsmall; ii++) {
        uint64_t q = small[ii], r, t;
        long n;
        if (q == 2) continue;
        if (q <= pmax && insup[q]) continue;      /* cannot divide M - 2^n */
        r = mpz_fdiv_ui(M, (unsigned long) q);
        if (r == 0) continue;
        t = 1 % q;
        for (n = 0; n <= nmax; n++) {
            if (t == r) alive[n] = 0;             /* q | M - 2^n            */
            t += t; if (t >= q) t -= q;           /* one add, no multiply   */
        }
    }
    { long n; for (n = 1; n <= nmax; n++) alive_cnt += alive[n]; }
    printf("sieve to %llu: %.1f s -> %ld survivors of %ld  (expected primes %.2f)\n",
           (unsigned long long) B, now() - t0, alive_cnt, nmax,
           alive_cnt * EG * log((double) B) / (mpz_sizeinbase(M,10) * log(10.0)));

    /* ---- test ---- */
    out = fopen("pow2_primes.txt", "a");
    t0 = now();
#pragma omp parallel
    {
        mpz_t c2, pw;
        mpz_init(c2); mpz_init(pw);
#pragma omp for schedule(dynamic, 1)
        for (long n = 1; n <= nmax; n++) {
            if (!alive[n]) continue;
            mpz_ui_pow_ui(pw, 2, (unsigned long) n);
            mpz_sub(c2, M, pw);
            mpz_abs(c2, c2);
            if (mpz_cmp_ui(c2, 2) < 0) continue;
            if (!mpz_probab_prime_p(c2, 25)) continue;
#pragma omp critical
            {
                found++;
                printf("  PRIME  M - 2^%ld   (%zu digits)\n", n, mpz_sizeinbase(c2,10));
                gmp_fprintf(out, "# scheme=%c p<=%u n=%ld %zu digits\n%Zd\n",
                            "ABC"[scheme], pmax, n, mpz_sizeinbase(c2,10), c2);
                fflush(stdout); fflush(out);
            }
        }
        mpz_clear(c2); mpz_clear(pw);
    }
    fclose(out);
    printf("\n%ld primes from %ld tests in %.1f s\n", found, alive_cnt, now() - t0);
    { FILE *f = fopen("pow2.csv", "a");
      fprintf(f, "%c,%u,%ld,%zu,%ld,%ld,%ld\n", "ABC"[scheme], pmax, support,
              mpz_sizeinbase(M,10), nmax, alive_cnt, found);
      fclose(f); }
    mpz_clear(M); mpz_clear(cand);
    free(alive); free(insup); free(odd); free(small);
    return 0;
}
