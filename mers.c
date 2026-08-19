/*
 * mers.c  --  the pseudo-Mersenne branch of the family
 *
 *      M(e) = product over odd primes q <= p of q^(e_q)
 *      N    = 2^n  -/+ M  -/+ l ,   n0 <= n <= n0 + J,  l in {0,2,...,20,30}
 *
 * where n0 is the bit length of M, so 2^n is at least as large as M.
 *
 * Why this belongs in the same paper
 * ----------------------------------
 * The protection that defines the family survives unchanged: for every prime
 * q dividing M we have N = 2^n -/+ l (mod q), and that is non-zero whenever
 * q does not divide 2^n -/+ l.  So this is the same family seen from the other
 * side -- instead of subtracting a power of two from the product, we subtract
 * the product from a power of two.  Measured against the Bateman-Horn
 * prediction, all four sign structures sit on the same density:
 *
 *      M - 2^k   1.17        2^n - M   1.04
 *      M + 2^k   1.05        2^n + M   1.00
 *
 * What it adds that the other branch cannot do
 * --------------------------------------------
 * In the branch N = M -/+ 2^k the size of N is fixed by p, and p moves in
 * coarse steps.  Here n is a free parameter, so the size can be set BIT BY
 * BIT.  That is exactly what an equal-size comparison between exponent
 * strategies needs: strategy A, B and C can all be brought to the same digit
 * count instead of merely a similar one.
 *
 * Admissibility differs from the other branch and the difference matters.
 * There, an offset was inadmissible when its odd part met a prime of the
 * support.  Here the relevant quantity is 2^n -/+ l, and whether a prime of
 * the support divides it has to be tested, not assumed.
 *
 * Build:  gcc -O3 -march=native -fopenmp -static mers.c -lgmp -o mers.exe
 * Run:    ./mers.exe -D 2000 -scheme A -J 40 -t 4
 *         ./mers.exe -bits 6644 -scheme C -r 5 -t 4     (exact size targeting)
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

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = (unsigned __int128) r * a % m;
                a = (unsigned __int128) a * a % m; e >>= 1; }
    return r;
}

int main(int argc, char **argv)
{
    int scheme = 0, threads = 4, J = 40, trials = 1, opt, r;
    long D = 2000, bits = 0;
    uint64_t seed = 12345, B = 2000000;
    size_t nodd, nsmall, i;
    uint32_t *odd, *small;
    FILE *csv;

    for (opt = 1; opt < argc; opt++) {
        if      (!strcmp(argv[opt], "-D")      && opt+1 < argc) D = atol(argv[++opt]);
        else if (!strcmp(argv[opt], "-bits")   && opt+1 < argc) bits = atol(argv[++opt]);
        else if (!strcmp(argv[opt], "-J")      && opt+1 < argc) J = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-t")      && opt+1 < argc) threads = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-B")      && opt+1 < argc) B = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-s")      && opt+1 < argc) seed = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-r")      && opt+1 < argc) trials = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-scheme") && opt+1 < argc) {
            char c = argv[++opt][0];
            scheme = (c == 'B' || c == 'b') ? 1 : (c == 'C' || c == 'c') ? 2 : 0;
        }
    }
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    odd = sieve_primes(4000000, &nodd);
    { size_t k = 0; for (i = 0; i < nodd; i++) if (odd[i] > 2) odd[k++] = odd[i]; nodd = k; }
    small = sieve_primes((uint32_t) (B < 4000000000ULL ? B : 4000000000ULL), &nsmall);

    csv = fopen("mers.csv", "a");
    fprintf(csv, "scheme,trial,pmax,support,n0,bits,digits,candidates,survivors,primes,seconds,max_digits\n");

    for (r = 0; r < trials; r++) {
        mpz_t M, N;
        uint32_t *sup; size_t nsup = 0, pmax = 0;
        uint8_t *insup_tab;
        double t0 = now();
        long n0, nn, li, sM, sl, cand_n = 0, surv = 0, found = 0;
        uint8_t *alive;
        size_t maxdig = 0, digits;
        double target = bits ? (double) bits * log10(2.0) : (double) D;

        /* ---- build M with the chosen exponent scheme ---- */
        mpz_init(M); mpz_init(N);
        mpz_set_ui(M, 1);
        sup = malloc(nodd * sizeof(uint32_t));
        rng_seed(seed + 1000 * r);
        { double acc = 0.0;
          for (i = 0; i < nodd; i++) {
              int ex = (scheme == 0) ? 1 : (scheme == 1) ? 2 : (int) (rng_next() % 3);
              if (ex) {
                  mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, odd[i], ex);
                  mpz_mul(M, M, t); mpz_clear(t);
                  sup[nsup++] = odd[i];
                  acc += ex * log10((double) odd[i]);
              }
              pmax = odd[i];
              if (acc >= target) break;
          } }
        insup_tab = calloc((size_t) sup[nsup - 1] + 2, 1);
        for (i = 0; i < nsup; i++) insup_tab[sup[i]] = 1;
        n0 = (long) mpz_sizeinbase(M, 2);
        if (bits && n0 < bits) n0 = bits;          /* exact size targeting  */

        /* ---- candidates: (n, l, sM, sl) ---- */
        alive = calloc((size_t) (J + 1) * N_L * 4, 1);
#define IDX(nn,li,sM,sl) ((((nn)*N_L+(li))*2+(sM))*2+(sl))
        for (nn = 0; nn <= J; nn++)
            for (li = 0; li < N_L; li++)
                for (sM = 0; sM < 2; sM++)
                    for (sl = 0; sl < 2; sl++) {
                        if (L_LIST[li] == 0 && sl == 1) continue;   /* dedup */
                        alive[IDX(nn,li,sM,sl)] = 1; cand_n++;
                    }

        /* ---- kill candidates that a prime of the support divides ------- */
        for (i = 0; i < nsup; i++) {
            uint64_t q = sup[i], pw = powmod(2, (uint64_t) n0, q);
            for (nn = 0; nn <= J; nn++) {
                for (li = 0; li < N_L; li++)
                    for (sl = 0; sl < 2; sl++) {
                        long l = (sl ? -1 : 1) * L_LIST[li];
                        long v = (long) (pw % q) + l;
                        v %= (long) q; if (v < 0) v += q;
                        if (v == 0) { alive[IDX(nn,li,0,sl)] = 0;
                                      alive[IDX(nn,li,1,sl)] = 0; }
                    }
                pw = pw * 2 % q;
            }
        }

        /* ---- sieve with the primes above the support -------------------- */
#pragma omp parallel for schedule(static)
        for (long ii = 0; ii < (long) nsmall; ii++) {
            uint64_t q = small[ii], pw, rM;
            long nn2, li2, sM2, sl2;
            if (q > B) continue;
            if (q == 2 || (q <= sup[nsup - 1] && insup_tab[q])) continue;
            rM = mpz_fdiv_ui(M, (unsigned long) q);
            pw = powmod(2, (uint64_t) n0, q);
            for (nn2 = 0; nn2 <= J; nn2++) {
                for (li2 = 0; li2 < N_L; li2++)
                    for (sM2 = 0; sM2 < 2; sM2++)
                        for (sl2 = 0; sl2 < 2; sl2++) {
                            long l = (sl2 ? -1 : 1) * L_LIST[li2];
                            long v;
                            if (!alive[IDX(nn2,li2,sM2,sl2)]) continue;
                            v = (long) (pw % q) + (sM2 ? (long) rM : -(long) rM) + l;
                            v %= (long) q; if (v < 0) v += q;
                            if (v == 0) alive[IDX(nn2,li2,sM2,sl2)] = 0;
                        }
                pw = pw * 2 % q;
            }
        }

        /* ---- test the survivors ---------------------------------------- */
#pragma omp parallel
        {
            mpz_t T, PW;
            mpz_init(T); mpz_init(PW);
#pragma omp for schedule(dynamic, 1)
            for (nn = 0; nn <= J; nn++) {
                long li2, sM2, sl2;
                mpz_ui_pow_ui(PW, 2, (unsigned long) (n0 + nn));
                for (li2 = 0; li2 < N_L; li2++)
                    for (sM2 = 0; sM2 < 2; sM2++)
                        for (sl2 = 0; sl2 < 2; sl2++) {
                            long l;
                            if (!alive[IDX(nn,li2,sM2,sl2)]) continue;
                            l = (sl2 ? -1 : 1) * L_LIST[li2];
                            mpz_set(T, PW);
                            if (sM2) mpz_add(T, T, M); else mpz_sub(T, T, M);
                            if (l >= 0) mpz_add_ui(T, T, (unsigned long) l);
                            else        mpz_sub_ui(T, T, (unsigned long) (-l));
                            mpz_abs(T, T);
#pragma omp atomic
                            surv++;
                            if (mpz_cmp_ui(T, 2) < 0) continue;
                            if (!mpz_probab_prime_p(T, 25)) continue;
#pragma omp critical
                            {
                                found++;
                                if (mpz_sizeinbase(T, 10) > maxdig)
                                    maxdig = mpz_sizeinbase(T, 10);
                                gmp_printf("   PRIME 2^%ld %c M %+ld  (%zu digits)\n",
                                           n0 + nn, sM2 ? '+' : '-', l,
                                           mpz_sizeinbase(T, 10));
                            }
                        }
            }
            mpz_clear(T); mpz_clear(PW);
        }
        mpz_ui_pow_ui(N, 2, (unsigned long) (n0 + J));
        digits = mpz_sizeinbase(N, 10);
        printf("%c trial %d: p<=%zu supp=%zu  n0=%ld  ~%zu digits  "
               "cand %ld  surv %ld  primes %ld  %.2f s\n",
               "ABC"[scheme], r, (size_t) pmax, nsup, n0, digits,
               cand_n, surv, found, now() - t0);
        fprintf(csv, "%c,%d,%zu,%zu,%ld,%ld,%zu,%ld,%ld,%ld,%.3f,%zu\n",
                "ABC"[scheme], r, (size_t) pmax, nsup, n0, n0 + J, digits,
                cand_n, surv, found, now() - t0, maxdig);
        fflush(csv);
        free(alive); free(sup); free(insup_tab); mpz_clear(M); mpz_clear(N);
    }
    fclose(csv);
    printf("results appended to mers.csv\n");
    free(odd); free(small);
    return 0;
}
