/*
 * family_sweep4.c  --  the sweep that produces the paper's dataset
 *
 *        N = | P^e - c | ,   P = 3*5*7*...*p ,   c even and admissible
 *
 * for every odd prime p in a range, both signs, all admissible c <= C.
 * The subfamily c = 2^k is included automatically (it is the case where the
 * odd part of c is 1), so one run gives both the distinctive slice and the
 * general family, measured under identical conditions -- which is what the
 * comparison between them requires.
 *
 * What is new compared with family_sweep2:
 *
 *  1. PER-FAMILY AUTO SIEVE BOUND.  A fixed bound is wrong: the right depth
 *     depends on how expensive one test is, and that grows with p.  Each
 *     family times one reduction and uses the measured cost of a test to set
 *     B = T_test * (live candidates) / T_reduction, which is where sieving
 *     one more prime stops paying for itself.  Small families get a shallow
 *     sieve, large ones a deep one.
 *
 *  2. ADMISSIBILITY PRECOMPUTED ONCE.  The smallest prime factor of the odd
 *     part of every c is computed a single time; each family then only has to
 *     compare it with p.
 *
 *  3. TIMING PER FAMILY, written to sweep4_time.csv along with the number of
 *     admissible offsets, survivors, and primes.  This is what lets you say
 *     in the paper how long the computation took and how it scales.
 *
 *  4. RESUME.  -f sets the first family index, so an interrupted run can be
 *     continued, or a run split across machines.  Output is appended.
 *
 *  5. DEGENERATE FAMILIES FLAGGED.  When P^e is not much larger than C the
 *     numbers tested are no longer "near P^e" and should be excluded from the
 *     analysis; those rows get degenerate=1.
 *
 * Build:
 *   gcc -O3 -march=native -fopenmp -static family_sweep4.c -lgmp -o family_sweep4.exe
 * Run:
 *   ./family_sweep4.exe -N 1000 -e 2 -C 50000 -t 4
 *   ./family_sweep4.exe -N 1000 -e 1 -C 50000 -t 4
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

static void product_tree(mpz_t out, const uint32_t *v, size_t n)
{
    if (n == 0) { mpz_set_ui(out, 1); return; }
    if (n == 1) { mpz_set_ui(out, v[0]); return; }
    { mpz_t l, r; mpz_init(l); mpz_init(r);
      product_tree(l, v, n / 2); product_tree(r, v + n / 2, n - n / 2);
      mpz_mul(out, l, r); mpz_clear(l); mpz_clear(r); }
}

int main(int argc, char **argv)
{
    int N = 1000, e = 2, threads = 4, first = 1, opt, quiet = 0, savetxt = 1;
    uint32_t C = 50000;
    uint64_t Bmax = 200000000ULL;
    size_t nprimes, nsmall, i;
    uint32_t *primes, *small, *oddspf;
    uint8_t *ppow_minus, *ppow_plus;
    long n_c;
    FILE *csv, *tcsv, *fp;

    for (opt = 1; opt < argc; opt++) {
        if      (!strcmp(argv[opt], "-N") && opt+1 < argc) N       = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-e") && opt+1 < argc) e       = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-C") && opt+1 < argc) C       = strtoul(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-t") && opt+1 < argc) threads = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-f") && opt+1 < argc) first   = atoi(argv[++opt]);
        else if (!strcmp(argv[opt], "-B") && opt+1 < argc) Bmax    = strtoull(argv[++opt],0,10);
        else if (!strcmp(argv[opt], "-q")) quiet = 1;
        /* a full sweep can emit millions of primes; writing them all out costs
           gigabytes and slows the run.  The CSV keeps p, sign and c, and every
           number is reconstructible from those three. */
        else if (!strcmp(argv[opt], "-nosave")) savetxt = 0;
    }
    if (C & 1) C--;
    n_c = C / 2;
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif

    /* ---- primes ------------------------------------------------------- */
    {   double n = N + 2;
        uint32_t bound = (uint32_t)(n * (log(n) + log(log(n))) + 20);
        primes = sieve_primes(bound, &nprimes);
        while ((int) nprimes < N + 1) {
            free(primes); bound *= 2; primes = sieve_primes(bound, &nprimes);
        }
    }
    small = sieve_primes((uint32_t) (Bmax < 4294967295ULL ? Bmax : 4294967295ULL),
                         &nsmall);

    /* ---- smallest prime factor of the odd part of every even c -------- */
    {
        uint32_t *spf = malloc(((size_t) C + 2) * sizeof(uint32_t));
        uint32_t a, b;
        long j;
        for (a = 0; a <= C + 1; a++) spf[a] = a;
        for (a = 2; (uint64_t) a * a <= C + 1; a++)
            if (spf[a] == a)
                for (b = a * a; b <= C + 1; b += a) if (spf[b] == b) spf[b] = a;
        oddspf = malloc((n_c + 1) * sizeof(uint32_t));
        for (j = 1; j <= n_c; j++) {
            uint32_t m = (uint32_t) (2 * j);
            while ((m & 1) == 0) m >>= 1;
            oddspf[j] = (m == 1) ? 0xFFFFFFFFu : spf[m];   /* 0xFFFF.. = 2^k */
        }
        free(spf);
    }

    /* ---- perfect-power filter -----------------------------------------
       If d > 1 divides e and c is a perfect d-th power, then
           P^e - c = (P^(e/d))^d - (c^(1/d))^d
       is divisible by P^(e/d) - c^(1/d), and for odd d the same identity
       with a plus sign divides P^e + c.  These offsets are composite by
       algebra and must not be tested at all.  (c = 2^k with even k and
       e = 2 is only the smallest case of this; every even perfect square
       is excluded too.)                                                 */
    {
        long j;
        int d;
        ppow_minus = calloc(n_c + 1, 1);
        ppow_plus  = calloc(n_c + 1, 1);
        for (d = 2; d <= e; d++) {
            if (e % d) continue;
            for (j = 2; ; j += 1) {
                double v = pow((double) j, (double) d);
                unsigned long c;
                if (v > (double) C + 1) break;
                c = (unsigned long) (v + 0.5);
                if ((c & 1) || c < 2 || c > (unsigned long) C) continue;
                ppow_minus[c / 2] = 1;
                if (d & 1) ppow_plus[c / 2] = 1;
            }
        }
    }

    csv  = fopen("sweep4.csv",  "a");
    tcsv = fopen("sweep4_time.csv", "a");
    fp   = savetxt ? fopen("sweep4_primes.txt", "a") : NULL;
    if (first == 1) {
        fprintf(csv,  "p,digits,sign,c,pow2\n");
        fprintf(tcsv, "p,digits,degenerate,admissible,sieve_bound,survivors,primes,seconds\n");
    }
    printf("N = |P^%d -/+ c|,  families %d..%d (p <= %u),  c <= %u,  %d threads\n",
           e, first, N, primes[N], C, threads);
    if (!quiet)
        printf("%-8s %8s %5s %7s %11s %7s %6s %9s\n", "p", "digits", "deg",
               "adm", "sieve B", "surv", "prime", "seconds");

#pragma omp parallel for schedule(dynamic, 1) ordered
    for (i = (size_t) first; i <= (size_t) N; i++) {
        mpz_t P, Q, cand;
        uint32_t p = primes[i];
        uint8_t *alive[2];
        long j, adm = 0, surv = 0, nprime = 0;
        double t0 = now(), t_test, t_red, secs;
        uint64_t B;
        size_t digits, k;
        int degenerate;

        mpz_init(P); mpz_init(Q); mpz_init(cand);
        product_tree(P, primes + 1, i);
        mpz_pow_ui(Q, P, e);
        digits = mpz_sizeinbase(Q, 10);
        degenerate = (mpz_sizeinbase(Q, 2) < 2 * (size_t) (log((double) C) / log(2.0) + 1));

        alive[0] = malloc(n_c + 1);
        alive[1] = malloc(n_c + 1);
        for (j = 1; j <= n_c; j++) {
            uint8_t ok = (oddspf[j] > p);           /* admissible for this p */
            /* the algebraic filter is skipped for degenerate families,
               where the cofactor can legitimately be 1                  */
            alive[0][j] = ok && !(ppow_minus[j] && !degenerate);
            alive[1][j] = ok && !(ppow_plus[j]  && !degenerate);
            adm += alive[0][j] + alive[1][j];
        }

        /* ---- calibrate this family, then pick the sieve depth ---------- */
        { double s = now();
          mpz_sub_ui(cand, Q, 2);
          mpz_probab_prime_p(cand, 1);
          t_test = now() - s;
          s = now();
          { uint32_t q; volatile unsigned long acc = 0;
            for (q = 1000003; q < 1000003 + 500; q += 2) acc += mpz_fdiv_ui(Q, q); }
          t_red = (now() - s) / 500.0;
          if (t_red <= 0) t_red = 1e-7;
          { double b = t_test * (double) adm / t_red;
            if (b < 100000.0) b = 100000.0;
            if (b > (double) Bmax) b = (double) Bmax;
            B = (uint64_t) b; } }

        /* ---- sieve the offset axis ------------------------------------ */
        {   uint64_t q_ceiling = 0xFFFFFFFFFFFFFFFFull;
            if (mpz_fits_ulong_p(Q)) {
                unsigned long qv = mpz_get_ui(Q);
                q_ceiling = (qv > (unsigned long) C + 2) ? qv - C : 0;
            }
            for (k = 0; k < nsmall; k++) {
                uint32_t q = small[k];
                uint64_t t, c0;
                int s;
                if (q > B || q >= q_ceiling) break;
                if (q <= p) continue;
                { uint64_t r = mpz_fdiv_ui(Q, q);
                  for (s = 0; s < 2; s++) {
                      t = s ? (q - r) % q : r;
                      if (t == 0) t = q;
                      c0 = (t & 1) ? t + q : t;
                      for (; c0 <= C; c0 += 2 * q) alive[s][c0 / 2] = 0;
                  } }
            }
        }

        /* ---- test the survivors --------------------------------------- */
        for (j = 1; j <= n_c; j++) {
            int s;
            for (s = 0; s < 2; s++) {
                if (!alive[s][j]) continue;
                surv++;
                if (s) mpz_add_ui(cand, Q, (unsigned long) (2 * j));
                else   mpz_sub_ui(cand, Q, (unsigned long) (2 * j));
                mpz_abs(cand, cand);
                if (mpz_cmp_ui(cand, 2) < 0) continue;
                if (!mpz_probab_prime_p(cand, 25)) continue;
                nprime++;
#pragma omp critical
                {
                    fprintf(csv, "%u,%zu,%c,%ld,%d\n", p, digits, s ? '+' : '-',
                            2 * j, oddspf[j] == 0xFFFFFFFFu ? 1 : 0);
                    if (fp)
                        gmp_fprintf(fp, "# p=%u e=%d c=%c%ld %zu digits\n%Zd\n",
                                    p, e, s ? '+' : '-', 2 * j,
                                    mpz_sizeinbase(cand, 10), cand);
                }
            }
        }
        secs = now() - t0;

#pragma omp ordered
        {
            fprintf(tcsv, "%u,%zu,%d,%ld,%llu,%ld,%ld,%.3f\n", p, digits,
                    degenerate, adm, (unsigned long long) B, surv, nprime, secs);
            if (!quiet)
                printf("%-8u %8zu %5d %7ld %11llu %7ld %6ld %9.3f\n",
                       p, digits, degenerate, adm, (unsigned long long) B,
                       surv, nprime, secs);
            fflush(stdout); fflush(csv); fflush(tcsv);
            if (fp) fflush(fp);
        }
        free(alive[0]); free(alive[1]);
        mpz_clear(P); mpz_clear(Q); mpz_clear(cand);
    }
    fclose(csv); fclose(tcsv); if (fp) fclose(fp);
    free(primes); free(small); free(oddspf);
    free(ppow_minus); free(ppow_plus);
    return 0;
}
