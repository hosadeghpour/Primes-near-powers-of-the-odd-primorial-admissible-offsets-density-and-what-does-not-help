# Primes near powers of the odd primorial

Software and data for the paper *Primes near powers of the odd primorial:
admissible offsets, density, and what does not help*.

For an odd prime `p` and a vector of non-negative exponents
`e = (e_q)` over the odd primes `q <= p`, write

```
M(e) = prod q^(e_q)          N = | M(e) +- c |
```

The paper characterises which offsets `c` can accompany `M(e)`, works out the
density of primes in the family, and reports several things that look as though
they should help and do not. Everything below reproduces a specific table or
figure in the paper.

Three exponent strategies appear throughout:

| | |
|---|---|
| **A** | `e_q = 1` for every `q` — the plain primorial |
| **B** | `e_q = 2` for every `q` |
| **C** | `e_q` drawn uniformly from `{0,1,2}` |

## Building

Everything in C needs only GMP and, for the parallel programs, OpenMP:

```
gcc -O3 -march=native -fopenmp -static <file>.c -lgmp -o <file>.exe
```

`indices.c` needs no libraries at all. The Python scripts use only the
standard library, except `pairtest.py`, which uses `gmpy2`.

The programs were run under MSYS2/MinGW-w64 on Windows and under Linux; no
run in the paper took more than ten hours on a four-core laptop.

## The programs

| file | what it does |
|---|---|
| `family_sweep4.c` | The main sweep. For every odd prime up to the *N*th, tests every admissible offset up to `C` on both signs. Chooses the sieve bound per family by timing the host machine. Produces the complete determinations of Section 9. |
| `hunt2.c` | Budget-driven record hunt (Algorithm 2). Given a time budget it measures the machine, solves for the largest size with a realistic chance of a prime inside that budget, sizes the offset range accordingly, and runs until the budget expires. |
| `pow2.c` | The `Pr = 1` branch, sieved on the **exponent** axis (Algorithm 3). The exponent runs to `log2(M)`, several thousand values, not the twenty we first allowed; this program is what that range requires. |
| `mers.c` | The reflected branch `N = 2^n +- M(e) +- l`, where `n` tunes the size bit by bit. Used for the equal-size comparison of the three strategies in Section 7. |
| `compare.c` | Compares strategies A, B and C under three regimes: equal digit count, equal `p`, and equal wall-clock budget. |
| `indices.c` | For every `n` up to a limit: `omega(n)`, `Omega(n)`, the combined index `omega + (1 - omega/Omega)`, and the Mertens boost `prod q/(q-1)`. Sieve-based, two bytes per `n`. |

## The analysis scripts

| file | what it does |
|---|---|
| `analyze.py` | Reads a sweep CSV and produces the sign balance, the distribution over powers of two, the Fortunate-type function, and the observed/predicted ratio. |
| `features.py` | Tests whether any property of the offset predicts its yield. Fourteen features are built in; adding another means writing one function. |
| `pairtest.py` | The matched-pairs experiment of Section 6, which is what showed that the apparent 23% advantage of the `Pr = 1` offsets does not replicate. |
| `make_figures.py`, `mkfig3.py` | Regenerate the five figures of the paper from the result files. |

## Reproducing the paper

```
# Section 9, complete determination (about six hours on four cores)
./family_sweep4.exe -N 700 -e 1 -C 50000 -t 4 -B 20000000 -nosave -q
python analyze.py sweep4.csv --cmax 50000

# Section 5 and 6, the feature scan and the power-of-two question
python features.py sweep4.csv --cmax 50000 --lo 100 --hi 2003
python pairtest.py 100 800 1

# Section 7, three strategies at exactly matched bit lengths
for b in $(seq 1500 25 6000); do
  for s in A B C; do ./mers.exe -bits $b -scheme $s -J 60 -t 4 -B 5000000 -s $b; done
done

# Section 9, a large run
./hunt2.exe -scheme B -digits 20500 -budget 36000 -t 4 -B 1000000000
```

Runs are reproducible: `-digits` fixes the target size rather than letting the
calibration choose it, and `-s` fixes the seed for strategy C.

## Data

`b399051.txt` and `b399052.txt` are the OEIS b-files for the Fortunate-type
functions of Section 10, giving all 699 and 701 computed values respectively,
for `p` from 5 to 5281. Each line is an index and a value, in the format the
OEIS expects.

`run_B_10h.sh` and `run_rented.sh` are the driver scripts for the two large
runs reported in Section 9.

## Sequences

- [A399020](https://oeis.org/A399020) — primes `p` such that `(p#/2)^2 - 2` is prime
- [A399051](https://oeis.org/A399051) — least even `c` with `(prime(n)#/2) - c` prime
- [A399052](https://oeis.org/A399052) — least even `c` with `(prime(n)#/2) + c` prime

## A note on the results

Numbers reported as prime passed 25 rounds of a strong probable-prime test and
are probable primes, not proved primes, except where the paper says otherwise.
The largest is 20501 digits.

The one positive claim an earlier draft of this work made — that offsets with
`Pr = 1` produce about 23% more primes than offsets whose odd part is prime —
did not survive replication. `pairtest.py` is the experiment that killed it,
and the episode is kept in the paper deliberately.

## Acknowledgements

The correspondent posting as firejuggler on mersenneforum pointed out that the
exponent need not be capped, wrote an independent implementation that sieves
the exponent axis, and rediscovered the difference-of-squares obstruction along
the way. The OEIS editors, and Michael S. Branicky in particular, corrected the
three sequence submissions in several places.

## License

MIT for the code. The b-files and result data are released into the public
domain.
