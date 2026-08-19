"""Matched-pairs test: does an offset that is a pure power of two do better?

For each family M = product of odd primes <= p, and for each bit length b,
arm A takes c = 2^(b-1) and arm B takes c = 2*Pr with Pr a prime of b-1 bits.
The two arms are therefore matched in size, in family, and in count -- the
only difference is the odd part of the offset.
"""
import sys, random, time
from gmpy2 import mpz, is_prime, next_prime

def primes_upto(n):
    f = bytearray([1])*(n+1); f[0]=f[1]=0
    i=2
    while i*i<=n:
        if f[i]: f[i*i::i]=bytearray(len(f[i*i::i]))
        i+=1
    return [i for i in range(n+1) if f[i]]

def run(plo, phi, seed=1, verbose=True):
    random.seed(seed)
    pr = primes_upto(phi)
    M = mpz(1)
    for q in pr:
        if q > 2 and q <= plo: M *= q
    a_cand=a_hit=b_cand=b_hit=0
    t0=time.time()
    for p in [q for q in pr if plo < q <= phi]:
        M *= p
        bits = M.bit_length()
        nlo = max(2, p.bit_length()+1)
        nhi = bits-2
        for n in range(nlo, nhi+1):
            cA = mpz(1) << n
            NA = M - cA
            a_cand += 1
            if is_prime(NA): a_hit += 1
            # matched control: same bit length, odd part an odd prime > p
            lo, hi = mpz(1)<<(n-2), mpz(1)<<(n-1)
            Pr = next_prime(mpz(random.randint(int(lo), int(hi)-1)))
            if Pr <= p or Pr >= hi: continue
            cB = 2*Pr
            b_cand += 1
            if is_prime(M - cB): b_hit += 1
        if verbose:
            print("p=%-5d digits=%-6d A %d/%d   B %d/%d   %.0fs"
                  %(p, M.num_digits(), a_hit, a_cand, b_hit, b_cand, time.time()-t0))
            sys.stdout.flush()
    return a_cand,a_hit,b_cand,b_hit

if __name__ == "__main__":
    plo, phi, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ac,ah,bc,bh = run(plo, phi, seed)
    print("TOTAL  pow2 %d/%d = %.5f   control %d/%d = %.5f"
          %(ah,ac,ah/ac,bh,bc,bh/bc))
