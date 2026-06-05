#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m=5 test of the tuple-size law. Counts the prime quintuplet
(6N-1, 6N+1, 6N+5, 6N+7, 6N+11)  [pattern (0,2,6,8,12)]  AND the twin
(6N-1, 6N+1) in the SAME segmented sieve pass, so the deep-slope ratio
quintuplet/twin -> 5:2 is self-contained and from scratch.

  rho = (# pattern centres) / (# centres), centres per shell = 1.5e(k-1)
  ell = ln ln(10^k);  target  d ln rho_5 / d ell -> -5 ,  ratio to twin -> 2.5

WHY S11. At S10 the quintuplet has only ~16,570 members (2,936 at S9), so the
S9->S10 ratio is 2.62 +/- 0.09 -- consistent with 2.5 but outside the +/-0.06
bound. S11 yields ~1.0e5 quintuplets; the S10->S11 ratio then has noise ~+/-0.04
and the (deeper, smaller) tail brings it to ~2.51 -- a clean confirmation.

Usage (PowerShell):
    $env:SHELLS="10,11"; python quintuplet_deep.py
S10 ~ minutes; S11 ~ 1-2 hours (pure segmented sieve, no factorization).
Optionally add S12 to pin the count-aware intercept on -5. Requires: numpy.
"""
import os, math, time
import numpy as np

SHELLS = [int(x) for x in os.environ.get("SHELLS", "10,11").split(",")]
SEG    = int(os.environ.get("SEG", 12_000_000))
QUINT  = ((0,-1),(0,1),(1,-1),(1,1),(2,-1))   # value = 6(N+off)+s

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].astype(np.int64)

def shell_counts(k):
    """Return (twin_count, quint_count, centres) for shell S_k, one sieve pass."""
    base = primes_upto(int(math.isqrt(10**k)) + 2)
    Nlo = 10**(k-1)//6 + 1
    Nhi = (10**k - 1)//6
    tw = qn = 0; ns = Nlo
    while ns <= Nhi:
        ne  = min(ns + SEG, Nhi + 1)
        vlo = 6*ns - 1; vhi = 6*(ne - 1) + 11
        comp = np.zeros(vhi - vlo + 1, bool)
        for p in base:
            if p*p > vhi: break
            st = max(p*p, ((vlo + p - 1)//p)*p)
            comp[st - vlo::p] = True
        isp = ~comp
        N = np.arange(ns, ne, dtype=np.int64)
        m1 = isp[6*N - 1 - vlo]; p1 = isp[6*N + 1 - vlo]
        tw += int(np.count_nonzero(m1 & p1))
        ok = m1 & p1
        for off, s in QUINT[2:]:
            ok &= isp[6*N + (6*off + s) - vlo]
        qn += int(np.count_nonzero(ok)); ns = ne
    return tw, qn, (Nhi - Nlo + 1)

def main():
    print("== m=5 quintuplet vs twin, same-pass recompute ==  shells", SHELLS, "\n")
    rows = []
    for k in SHELLS:
        t = time.time(); tw, qn, nc = shell_counts(k)
        ell = math.log(math.log(10.0**k))
        rows.append((k, tw, qn, nc, ell))
        print(f"S{k}: twin={tw}  quint={qn}  centres={nc}  "
              f"ln_rho_twin={math.log(tw/nc):.5f}  ln_rho_quint={math.log(qn/nc):.5f}  "
              f"ell={ell:.5f}  [{time.time()-t:.0f}s]")
    print()
    for (k0,t0,q0,n0,e0),(k1,t1,q1,n1,e1) in zip(rows, rows[1:]):
        s2 = (math.log(t1/n1)-math.log(t0/n0))/(e1-e0)
        s5 = (math.log(q1/n1)-math.log(q0/n0))/(e1-e0)
        sd = math.sqrt(1/q0+1/q1)/(e1-e0)
        print(f"S{k0}->S{k1}:  twin slope {s2:.3f}   quint slope {s5:.3f}+/-{sd:.2f}   "
              f"ratio {s5/s2:.4f}+/-{abs(s5/s2)*sd/abs(s5):.3f}   (target -5 / 2.5)")

if __name__ == "__main__":
    main()
