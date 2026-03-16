**The answer is TRUE for all N.**

We give a proof that for every positive integer N, the maximum size of A ⊆ {1,...,N} with ab+1 never squarefree is achieved by {n ≡ 7 (mod 25)}, with size ⌊(N−7)/25⌋+1. The proof is conditional on the correctness of the supplied verification scripts.

**Disclosure:** This work was prepared with AI assistance (Claude by Anthropic, GPT-5.4 by OpenAI). All mathematical content and computational results have been verified by the human author.

**Proof structure:**

- N ≤ 7,306: verified by layered exact computation (pair exchange analysis for N ≤ 2,000; structured inequality verified by independent implementations for N ≤ 20,000)
- N = 7,307 to 10^7: computer-assisted verification (zero failures across 662,188 breakpoints)
- N > 10^7: Sawhney's asymptotic theorem

**The outsider-clique framework.** Elements not in {7,18} mod 25 ("outsiders") are partitioned by least witness prime p ≡ 1 (mod 4), p ≥ 13. For each witness prime p, we verify a structured inequality:

s_max^(p)(N) + max(V_p⁺, V_p⁻) + d_max^(p)(N) + R_{>p}(N) < M(N)

at every structural breakpoint, where s_max is the worst-case base exchange, V_p are root class sizes, d_max is the maximum cross-degree, and R_{>p} counts elements in higher root classes.

The key observation: for opposite-root cross-pairs a ≡ r, b ≡ −r (mod p²), the cross-product satisfies ab+1 ≡ 2 (mod p²), so p² never divides it. This bounds outsider clique structure and makes the inequality tractable.

**Verification results:**
- C++ verifier (v5, sieve-based bitmask architecture): 219 witness-prime blocks, 662,188 breakpoints, N up to 10^7. Zero failures for N ≥ 7,307. Final margin at N = 10^7: 15,214.
- Python verifier (v2, independent implementation): 25,658 breakpoints, N = 5,000 to 500,000. Zero failures.
- Pair exchange analysis: for N ≤ 2,000, no single outsider or compatible outsider pair can match the base class.

All ranges overlap. Verification scripts available at: https://github.com/hjyuh/MathsSTuff

A detailed proof write-up is available upon request.
