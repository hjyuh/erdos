# Erdős Problem Solutions

Contributions to open problems from the [Erdős Problems](https://www.erdosproblems.com/) collection.

**Author:** Mahmoud  
**AI Assistance:** Claude (Anthropic), GPT-5.4 (OpenAI), Codex (OpenAI)  
**Period:** March 12–15, 2026

---

## Problem 848 — SOLVED ✓

**Statement:** Is the maximum size of A ⊆ {1,...,N} with ab+1 never squarefree achieved by {n ≡ 7 (mod 25)}?

**Answer:** Yes, for all N.

**Method:** Computer-assisted verification (N ≤ 10⁷) via an outsider-clique framework + Sawhney's asymptotic theorem (N > 10⁷).

**Key result:** 219 witness-prime blocks, 662,188 breakpoints, zero failures for N ≥ 7,307. Margin = 15,214 at N = 10⁷.

→ [Full proof](848/proof.md) · [Verifier (C++)](848/erdos848_verifier_v5.cpp) · [Verifier (Python)](848/erdos848_v2.py) · [Certificate](848/erdos848_certificate_v5.tsv)

---

## Problem 388 — Partial Result

**Statement:** Are there only finitely many solutions to f_{k₁}(x) = f_{k₂}(y)?

**Result:** For fixed (k₁, k₂) with both ≥ 4 and k₁ ≠ k₂, only finitely many solutions exist. Proof via Kulkarni-Sury Theorem C.

→ [Full note](388/note.md)

---

## Problem 931 — Partial Result

**Statement (case (4,3)):** Are there infinitely many pairs with equal prime support?

**Results:** Local finiteness (fixed n₁) and gap-fixed finiteness proved. 26 valid pairs found computationally, all with n₁ ≤ 636.

→ [Full note](931/note.md)

---

## Reproduction

### Problem 848

```bash
# Compile the verifier
g++ -O3 -march=native -std=c++20 -DNDEBUG 848/erdos848_verifier_v5.cpp -o erdos848_v5

# Run (takes ~10 hours on a modern CPU)
./erdos848_v5 -n 10000000

# Check results
grep "FAIL" erdos848_certificate_v5.tsv | tail -3
# Expected: last FAIL at N = 7306

# Independent verification (Python, covers N = 5000 to 500000)
python3 848/erdos848_v2.py 500000
```

---

## Disclosure

All work conducted with AI assistance. Claude (Anthropic, Opus 4.6) provided analysis, code, and orchestration. GPT-5.4 (OpenAI) provided mathematical proofs, algorithm design, and adversarial review. Codex (OpenAI) implemented multiple verifier iterations. All mathematical content and computational results verified by the human author.
