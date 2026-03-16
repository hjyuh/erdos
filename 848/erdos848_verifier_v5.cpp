// erdos848_verifier_v5.cpp — Optimized with cross-masks and 64-bit modinv
// Based on v4. Key changes:
// 1) 64-bit modinv (5-10x faster than 128-bit for our prime range)
// 2) Sieve-based cross-compatibility masks (eliminates trial-division bottleneck)
// 3) All other v4 features preserved (breakpoint sweep, incremental cert, resume)

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using i64 = long long;
using u64 = std::uint64_t;
using i128 = __int128_t;

struct Options {
  int n_max = 10'000'000;
  std::string certificate_path = "erdos848_certificate_v5.tsv";
  bool quiet = false;
  int resume_n = 0;
};

struct RootInfo {
  int p = 0;
  int p2 = 0;
  int r_plus = 0;
  int r_minus = 0;
};

struct RootEvent {
  int n = 0;
  int p_index = -1;
};

struct BitMatrix {
  int rows = 0;
  int cols = 0;
  int words_per_row = 0;
  std::vector<u64> words;

  BitMatrix() = default;

  BitMatrix(int row_count, int col_count)
      : rows(row_count), cols(col_count),
        words_per_row((col_count + 63) / 64),
        words(static_cast<std::size_t>(row_count) *
              static_cast<std::size_t>((col_count + 63) / 64), 0) {}

  u64 *row_ptr(int row) {
    return words.data() +
           static_cast<std::size_t>(row) *
               static_cast<std::size_t>(words_per_row);
  }
  const u64 *row_ptr(int row) const {
    return words.data() +
           static_cast<std::size_t>(row) *
               static_cast<std::size_t>(words_per_row);
  }

  void set_bit(int row, int col) {
    row_ptr(row)[static_cast<std::size_t>(col) >> 6] |=
        (u64{1} << (col & 63));
  }

  bool test_bit(int row, int col) const {
    return ((row_ptr(row)[static_cast<std::size_t>(col) >> 6] >> (col & 63)) &
            u64{1}) != 0;
  }

  int prefix_popcount(int row, int length) const {
    if (length <= 0) return 0;
    if (length > cols) length = cols;
    const u64 *ptr = row_ptr(row);
    const int full_words = length >> 6;
    const int trailing_bits = length & 63;
    int total = 0;
    for (int i = 0; i < full_words; ++i)
      total += std::popcount(ptr[i]);
    if (trailing_bits != 0)
      total += std::popcount(ptr[full_words] & ((u64{1} << trailing_bits) - 1));
    return total;
  }

  // Scan row, call func(col) for each set bit in [0, length)
  template <typename F>
  void scan_row(int row, int length, F &&func) const {
    if (length <= 0) return;
    if (length > cols) length = cols;
    const u64 *ptr = row_ptr(row);
    const int full_words = length >> 6;
    const int trailing_bits = length & 63;
    for (int w = 0; w < full_words; ++w) {
      u64 bits = ptr[w];
      while (bits) {
        int bit = __builtin_ctzll(bits);
        func(w * 64 + bit);
        bits &= bits - 1;
      }
    }
    if (trailing_bits != 0) {
      u64 bits = ptr[full_words] & ((u64{1} << trailing_bits) - 1);
      while (bits) {
        int bit = __builtin_ctzll(bits);
        func(full_words * 64 + bit);
        bits &= bits - 1;
      }
    }
  }
};

struct PrimeBlockData {
  RootInfo root;
  std::vector<int> plus_values;
  std::vector<int> minus_values;
  BitMatrix plus_base_masks;  // plus outsiders × base elements
  BitMatrix minus_base_masks; // minus outsiders × base elements
  BitMatrix cross_pm;         // plus outsiders × minus outsiders (compatibility)
  BitMatrix cross_mp;         // minus outsiders × plus outsiders (compatibility)
};

bool is_base_residue(int x) {
  const int r = x % 25;
  return r == 7 || r == 18;
}

// Fast 64-bit modular inverse (for moduli < 2^62)
std::optional<i64> modinv64(i64 a, i64 m) {
  if (m == 1) return 0;
  a %= m;
  if (a < 0) a += m;
  if (a == 0) return std::nullopt;
  i64 t = 0, new_t = 1;
  i64 r = m, new_r = a;
  while (new_r != 0) {
    i64 q = r / new_r;
    i64 tmp_t = t - q * new_t;
    t = new_t;
    new_t = tmp_t;
    i64 tmp_r = r - q * new_r;
    r = new_r;
    new_r = tmp_r;
  }
  if (r != 1) return std::nullopt;
  if (t < 0) t += m;
  return t;
}

i64 mul_mod(i64 a, i64 b, i64 mod) {
  return static_cast<i64>((static_cast<i128>(a) * static_cast<i128>(b)) %
                           static_cast<i128>(mod));
}

std::vector<int> sieve_primes(int limit) {
  std::vector<bool> is_prime(limit + 1, true);
  if (limit >= 0) is_prime[0] = false;
  if (limit >= 1) is_prime[1] = false;
  for (int i = 2; 1LL * i * i <= limit; ++i) {
    if (!is_prime[i]) continue;
    for (int j = i * i; j <= limit; j += i) is_prime[j] = false;
  }
  std::vector<int> primes;
  for (int i = 2; i <= limit; ++i)
    if (is_prime[i]) primes.push_back(i);
  return primes;
}

RootInfo hensel_root_prime_square(int p) {
  int r0 = -1;
  for (int r = 1; r < p; ++r) {
    if ((1LL * r * r + 1) % p == 0) { r0 = r; break; }
  }
  if (r0 < 0) throw std::runtime_error("no root mod p");
  const i64 p2 = 1LL * p * p;
  const i64 t = ((1LL * r0 * r0 + 1) / p) % p;
  auto inv = modinv64((2LL * r0) % p, p);
  if (!inv) throw std::runtime_error("hensel invert failed");
  i64 k = (-t * *inv) % p;
  if (k < 0) k += p;
  i64 r = (r0 + k * p) % p2;
  if (r <= 0) r += p2;
  if (((r * r) + 1) % p2 != 0)
    throw std::runtime_error("hensel verify failed");
  return {p, (int)p2, (int)r, (int)(p2 - r)};
}

std::vector<RootInfo> build_root_infos(int n_max) {
  const int limit = (int)std::sqrt((double)n_max) + 5;
  const auto primes = sieve_primes(limit);
  std::vector<RootInfo> roots;
  for (int p : primes) {
    if (p < 13 || p % 4 != 1) continue;
    if (1LL * p * p > n_max) break;
    roots.push_back(hensel_root_prime_square(p));
  }
  return roots;
}

int least_witness_prime(int x, const std::vector<RootInfo> &roots,
                        std::size_t up_to) {
  for (std::size_t i = 0; i <= up_to; ++i) {
    int mod = x % roots[i].p2;
    if (mod == roots[i].r_plus || mod == roots[i].r_minus) return roots[i].p;
  }
  return 0;
}

std::vector<int> build_witness_values(int residue, int p2, int n_max,
                                       const std::vector<RootInfo> &roots,
                                       std::size_t p_index) {
  std::vector<int> vals;
  for (int x = residue; x <= n_max; x += p2) {
    if (is_base_residue(x)) continue;
    if (least_witness_prime(x, roots, p_index) == roots[p_index].p)
      vals.push_back(x);
  }
  return vals;
}

std::vector<int> build_mask_primes(int n_max) {
  const i64 m_max = (n_max >= 7) ? ((n_max - 7) / 25 + 1LL) : 0LL;
  const i64 y_max = (m_max > 0) ? (7LL + 25LL * (m_max - 1)) : 0LL;
  const i64 max_val = (i64)n_max * y_max + 1;
  const int limit = (int)std::sqrt((double)max_val) + 5;
  return sieve_primes(limit);
}

// Build cross-mask primes: need primes up to sqrt(n_max^2)=n_max
std::vector<int> build_cross_primes(int n_max) {
  return sieve_primes(n_max);
}

void mark_progression_bits(u64 *row_words, int cols, i64 start, i64 step) {
  for (i64 idx = start; idx < cols; idx += step)
    row_words[static_cast<std::size_t>(idx) >> 6] |= (u64{1} << (idx & 63));
}

// Build base-exchange masks (outsider × base-class compatibility)
BitMatrix build_base_masks(const std::vector<int> &values, int m_max,
                            const std::vector<int> &mask_primes, bool quiet,
                            int p, const char *label) {
  BitMatrix mat((int)values.size(), m_max);
  if (values.empty() || m_max == 0) return mat;
  std::atomic<int> done{0};

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
  for (int row = 0; row < (int)values.size(); ++row) {
    const i64 x = values[row];
    const i64 a = 25LL * x;
    const i64 b = 7LL * x + 1;
    const i64 max_val = x * (7LL + 25LL * (m_max - 1)) + 1;
    u64 *rw = mat.row_ptr(row);
    for (int q : mask_primes) {
      const i64 qq = 1LL * q * q;
      if (qq > max_val) break;
      const i64 g = std::gcd(a, qq);
      if (b % g != 0) continue;
      const i64 a1 = a / g, mod = qq / g, b1 = b / g;
      auto inv = modinv64(a1 % mod, mod);
      if (!inv) continue;
      i64 r = ((-b1) % mod + mod) % mod;
      r = mul_mod(r, *inv, mod);
      mark_progression_bits(rw, m_max, r, mod);
    }
    if (!quiet) {
      int d = done.fetch_add(1) + 1;
      if (d % 10000 == 0 || d == (int)values.size()) {
#ifdef _OPENMP
#pragma omp critical
#endif
        std::cerr << "[p=" << p << "] base masks " << label << ": " << d
                  << "/" << values.size() << " ("
                  << (100 * d / (int)values.size()) << "%)\n";
      }
    }
  }
  if (!quiet)
    std::cerr << "Built " << label << " base masks p=" << p << " rows="
              << values.size() << " cols=" << m_max << "\n";
  return mat;
}

// Build cross-compatibility masks using sieve (the key optimization)
// cross[i][j] = 1 iff from_values[i] * to_values[j] + 1 is non-squarefree
void build_cross_masks(BitMatrix &cross, const std::vector<int> &from_vals,
                       const std::vector<int> &to_vals,
                       const std::vector<int> &cross_primes, int n_max,
                       bool quiet, int p, const char *label) {
  const int R = (int)from_vals.size();
  const int C = (int)to_vals.size();
  cross = BitMatrix(R, C);
  if (R == 0 || C == 0) return;

  // Build a map from to-value to its index
  std::unordered_map<int, int> to_index;
  to_index.reserve(C);
  for (int j = 0; j < C; ++j) to_index[to_vals[j]] = j;

  int primes_done = 0;
  const int total_primes = (int)cross_primes.size();

  for (int qi = 0; qi < total_primes; ++qi) {
    const int q = cross_primes[qi];
    const i64 qq = 1LL * q * q;

    // Skip if qq > max possible cross product
    if (qq > 1LL * n_max * n_max + 1) break;

    if (qq <= C) {
      // Phase 1: small primes — build buckets
      std::vector<std::vector<int>> buckets(qq);
      for (int j = 0; j < C; ++j)
        buckets[to_vals[j] % qq].push_back(j);

      for (int i = 0; i < R; ++i) {
        i64 x_mod = from_vals[i] % qq;
        auto inv = modinv64(x_mod, qq);
        if (!inv) continue;
        i64 target = (qq - *inv) % qq;
        for (int j : buckets[target])
          cross.set_bit(i, j);
      }
    } else if (qq <= n_max) {
      // Phase 2: medium primes — iterate progression + hash lookup
      for (int i = 0; i < R; ++i) {
        i64 x_mod = from_vals[i] % qq;
        auto inv = modinv64(x_mod, qq);
        if (!inv) continue;
        i64 target = (qq - *inv) % qq;
        // Find all to_vals ≡ target (mod qq) up to n_max
        for (i64 y = target; y <= n_max; y += qq) {
          auto it = to_index.find((int)y);
          if (it != to_index.end())
            cross.set_bit(i, it->second);
        }
      }
    } else {
      // Phase 3: large primes (qq > n_max) — at most one match per from-value
      for (int i = 0; i < R; ++i) {
        i64 x_mod = from_vals[i] % qq;
        auto inv = modinv64(x_mod, qq);
        if (!inv) continue;
        i64 target = (qq - *inv) % qq;
        if (target <= n_max) {
          auto it = to_index.find((int)target);
          if (it != to_index.end())
            cross.set_bit(i, it->second);
        }
      }
    }

    ++primes_done;
    if (!quiet && (primes_done % 50000 == 0 || qi + 1 == total_primes)) {
      std::cerr << "[p=" << p << "] cross " << label << ": prime " << primes_done
                << "/" << total_primes << " (" << (100 * primes_done / total_primes)
                << "%) q=" << q << "\n";
    }
  }

  if (!quiet)
    std::cerr << "Built cross masks " << label << " p=" << p << " " << R
              << "x" << C << "\n";
}

PrimeBlockData build_block(int n_max, int m_max,
                            const std::vector<RootInfo> &roots,
                            std::size_t p_index,
                            const std::vector<int> &mask_primes,
                            const std::vector<int> &cross_primes,
                            bool quiet) {
  PrimeBlockData b;
  b.root = roots[p_index];
  b.plus_values = build_witness_values(b.root.r_plus, b.root.p2, n_max, roots, p_index);
  b.minus_values = build_witness_values(b.root.r_minus, b.root.p2, n_max, roots, p_index);

  // Base masks
  b.plus_base_masks = build_base_masks(b.plus_values, m_max, mask_primes, quiet, b.root.p, "+");
  b.minus_base_masks = build_base_masks(b.minus_values, m_max, mask_primes, quiet, b.root.p, "-");

  // Cross masks (sieve-based, replaces trial division in sweep)
  build_cross_masks(b.cross_pm, b.plus_values, b.minus_values, cross_primes,
                    n_max, quiet, b.root.p, "+->-");
  build_cross_masks(b.cross_mp, b.minus_values, b.plus_values, cross_primes,
                    n_max, quiet, b.root.p, "-->+");

  return b;
}

std::vector<RootEvent> build_global_root_events(int n_max,
                                                 const std::vector<RootInfo> &roots) {
  std::vector<RootEvent> events;
  for (int pidx = 0; pidx < (int)roots.size(); ++pidx) {
    const auto &rt = roots[pidx];
    for (int x = rt.r_plus; x <= n_max; x += rt.p2)
      events.push_back({x, pidx});
    for (int x = rt.r_minus; x <= n_max; x += rt.p2)
      events.push_back({x, pidx});
  }
  std::sort(events.begin(), events.end(),
            [](const RootEvent &a, const RootEvent &b) {
              return a.n != b.n ? a.n < b.n : a.p_index < b.p_index;
            });
  return events;
}

std::vector<int> build_breakpoints(int n_max, int m_max,
                                    const std::vector<RootEvent> &root_events) {
  std::vector<int> bp;
  bp.reserve(m_max + root_events.size() + 8);
  for (int m = 0; m < m_max; ++m) bp.push_back(7 + 25 * m);
  for (const auto &ev : root_events)
    if (ev.n >= 2 && ev.n <= n_max) bp.push_back(ev.n);
  std::sort(bp.begin(), bp.end());
  bp.erase(std::unique(bp.begin(), bp.end()), bp.end());
  return bp;
}

void write_cert_interval(std::ofstream &out, int start, int end, int margin,
                          int wp) {
  for (int n = start; n <= end; ++n) {
    int bc = (n >= 7) ? ((n - 7) / 25 + 1) : 0;
    out << n << '\t' << ((margin > 0) ? "PASS" : "FAIL") << '\t' << margin
        << '\t' << wp << '\t' << bc << '\n';
    if (n % 10000 == 0 || n == end) out.flush();
  }
}

int read_last_n(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("resume: no file");
  std::string line;
  int last = -1;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] < '0' || line[0] > '9') continue;
    auto tab = line.find('\t');
    if (tab != std::string::npos) last = std::stoi(line.substr(0, tab));
  }
  if (last < 0) throw std::runtime_error("resume: no data");
  return last;
}

Options parse_options(int argc, char **argv) {
  Options opts;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-n" || arg == "--n-max") && i + 1 < argc)
      opts.n_max = std::stoi(argv[++i]);
    else if ((arg == "-o" || arg == "--output") && i + 1 < argc)
      opts.certificate_path = argv[++i];
    else if (arg == "--resume" && i + 1 < argc)
      opts.resume_n = std::stoi(argv[++i]);
    else if (arg == "-q" || arg == "--quiet")
      opts.quiet = true;
    else
      throw std::runtime_error("unknown argument: " + arg);
  }
  return opts;
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Options opts = parse_options(argc, argv);
    if (opts.n_max < 2) throw std::runtime_error("n_max must be >= 2");

    const int m_max = (opts.n_max >= 7) ? ((opts.n_max - 7) / 25 + 1) : 0;
    if (!opts.quiet) std::cerr << "Building primes...\n";
    const auto mask_primes = build_mask_primes(opts.n_max);
    const auto cross_primes = build_cross_primes(opts.n_max);
    const auto roots = build_root_infos(opts.n_max);
    const auto global_root_events = build_global_root_events(opts.n_max, roots);
    const auto breakpoints = build_breakpoints(opts.n_max, m_max, global_root_events);

    if (!opts.quiet)
      std::cerr << "Primes: " << roots.size() << " witness blocks, "
                << breakpoints.size() << " breakpoints\n";

    const auto resume_bp =
        std::upper_bound(breakpoints.begin(), breakpoints.end(), opts.resume_n);
    const std::size_t resume_idx = resume_bp - breakpoints.begin();

    std::vector<int> best_margin(breakpoints.size() + 1,
                                  std::numeric_limits<int>::max() / 4);
    std::vector<int> best_p(breakpoints.size() + 1, 0);

    for (std::size_t pi = 0; pi < roots.size(); ++pi) {
      const auto &root = roots[pi];
      if (!opts.quiet)
        std::cerr << "\nPreparing block p=" << root.p << " (" << (pi + 1)
                  << "/" << roots.size() << ")...\n";

      PrimeBlockData block =
          build_block(opts.n_max, m_max, roots, pi, mask_primes, cross_primes,
                      opts.quiet);

      if (!opts.quiet)
        std::cerr << "Sweeping p=" << root.p << " +=" << block.plus_values.size()
                  << " -=" << block.minus_values.size() << "\n";

      std::vector<int> pex(block.plus_values.size(), 0);
      std::vector<int> mex(block.minus_values.size(), 0);
      std::vector<int> pdeg(block.plus_values.size(), 0);
      std::vector<int> mdeg(block.minus_values.size(), 0);

      int raw_plus = 0, raw_minus = 0, smax = 0, dmax = 0;
      int active = 0, r_gt = 0, base_count = 0;
      std::size_t np_out = 0, nm_out = 0, ng_root = 0;
      int base_idx = 0;
      int next_base = (m_max > 0) ? 7 : std::numeric_limits<int>::max();
      int next_plus_raw = root.r_plus;
      int next_minus_raw = root.r_minus;
      if (next_plus_raw > opts.n_max)
        next_plus_raw = std::numeric_limits<int>::max();
      if (next_minus_raw > opts.n_max)
        next_minus_raw = std::numeric_limits<int>::max();

      int progress_counter = 0;

      for (std::size_t bi = 0; bi < breakpoints.size(); ++bi) {
        const int n = breakpoints[bi];

        // Base event
        if (n == next_base) {
          ++base_count;
          int bidx = base_count - 1;
          for (std::size_t i = 0; i < np_out; ++i) {
            if (block.plus_base_masks.test_bit((int)i, bidx)) {
              ++pex[i];
              smax = std::max(smax, pex[i]);
            }
          }
          for (std::size_t i = 0; i < nm_out; ++i) {
            if (block.minus_base_masks.test_bit((int)i, bidx)) {
              ++mex[i];
              smax = std::max(smax, mex[i]);
            }
          }
          ++base_idx;
          next_base = (base_idx < m_max) ? (7 + 25 * base_idx)
                                          : std::numeric_limits<int>::max();
        }

        // Plus outsider event — use cross masks instead of trial division
        if (n == next_plus_raw) {
          ++raw_plus;
          while (np_out < block.plus_values.size() &&
                 block.plus_values[np_out] == n) {
            int row = (int)np_out;
            pex[np_out] =
                block.plus_base_masks.prefix_popcount(row, base_count);
            smax = std::max(smax, pex[np_out]);

            // Cross degree from precomputed mask
            int degree = block.cross_pm.prefix_popcount(row, (int)nm_out);
            pdeg[np_out] = degree;
            dmax = std::max(dmax, degree);

            // Update minus degrees by scanning this row's cross mask
            block.cross_pm.scan_row(row, (int)nm_out, [&](int j) {
              ++mdeg[j];
              dmax = std::max(dmax, mdeg[j]);
            });

            ++active;
            ++np_out;
          }
          next_plus_raw += root.p2;
          if (next_plus_raw > opts.n_max)
            next_plus_raw = std::numeric_limits<int>::max();
        }

        // Minus outsider event — use cross masks instead of trial division
        if (n == next_minus_raw) {
          ++raw_minus;
          while (nm_out < block.minus_values.size() &&
                 block.minus_values[nm_out] == n) {
            int row = (int)nm_out;
            mex[nm_out] =
                block.minus_base_masks.prefix_popcount(row, base_count);
            smax = std::max(smax, mex[nm_out]);

            // Cross degree from precomputed mask
            int degree = block.cross_mp.prefix_popcount(row, (int)np_out);
            mdeg[nm_out] = degree;
            dmax = std::max(dmax, degree);

            // Update plus degrees
            block.cross_mp.scan_row(row, (int)np_out, [&](int i) {
              ++pdeg[i];
              dmax = std::max(dmax, pdeg[i]);
            });

            ++active;
            ++nm_out;
          }
          next_minus_raw += root.p2;
          if (next_minus_raw > opts.n_max)
            next_minus_raw = std::numeric_limits<int>::max();
        }

        // R_{>p} events
        while (ng_root < global_root_events.size() &&
               global_root_events[ng_root].n == n) {
          if (global_root_events[ng_root].p_index > (int)pi) ++r_gt;
          ++ng_root;
        }

        int margin;
        if (active == 0) {
          margin = std::numeric_limits<int>::max() / 4;
        } else {
          int vp = std::max(raw_plus, raw_minus);
          margin = base_count - (smax + vp + dmax + r_gt);
        }

        if (bi + 1 >= resume_idx && margin < best_margin[bi + 1]) {
          best_margin[bi + 1] = margin;
          best_p[bi + 1] = root.p;
        }

        if (n > opts.resume_n) {
          ++progress_counter;
          if (!opts.quiet && (progress_counter % 10000 == 0 ||
                              bi + 1 == breakpoints.size())) {
            std::cerr << "[p=" << root.p << "] bp " << progress_counter
                      << " N=" << n << " margin=" << margin << "\n";
          }
        }
      }
    }

    // Write certificate
    const bool resuming = opts.resume_n > 0;
    if (resuming && read_last_n(opts.certificate_path) != opts.resume_n)
      throw std::runtime_error("resume mismatch");

    std::ofstream out(opts.certificate_path,
                      std::ios::binary |
                          (resuming ? std::ios::app : std::ios::trunc));
    if (!out) throw std::runtime_error("can't open output");
    if (!resuming) { out << "N\tstatus\tmargin\tworst_p\tM\n"; out.flush(); }

    const int write_start = std::max(2, opts.resume_n + 1);
    if (write_start <= opts.n_max && !breakpoints.empty()) {
      std::size_t si = (std::size_t)(std::upper_bound(breakpoints.begin(),
                                                       breakpoints.end(),
                                                       opts.resume_n) -
                                      breakpoints.begin());
      while (si <= breakpoints.size()) {
        int seg_start = (si == 0) ? 2 : breakpoints[si - 1];
        int seg_end = (si < breakpoints.size()) ? breakpoints[si] - 1 : opts.n_max;
        int istart = std::max(write_start, seg_start);
        if (istart <= seg_end && istart <= opts.n_max) {
          int margin, wp;
          if (si == 0) { margin = 1; wp = 0; }
          else {
            margin = (best_margin[si] >= std::numeric_limits<int>::max() / 8)
                         ? ((istart >= 7 ? (istart - 7) / 25 + 1 : 0) + 1)
                         : best_margin[si];
            wp = best_p[si];
          }
          write_cert_interval(out, istart, std::min(seg_end, opts.n_max),
                              margin, wp);
        }
        if (si == breakpoints.size()) break;
        ++si;
      }
    }

    int fbc = (opts.n_max >= 7) ? ((opts.n_max - 7) / 25 + 1) : 0;
    int fs = breakpoints.empty() ? 0 : (int)breakpoints.size();
    int fm = (best_margin[fs] >= std::numeric_limits<int>::max() / 8)
                 ? (fbc + 1)
                 : best_margin[fs];
    int fp = best_p[fs];

    out.flush();
    std::cout << "Certificate written to " << opts.certificate_path << "\n";
    std::cout << "Final base count M(" << opts.n_max << ") = " << fbc << "\n";
    std::cout << "Worst active prime at N=" << opts.n_max << ": p=" << fp
              << ", margin=" << fm << "\n";
    return (fm > 0) ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
