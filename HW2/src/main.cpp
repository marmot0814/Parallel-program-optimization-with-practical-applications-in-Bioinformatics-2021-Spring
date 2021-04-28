#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "types/sequence.h"
#include "types/fasta.h"
#include "types/blast.h"

#include <simdpp/simd.h>

using namespace simdpp::arch_avx2;

class Table : public std::vector<std::vector<int>> {
 protected:
  int n, m;
 public:
  Table(
    int n, int m
  ) : std::vector<std::vector<int>>(n, std::vector<int>(m, 0)), n(n), m(m) {
    assert(n > 0 && m > 0);
  }
  int get(int i, int j) const {
    if (i < 0 || j < 0) return 0;
    return (*this)[i][j];
  }
  int* ptr(int i, int j) {
    assert(0 <= i && i < n);
    assert(0 <= j && j < m);
    return (*this)[i].data() + j;
  };
};

class StripedTable : public Table {
  int striped(int i) const {
    if (i < 0) return i;
    const int p = SIMDPP_FAST_INT32_SIZE;
    const int t = m / p;
    return i % t * p + i / t;
  }

 public:
  StripedTable(int n, int m)
    : Table(n, (m + SIMDPP_FAST_INT32_SIZE - 1) / SIMDPP_FAST_INT32_SIZE * SIMDPP_FAST_INT32_SIZE) {}

  int get(int i, int j) const {
    return Table::get(i, striped(j));
  }
  int* ptr(int i, int j) {
    return Table::ptr(i, j);
  }
};

struct StripedSequence : public Sequence {

  int striped(int i) const {
    if (i < 0) return i;
    const int p = SIMDPP_FAST_INT32_SIZE;
    const int t = size() / p;
    return i % p * t + i / p;
  }

 public:
  StripedSequence(const Sequence &seq) : Sequence(seq) {}
   
  char striped_get(int i) const {
    return Sequence::get(striped(i));
  }

  int size() const {
    const int p = SIMDPP_FAST_INT32_SIZE;
    return (Sequence::size() + p - 1) / p * p;
  }
};

class BaseSmithWaterman {
 protected:
  int match, mismatch, indel;

  int c2i(char c) {
    switch (c) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      default: assert(false);
    }
  }
  char i2c(int i) {
    switch (i) {
      case 0: return 'A';
      case 1: return 'C';
      case 2: return 'G';
      case 3: return 'T';
      default: assert(false);
    }
  }

 public:
  BaseSmithWaterman(
    int match, int mismatch, int indel
  ) : match(match), mismatch(mismatch), indel(indel) {
    assert(match > 0);
    assert(mismatch < 0);
    assert(indel < 0);
  }

  BLAST backtrack(
    int i, int j,
    const Sequence &D,
    const auto &Q,
    const auto &dp
  ) {
    if (dp.get(i, j) == 0)
      return BLAST(Sequence(D.name, ""), Sequence(Q.name, ""), "");

    if (dp.get(i, j) == dp.get(i - 1, j - 1) + match && D.get(i) == Q.get(j)) {
      auto ret = backtrack(i - 1, j - 1, D, Q, dp);
      ret.addM(D.get(i), Q.get(j));
      return ret;
    }

    if (dp.get(i, j) == dp.get(i - 1, j - 1) + mismatch) {
      auto ret = backtrack(i - 1, j - 1, D, Q, dp);
      ret.addR(D.get(i), Q.get(j));
      return ret;
    }
    if (dp.get(i, j) == dp.get(i, j - 1) + indel) {
      auto ret = backtrack(i, j - 1, D, Q, dp);
      ret.addI(Q.get(j));
      return ret;
    }
    if (dp.get(i, j) == dp.get(i - 1, j) + indel) {
      auto ret = backtrack(i - 1, j, D, Q, dp);
      ret.addD(D.get(i));
      return ret;
    }
    assert(false);
  }
};

class StripedSmithWaterman : public BaseSmithWaterman {
 public:
  StripedSmithWaterman(
    int match, int mismatch, int indel
  ) : BaseSmithWaterman(match, mismatch, indel) {}

  BLAST operator() (
    const Sequence &D,
    const Sequence &_Q
  ) {
    const auto Q = StripedSequence(_Q);

    auto dp = StripedTable(D.size(), Q.size());
  
    std::vector<std::vector<int>> W(4, std::vector<int>(Q.size()));
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < (int)Q.size(); j++)
        W[i][j] = (Q.striped_get(j) == i2c(i) ? match : mismatch);

    const int p = SIMDPP_FAST_INT32_SIZE;
    const int t = Q.size() / p;

    int32<p> all_zero = make_int(0);
    int32<p> all_indel = make_int(indel);
    std::vector<int32<p>> reg(t, all_zero);
    for (int i = 0; i < (int)D.size(); i++) {
      auto prev = std::vector<int32<p>>(t, all_zero);
      swap(prev, reg);

      for (int j = 0; j < t; j++)
        reg[j] = max(reg[j], prev[j] + all_indel);

      int *W_ptr = W[c2i(D[i])].data();
      for (int j = 1; j < t; j++) {
        int32<p> w = load_u(W_ptr + j * p);
        reg[j] = max(reg[j], prev[j - 1] + w);
      }

      int32<p> w = load_u(W_ptr);
      std::vector<int> tmp(p + 1, 0);
      store_u(tmp.data() + 1, prev[t - 1]);
      int32<p> checker = load_u(tmp.data());
      reg[0] = max(reg[0], checker + w);

      while (1) {
        for (int j = 1; j < t; j++)
          reg[j] = max(reg[j], reg[j - 1] + all_indel);

        std::vector<int> tmp(p + 1, 0);
        store_u(tmp.data() + 1, reg[t - 1]);
        int32<p> checker = load_u(tmp.data());

        int32<p> checked = reg[0];

        if (not reduce_or(max(all_zero, checker + all_indel - checked)))
          break;

        reg[0] = max(reg[0], checker + all_indel);
      }

      for (int j = 0; j < t; j++)
        store_u(dp.ptr(i, j * p), reg[j]);
    }

    int maxi = 0, maxv = 0;
    for (int i = 0; i < (int)D.size(); i++) {
      int cur = 0;
      for (int j = 0; j < t; j++) {
        int32<p> x = load_u(dp.ptr(i, j * p));
        cur = std::max(cur, reduce_max(x));
      }
      if (cur > maxv)
        maxi = i, maxv = cur;
    }

    int maxj = 0;
    for (int j = 0; j < (int)Q.size(); j++)
      if (dp.get(maxi, j) > dp.get(maxi, maxj))
        maxj = j;

    return backtrack(maxi, maxj, D, Q, dp);
  }
};

class BandedSmithWaterman : public BaseSmithWaterman {
 public:
  BandedSmithWaterman(
    int match, int mismatch, int indel
  ) : BaseSmithWaterman(match, mismatch, indel) {}

  BLAST operator() (
    const Sequence &D,
    const Sequence &Q
  ) {
    auto dp = Table(D.size(), Q.size());
    for (int i = 0; i < (int)D.size(); i++) {
      for (int j = 0; j < (int)Q.size(); j++) {
        *dp.ptr(i, j) = std::max(dp.get(i, j), dp.get(i - 1, j - 1) + (D.get(i) == Q.get(j) ? match : mismatch));
        *dp.ptr(i, j) = std::max(dp.get(i, j), std::max(dp.get(i, j - 1), dp.get(i - 1, j)) + indel);
      }
    }

    int maxi = 0, maxj = 0;
    for (int i = 0; i < (int)D.size(); i++)
      for (int j = 0; j < (int)Q.size(); j++)
        if (dp.get(i, j) > dp.get(maxi, maxj))
          maxi = i, maxj = j;

    return backtrack(maxi, maxj, D, Q, dp);
  }
};

#include <chrono>
int main() {
  FASTA fasta;
  std::cin >> fasta;

  auto seq_bg = std::chrono::high_resolution_clock::now();
  auto sequential = BandedSmithWaterman(3, -3, -2)(fasta[0], fasta[1]);
  auto seq_ed = std::chrono::high_resolution_clock::now();
  auto seq_time = std::chrono::duration_cast<std::chrono::duration<double>>(seq_ed - seq_bg).count();

  auto simd_bg = std::chrono::high_resolution_clock::now();
  auto simd = StripedSmithWaterman(3, -3, -2)(fasta[0], fasta[1]);
  auto simd_ed = std::chrono::high_resolution_clock::now();
  auto simd_time = std::chrono::duration_cast<std::chrono::duration<double>>(simd_ed - simd_bg).count();

  std::cout << "seq : " << seq_time << '\n';
  std::cout << sequential << '\n';
  std::cout << "simd: " << simd_time << '\n';
  std::cout << simd << '\n'; 
  std::cout << "speed up: " << seq_time / simd_time << "x\n";
}
