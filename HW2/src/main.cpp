#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "types/table/table.h"
#include "types/table/striped_table.h"
#include "types/sequence/sequence.h"
#include "types/sequence/striped_sequence.h"
#include "types/fasta.h"
#include "types/blast.h"

#include <simdpp/simd.h>

using namespace simdpp::arch_avx2;

class BaseSmithWaterman {
 protected:
  int match, mismatch, indel;

  int c2i(char c) const {
    switch (c) {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      default: assert(false);
    }
  }

  char i2c(int i) const {
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
  std::vector<std::vector<int>> prepare_match_mismatch(const StripedSequence &Q) {
    std::vector<std::vector<int>> W(4, std::vector<int>(Q.size()));
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < (int)Q.size(); j++)
        W[i][j] = (Q.striped_get(j) == i2c(i) ? match : mismatch);
    return W;
  }

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
    auto match_mismatch = prepare_match_mismatch(Q);

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

      int *match_mismatch_ptr = match_mismatch[c2i(D[i])].data();
      for (int j = 1; j < t; j++) {
        int32<p> all_match_mismatch = load_u(match_mismatch_ptr + j * p);
        reg[j] = max(reg[j], prev[j - 1] + all_match_mismatch);
      }

      auto reg_right_shift = [&](const int32<p> &reg) {
        std::vector<int> tmp(p + 1, 0);
        store_u(tmp.data() + 1, reg);
        int32<p> ret = load_u(tmp.data());
        return ret;
      };

      reg[0] = max(reg[0], reg_right_shift(prev[t - 1]) + int32<p>(load_u(match_mismatch_ptr)));

      while (1) {
        for (int j = 1; j < t; j++)
          reg[j] = max(reg[j], reg[j - 1] + all_indel);
        int32<p> checker = reg_right_shift(reg[t - 1]), checked = reg[0];
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
    int maxi = 0, maxj = 0, maxv = 0;
    auto dp = Table(D.size(), Q.size());
    for (int i = 0; i < (int)D.size(); i++) {
      for (int j = 0; j < (int)Q.size(); j++) {
        int curv = std::max(dp.get(i, j - 1), dp.get(i - 1, j)) + indel;
        if (D.get(i) == Q.get(j))
          curv = std::max(curv, dp.get(i - 1, j - 1) + match);
        else
          curv = std::max(curv, dp.get(i - 1, j - 1) + mismatch);
        if (curv > maxv)
          maxi = i, maxj = j;
        *dp.ptr(i, j) = curv;
      }
    }
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
