#pragma once

#include <simdpp/simd.h>

#include "base_smith_waterman.hpp"

using namespace simdpp::arch_avx2;

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
