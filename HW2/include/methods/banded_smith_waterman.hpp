#pragma once

#include "base_smith_waterman.hpp"

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

