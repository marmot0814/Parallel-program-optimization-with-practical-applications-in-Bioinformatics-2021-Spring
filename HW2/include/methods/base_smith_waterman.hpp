#pragma once

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
