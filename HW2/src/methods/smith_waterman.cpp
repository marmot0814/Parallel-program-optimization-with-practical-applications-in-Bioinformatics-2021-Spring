#include <cassert>

#include "methods/smith_waterman.h"

SmithWaterman::SmithWaterman(
  int M, int R, int I, int D
) : M(M), R(R), I(I), D(D) {
  assert(M > 0);
  assert(R < 0);
  assert(I < 0);
  assert(D < 0);
}

BLAST SmithWaterman::operator() (
  const Sequence &seq1,
  const Sequence &seq2
) {
  int len1 = seq1.size(), len2 = seq2.size();
  dp = std::vector<std::vector<int>>(len1 + 1, std::vector<int>(len2 + 1, 0));
  for (int i = 0; i < len1; i++) {
    for (int j = 0; j < len2; j++) {
      if (seq1[i] == seq2[j])
        dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j] + M);
      dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i + 1][j] + I);
      dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j + 1] + D);
      dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j] + R); 
    }
  }

  int maxi = 0, maxj = 0;
  for (int i = 0; i <= len1; i++)
    for (int j = 0; j <= len2; j++)
      if (dp[maxi][maxj] < dp[i][j])
        maxi = i, maxj = j;
  return answer(maxi, maxj, seq1, seq2);
}

BLAST SmithWaterman::answer(
  int i, int j,
  const Sequence &seq1,
  const Sequence &seq2
) {
  int len1 = dp.size() - 1, len2 = dp[0].size() - 1;
  assert(0 <= i && i <= len1);
  assert(0 <= j && j <= len2);
  if (dp[i][j] == 0)
    return BLAST(Sequence(seq1.name, "", i), Sequence(seq2.name, "", j), "");
  if (dp[i][j] == dp[i - 1][j - 1] + M && seq1[i - 1] == seq2[j - 1]) {
    auto ret = answer(i - 1, j - 1, seq1, seq2);
    ret.addM(seq1[i - 1], seq2[j - 1]);
    return ret;
  }
  if (dp[i][j] == dp[i - 1][j - 1] + R) {
    auto ret = answer(i - 1, j - 1, seq1, seq2);
    ret.addR(seq1[i - 1], seq2[j - 1]);
    return ret;
  }
  if (dp[i][j] == dp[i][j - 1] + I) {
    auto ret = answer(i, j - 1, seq1, seq2);
    ret.addI(seq2[j - 1]);
    return ret;
  }
  if (dp[i][j] == dp[i - 1][j] + D) {
    auto ret = answer(i - 1, j, seq1, seq2);
    ret.addD(seq1[i - 1]);
    return ret;
  }
  assert(false);
}
