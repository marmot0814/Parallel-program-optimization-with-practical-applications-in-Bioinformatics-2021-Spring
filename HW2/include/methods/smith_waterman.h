#pragma once

#include <vector>

#include "types/blast.h"
#include "types/sequence.h"

class SmithWaterman {
  int M, R, I, D;
  std::vector<std::vector<int>> dp;
 public:
  SmithWaterman(
    int, int, int, int
  );
  BLAST operator() (
    const Sequence &,
    const Sequence &
  );
  BLAST answer(
    int, int,
    const Sequence &,
    const Sequence &
  );
};
