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

#include "methods/striped_smith_waterman.hpp"
#include "methods/banded_smith_waterman.hpp"

#include "utils/clock.h"

#include <chrono>
int main() {
  FASTA fasta;
  std::cin >> fasta;

  Clock clock;
  double seq_time = 0, simd_time = 0;
  int kase = 100; while (kase--) {
    clock.tic();
    auto sequential = BandedSmithWaterman(3, -3, -2)(fasta[0], fasta[1]);
    seq_time += clock.toc();

    clock.tic();
    auto simd = StripedSmithWaterman(3, -3, -2)(fasta[0], fasta[1]);
    simd_time += clock.toc();

    if (kase == 0) {
      std::cout << "seq : " << seq_time << '\n';
      std::cout << sequential << '\n';
      std::cout << "simd: " << simd_time << '\n';
      std::cout << simd << '\n'; 
      std::cout << "speed up: " << seq_time / simd_time << "x\n";
    }
  }
}
