#pragma once
#include <iostream>

#include "types/sequence/sequence.h"

class BLAST {
  Sequence seq1, seq2;
  std::string alignment;
 public:
  BLAST(
    const Sequence &,
    const Sequence &,
    const std::string &
  );

  friend std::ostream& operator<< (
    std::ostream &,
    const BLAST &
  );

  void addM(char, char);
  void addD(char);
  void addI(char);
  void addR(char, char);
};
