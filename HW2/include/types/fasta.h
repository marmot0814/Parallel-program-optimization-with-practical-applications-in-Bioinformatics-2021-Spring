#pragma once

#include <vector>
#include <iostream>

#include "types/sequence.h"

class FASTA : public std::vector<Sequence> {
 public:
  friend std::istream& operator>> (std::istream &, FASTA &);
};
