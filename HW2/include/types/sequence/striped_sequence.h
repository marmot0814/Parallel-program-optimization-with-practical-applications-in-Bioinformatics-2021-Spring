#pragma once

#include "types/sequence/sequence.h"

struct StripedSequence : public Sequence {
  int striped(int) const;

 public:
  StripedSequence(const Sequence &);
  char striped_get(int) const;
  int size() const;
};
