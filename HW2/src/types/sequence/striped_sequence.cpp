#include <simdpp/simd.h>

#include "types/sequence/striped_sequence.h"

int StripedSequence::striped(int i) const {
  if (i < 0) return i;
  const int p = SIMDPP_FAST_INT32_SIZE;
  const int t = size() / p;
  return i % p * t + i / p;
}

StripedSequence::StripedSequence(const Sequence &seq) : Sequence(seq) {}
 
char StripedSequence::striped_get(int i) const {
  return Sequence::get(striped(i));
}

int StripedSequence::size() const {
  const int p = SIMDPP_FAST_INT32_SIZE;
  return (Sequence::size() + p - 1) / p * p;
}
