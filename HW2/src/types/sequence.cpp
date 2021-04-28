#include <iomanip>
#include <cassert>

#include "types/sequence.h"

Sequence::Sequence(
  const std::string &name,
  const std::string &seq,
  int bg
) : std::string(seq), name(name), bg(bg), ed(bg + seq.size()) {}

Sequence::Sequence(
  const std::string &name,
  const std::string &seq
) : Sequence(name, seq, 0) {}

std::string Sequence::verbose(
  const std::string &gene
) const {
  std::stringstream ss;
  ss << name << ": ";
  ss << std::setw(8) << std::right << bg << " ";
  ss << gene << " " << ed - 1;
  return ss.str();
}

std::ostream& operator<< (
  std::ostream &os,
  const Sequence &seq
) {
  return os << seq.verbose(seq), os;
}

std::string Sequence::align(
  const std::string &alignment,
  const char gap_ch
) const {
  std::string ret;
  int ptr = 0;
  for (auto &c : alignment) {
    if (c == gap_ch) {
      ret += "-";
      continue;
    }
    assert(ptr < size());
    ret += (*this)[ptr++];
  }
  return ret;
}

void Sequence::add(char c) {
  (*this) += c;
  ed++;
}

char Sequence::get(int i) const {
  if (i >= size() || i < 0)
    return '$';
  return (*this)[i];
}
