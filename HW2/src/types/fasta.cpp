#include <cassert>

#include "types/fasta.h"

std::istream& operator>> (std::istream &is, FASTA &fasta) {
  fasta.clear();
  std::string name, seq;
  while (is >> name) {
    assert(name[0] == '>');
    name = name.substr(1);
    is >> seq;
    fasta.push_back(Sequence(name, seq));
  }
  return is;
}
