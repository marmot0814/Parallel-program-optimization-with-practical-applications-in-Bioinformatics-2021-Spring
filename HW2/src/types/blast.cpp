#include <cassert>

#include "types/blast.h"

BLAST::BLAST(
  const Sequence &seq1,
  const Sequence &seq2,
  const std::string &alignment
) : seq1(seq1), seq2(seq2), alignment(alignment) {}

std::ostream& operator<< (std::ostream &os, const BLAST &blast) {
  assert(blast.seq1.align(blast.alignment, 'I').size() == blast.seq2.align(blast.alignment, 'D').size());
  os << std::string(std::max(0, (int)blast.seq2.name.size() - (int)blast.seq1.name.size()), ' ');
  os << blast.seq1.verbose(blast.seq1.align(blast.alignment, 'I')) << '\n';

  os << std::string(std::max(blast.seq1.name.size(), blast.seq2.name.size()) + 11, ' ');
  for (auto &c : blast.alignment) {
    if (c == 'I' || c == 'D')
      os << ' ';
    else if (c == 'R')
      os << '*';
    else if (c == 'M')
      os << "|";
  }
  os << '\n';

  os << std::string(std::max(0, (int)blast.seq1.name.size() - (int)blast.seq2.name.size()), ' ');
  os << blast.seq2.verbose(blast.seq2.align(blast.alignment, 'D'));
  return os;
}

void BLAST::addM(char c1, char c2) {
  assert(c1 == c2);
  seq1.add(c1);
  seq2.add(c2);
  alignment += "M";
}
void BLAST::addD(char c1) {
  seq1.add(c1);
  alignment += "D";
}
void BLAST::addI(char c2) {
  seq2.add(c2);
  alignment += "I";
}
void BLAST::addR(char c1, char c2) {
  assert(c1 != c2);
  seq1.add(c1);
  seq2.add(c2);
  alignment += "R";
}
