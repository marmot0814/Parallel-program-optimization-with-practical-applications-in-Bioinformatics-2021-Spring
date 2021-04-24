#include <vector>
#include <cassert>
#include <string>
#include <iomanip>
#include <string>
#include <iostream>
#include <sstream>

struct Sequence : public std::string {
  std::string name;
  int bg, ed;
  Sequence(
    const std::string &name,
    const std::string &seq,
    int bg
  ) : std::string(seq), name(name), bg(bg), ed(bg + seq.size() - 1) {
    assert((int)seq.size() > 0);
  }

  Sequence(
    const std::string &name,
    const std::string &seq
  ) : Sequence(name, seq, 0) {}

  std::string verbose(const std::string &gene) const {
    std::stringstream ss;
    ss << name << ": ";
    ss << std::setw(8) << std::right << bg << " ";
    ss << gene << " " << ed;
    return ss.str();
  }
  friend std::ostream& operator<< (std::ostream &os, const Sequence &seq) {
    return os << seq.verbose(seq), os;
  }

  std::string align(const std::string &alignment, const char gap_ch) const {
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
};

class FASTA : public std::vector<Sequence> {
 public:
};

class BLAST {
  Sequence seq1, seq2;
  std::string alignment;
 public:
  BLAST(
    const Sequence &seq1,
    const Sequence &seq2,
    const std::string &alignment
  ) : seq1(seq1), seq2(seq2), alignment(alignment) {

  }
  friend std::ostream& operator<< (std::ostream &os, const BLAST &blast) {
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
};

int main() {
  std::cout << BLAST(
    Sequence("Seq1", "CCAATGCCACAAAACATCTGTCTCTAACTGGTGTGTGTGT", 453),
    Sequence("Seq2", "CCAGCCCAAAATCTGTTTTAATGGTGGATTTGTGT", 17),
    "MMMDDMMMDMMMMDDMMMMMMDMRMMMDMMMMMIIMRMMMMM"
  ) << '\n';
}
