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
  ) : std::string(seq), name(name), bg(bg), ed(bg + seq.size()) {}

  Sequence(
    const std::string &name,
    const std::string &seq
  ) : Sequence(name, seq, 0) {}

  std::string verbose(const std::string &gene) const {
    std::stringstream ss;
    ss << name << ": ";
    ss << std::setw(8) << std::right << bg << " ";
    ss << gene << " " << ed - 1;
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

  void add(char c) {
    (*this) += c;
    ed++;
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

  void addM(char c1, char c2) {
    assert(c1 == c2);
    seq1.add(c1);
    seq2.add(c2);
    alignment += "M";
  }
  void addD(char c1) {
    seq1.add(c1);
    alignment += "D";
  }
  void addI(char c2) {
    seq2.add(c2);
    alignment += "I";
  }
  void addR(char c1, char c2) {
    assert(c1 != c2);
    seq1.add(c1);
    seq2.add(c2);
    alignment += "R";
  }
};

class SmithWaterman {
  int M, R, I, D;
  std::vector<std::vector<int>> dp;
 public:
  SmithWaterman(int M, int R, int I, int D) : M(M), R(R), I(I), D(D) {
    assert(M > 0);
    assert(R < 0);
    assert(I < 0);
    assert(D < 0);
  }
  BLAST operator() (const Sequence &seq1, const Sequence &seq2) {
    int len1 = seq1.size(), len2 = seq2.size();
    dp = std::vector<std::vector<int>>(len1 + 1, std::vector<int>(len2 + 1, 0));
    for (int i = 0; i < len1; i++) {
      for (int j = 0; j < len2; j++) {
        if (seq1[i] == seq2[j])
          dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j] + M);
        dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i + 1][j] + I);
        dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j + 1] + D);
        dp[i + 1][j + 1] = std::max(dp[i + 1][j + 1], dp[i][j] + R); 
      }
    }

    int maxi = 0, maxj = 0;
    for (int i = 0; i <= len1; i++)
      for (int j = 0; j <= len2; j++)
        if (dp[maxi][maxj] < dp[i][j])
          maxi = i, maxj = j;
    return answer(maxi, maxj, seq1, seq2);
  }
  BLAST answer(int i, int j, const Sequence &seq1, const Sequence &seq2) {
    int len1 = dp.size() - 1, len2 = dp[0].size() - 1;
    assert(0 <= i && i <= len1);
    assert(0 <= j && j <= len2);
    if (dp[i][j] == 0)
      return BLAST(Sequence(seq1.name, "", i), Sequence(seq2.name, "", j), "");
    if (dp[i][j] == dp[i - 1][j - 1] + M && seq1[i - 1] == seq2[j - 1]) {
      auto ret = answer(i - 1, j - 1, seq1, seq2);
      ret.addM(seq1[i - 1], seq2[j - 1]);
      return ret;
    }
    if (dp[i][j] == dp[i - 1][j - 1] + R) {
      auto ret = answer(i - 1, j - 1, seq1, seq2);
      ret.addR(seq1[i - 1], seq2[j - 1]);
      return ret;
    }
    if (dp[i][j] == dp[i][j - 1] + I) {
      auto ret = answer(i, j - 1, seq1, seq2);
      ret.addI(seq2[j - 1]);
      return ret;
    }
    if (dp[i][j] == dp[i - 1][j] + D) {
      auto ret = answer(i - 1, j, seq1, seq2);
      ret.addD(seq1[i - 1]);
      return ret;
    }
    assert(false);
  }
};

int main() {

  std::cout << SmithWaterman(1, -3, -10, -2)(
    Sequence("Seq1", "CCAATGCCACAAAACATCTGTCTCTAACTGGTGTGTGTGT", 0),
    Sequence("Seq2", "CCAGCCCAAAATCTGTTTTAATGGTGGATTTGTGT", 0)
  ) << '\n';
}
