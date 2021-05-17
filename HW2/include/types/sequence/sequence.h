#pragma once

#include <string>

struct Sequence : public std::string {
  std::string name;
  int bg, ed;
  Sequence(
    const std::string &,
    const std::string &,
    int
  );

  Sequence(
    const std::string &,
    const std::string &
  );

  std::string verbose(
    const std::string &
  ) const;

  std::string align(
    const std::string &,
    const char gap_ch
  ) const;

  void add(char);
  char get(int i) const;
};

