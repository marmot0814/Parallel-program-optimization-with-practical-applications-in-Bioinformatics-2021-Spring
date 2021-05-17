#pragma once

#include <vector>

class Table : public std::vector<std::vector<int>> {
 protected:
  int n, m;
 public:
  Table(int, int);
  int get(int, int) const;
  int* ptr(int, int);
};
