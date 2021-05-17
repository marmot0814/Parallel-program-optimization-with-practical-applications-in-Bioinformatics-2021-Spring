#include <cassert>

#include "types/table/table.h"

Table::Table(
  int n, int m
) : std::vector<std::vector<int>>(n, std::vector<int>(m, 0)), n(n), m(m) {
  assert(n > 0 && m > 0);
}
int Table::get(int i, int j) const {
  if (i < 0 || j < 0) return 0;
  return (*this)[i][j];
}
int* Table::ptr(int i, int j) {
  assert(0 <= i && i < n);
  assert(0 <= j && j < m);
  return (*this)[i].data() + j;
};
