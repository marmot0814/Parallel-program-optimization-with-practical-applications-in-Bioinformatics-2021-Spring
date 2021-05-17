#include <simdpp/simd.h>

#include "types/table/striped_table.h"

int StripedTable::striped(int i) const {
  if (i < 0) return i;
  const int p = SIMDPP_FAST_INT32_SIZE;
  const int t = m / p;
  return i % t * p + i / t;
}

StripedTable::StripedTable(int n, int m)
  : Table(n, (m + SIMDPP_FAST_INT32_SIZE - 1) / SIMDPP_FAST_INT32_SIZE * SIMDPP_FAST_INT32_SIZE) {}

int StripedTable::get(int i, int j) const {
  return Table::get(i, striped(j));
}
int* StripedTable::ptr(int i, int j) {
  return Table::ptr(i, j);
}
