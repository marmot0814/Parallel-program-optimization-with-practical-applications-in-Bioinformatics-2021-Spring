#pragma once

#include "types/table/table.h"

class StripedTable : public Table {
  int striped(int) const;

 public:
  StripedTable(int, int);
  int get(int, int) const;
  int* ptr(int, int);
};

