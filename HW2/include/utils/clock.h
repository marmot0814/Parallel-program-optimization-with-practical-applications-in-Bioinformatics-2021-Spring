#pragma once

#include <chrono>

class Clock {
  std::chrono::high_resolution_clock::time_point prev;
  bool is_start;
 public:
  Clock();
  void tic();
  double toc();
};
