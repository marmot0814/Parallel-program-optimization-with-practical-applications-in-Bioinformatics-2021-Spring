#include <cassert>

#include "utils/clock.h"

Clock::Clock() : is_start(false) { }

void Clock::tic() {
  assert(not is_start);
  prev = std::chrono::high_resolution_clock::now();
  is_start = true;
}

double Clock::toc() {
  assert(is_start);
  is_start = false;
  auto curr = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::duration<double>>(curr - prev).count();
}
