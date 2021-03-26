#include <iostream>
#include "thread_pool.hpp"

#include <atomic>
#include <mutex>
#include <vector>

#include <random>
#include <cstdlib>
#include <ctime>

std::atomic<int> print_1_cnt = 0;
std::mt19937 gen;
std::uniform_int_distribution<int> dist(0, 9);
class print_1 {
 public:
  void operator()() {
    tscout << std::to_string(dist(gen) % 2) + "\n";
    print_1_cnt--;
  }
};

void print_2() {
  while (print_1_cnt);
  tscout << "2\n";
}

int main() {
  srand(time(0));
  ThreadPool thread_pool(5);

  std::vector<std::future<void>> res;
  for (int i = 0; i < 8; i++) {
    res.push_back(thread_pool.send(print_1()));
    print_1_cnt++;
  }

  for (int i = 0; i < 4; i++) {
    res.push_back(thread_pool.send(print_2));
  }

  for (auto &r : res)
    r.get();
}
