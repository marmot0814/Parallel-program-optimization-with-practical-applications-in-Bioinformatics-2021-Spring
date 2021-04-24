#pragma once

#include <thread>
#include <vector>
#include <future>
#include <atomic>
#include <queue>
#include <functional>
#include <sstream>
#include <condition_variable>

class ThreadSafeCout {
  std::mutex cout_mutex;
 public:
  template <typename T>
  ThreadSafeCout& operator<< (const T &t) {
    std::unique_lock<std::mutex> lock(cout_mutex);
    return std::cout << t, *this;
  }
} tscout;

class ThreadPool {
  bool stop;

  std::vector<std::thread> workers;
  std::condition_variable workers_cv;
  std::mutex jobs_mutex;

  std::queue<std::function<void()>> jobs;

  std::string exit_logger(auto id) {
    std::stringstream ss;
    ss << "[Thread 0x" << std::hex << std::this_thread::get_id() << " exited]\n";
    return ss.str();
  }

 public:
  ThreadPool(int num_of_workers) : stop(false) {
    while (num_of_workers--)
      workers.emplace_back([&](){
        while (1) {

          std::unique_lock<std::mutex> lock(jobs_mutex);
          workers_cv.wait(lock, [&](){
            return !jobs.empty() || stop;
          });

          if (jobs.empty()) {
            tscout << exit_logger(std::this_thread::get_id());
            return ;
          }

          auto job = std::move(jobs.front());
          jobs.pop();

          lock.unlock();
          job();
        }
      });
  }
  ~ThreadPool() {
    stop = true;
    workers_cv.notify_all();
    for (auto &worker : workers)
      worker.join();
  }

  template <typename F, typename... Args>
  std::future<std::invoke_result_t<F, Args...>> send(F&& func, Args&& ...args) {

    auto job_ptr = std::make_shared<
      std::packaged_task<
        std::invoke_result_t<F, Args...>()
      >
    >(
      std::bind(
        std::forward<F>(func),
        std::forward<Args>(args)...
      )
    );

    auto res = job_ptr->get_future();

    std::unique_lock<std::mutex> lock(jobs_mutex);
    jobs.emplace([job_ptr](){ (*job_ptr)(); });
    lock.unlock();

    workers_cv.notify_one();
    return res;
  }
};
