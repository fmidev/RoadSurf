// Mostly based on https://blind.guru/simple_cxx11_workqueue.html

#pragma once

#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

//namespace Queue
//{
template <typename Type, typename Queue = std::queue<Type>>
class WorkQueue : std::mutex, std::condition_variable
{
 public:
  template <typename Function>
  WorkQueue(Function &&function, unsigned int concurrency = std::thread::hardware_concurrency())
  {
    try
    {
      if (not concurrency)
        throw std::runtime_error("Concurrency must not be zero");

      for (unsigned int count{}; count < concurrency; ++count)
        threads.emplace_back(static_cast<void (WorkQueue::*)(Function)>(&WorkQueue::consume),
                             this,
                             std::forward<Function>(function));
    }
    catch (...)
    {
      throw std::runtime_error("Operation failed!");
    }
  }

  // disable move
  WorkQueue(WorkQueue &&) = delete;
  WorkQueue &operator=(WorkQueue &&) = delete;

  template <typename... Args>
  WorkQueue &operator()(Args &&...args)
  {
    try
    {
      unique_lock lock{*this};
      while (queue.size() == threads.size())
        wait(lock);
      queue.emplace(std::forward<Args>(args)...);
      notify_one();
      return *this;
    }
    catch (...)
    {
      throw std::runtime_error("Operation failed!");
    }
  }

  void join_all()
  {
    try
    {
      lock();
      if (!done)
      {
        done = true;
        notify_all();
        unlock();
        threads.join();
      }
    }
    catch (...)
    {
      std::runtime_error("Operation failed!");
    }
  }

  ~WorkQueue() { join_all(); }

 private:
  bool done = false;
  Queue queue;
  struct : std::vector<std::thread>
  {
    void join() { for (auto& t : *this) { t.join(); } }
  } threads;

  using lock_guard = std::lock_guard<std::mutex>;
  using unique_lock = std::unique_lock<std::mutex>;

  template <typename Function>
  void consume(Function process)
  {
    try
    {
      unique_lock lock{*this};
      while (true)
      {
        if (not queue.empty())
        {
#if 0
      // doesn't work in g++ 4.8.5
          Type item{std::move(queue.front())};
#else
          Type item;
          std::swap(item, queue.front());
#endif
          queue.pop();
          notify_one();
          lock.unlock();
          process(item);
          lock.lock();
        }
        else if (done)
        {
          break;
        }
        else
        {
          wait(lock);
        }
      }
    }
    catch (...)
    {
      std::runtime_error("Operation failed!");
    }
  }
};

//} 
