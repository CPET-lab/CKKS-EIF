#include "multicore.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace multicore
{

std::size_t default_worker_count() noexcept
{
    const unsigned int detected = std::thread::hardware_concurrency();
    return detected == 0 ? 1 : static_cast<std::size_t>(detected);
}

TaskExecutor::TaskExecutor(std::size_t worker_count) noexcept
    : worker_count_(worker_count == 0 ? default_worker_count() : worker_count)
{}

std::size_t TaskExecutor::worker_count() const noexcept
{
    return worker_count_;
}

void TaskExecutor::for_each(
    std::size_t task_count,
    const IndexedTask& task) const
{
    if(task_count == 0)
        return;
    if(!task)
        throw std::invalid_argument("TaskExecutor::for_each requires a task");

    const std::size_t active_worker_count = std::min(worker_count_, task_count);
    std::atomic<std::size_t> next_task = 0;
    std::atomic<bool> stop = false;
    std::exception_ptr first_exception;
    std::mutex exception_mutex;

    auto run_worker = [&](std::size_t worker_index) {
        while(!stop.load(std::memory_order_acquire))
        {
            const std::size_t task_index = next_task.fetch_add(
                1,
                std::memory_order_relaxed
            );
            if(task_index >= task_count)
                return;

            try
            {
                task(task_index, worker_index);
            }
            catch(...)
            {
                {
                    std::lock_guard<std::mutex> lock(exception_mutex);
                    if(!first_exception)
                        first_exception = std::current_exception();
                }
                stop.store(true, std::memory_order_release);
                return;
            }
        }
    };

    std::vector<std::jthread> workers;
    workers.reserve(active_worker_count);
    for(std::size_t worker_index = 0;
        worker_index < active_worker_count;
        ++worker_index)
    {
        workers.emplace_back(run_worker, worker_index);
    }

    for(auto& worker : workers)
        worker.join();

    if(first_exception)
        std::rethrow_exception(first_exception);
}

} // namespace multicore
