#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace multicore
{

using IndexedTask = std::function<void(std::size_t task_index, std::size_t worker_index)>;

[[nodiscard]] std::size_t default_worker_count() noexcept;

class TaskExecutor
{
public:
    // A value of zero selects std::thread::hardware_concurrency().
    explicit TaskExecutor(std::size_t worker_count = 0) noexcept;

    [[nodiscard]] std::size_t worker_count() const noexcept;

    // Each index in [0, task_count) is executed exactly once unless a task
    // throws. The first exception is rethrown after all active workers stop.
    // A worker processes only one task at a time, so worker_index can be used
    // to address caller-owned per-worker state without a lock.
    void for_each(std::size_t task_count, const IndexedTask& task) const;

    template<typename Function>
    [[nodiscard]] auto map(std::size_t task_count, Function&& function) const
    {
        using Result = std::remove_cvref_t<std::invoke_result_t<
            Function&,
            std::size_t,
            std::size_t
        >>;
        static_assert(!std::is_void_v<Result>, "TaskExecutor::map requires a result");

        std::vector<std::optional<Result>> slots(task_count);
        for_each(
            task_count,
            [&](std::size_t task_index, std::size_t worker_index) {
                slots[task_index].emplace(function(task_index, worker_index));
            }
        );

        std::vector<Result> results;
        results.reserve(task_count);
        for(auto& slot : slots)
            results.push_back(std::move(slot).value());
        return results;
    }

private:
    std::size_t worker_count_;
};

} // namespace multicore
