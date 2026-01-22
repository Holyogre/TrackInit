/*****************************************************************************
 * @file ObjectPool.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 对象池ObjectPool<T> 用于取代频繁的new/delete操作，提升性能
 * 因为由扩容需求。所以引入了deque和unordered_map，比vector实现慢不少
 * 不适用于超高频率申请释放场景，不能随便移植
 *
 * 主要功能：
 * T *acquire_target()：快速分配一个对象
 * bool release_target(T * ptr)：快速释放一个对象
 * pre_allocate(size_t count)：预分配指定数量的对象槽位，可以超过初始容量
 * const std::vector<T*> get_allocated_ptr()：获取当前已分配对象的指针列表
 * reclaim_all()：清除所有申请的对象，回收内存
 *
 * 特点：
 * 1. ⚠️类型T必须有clear方法，否则报错
 * 2. 不具备自动析构功能，只是取代new/delete实现快速分配和释放
 *
 * @version 0.3
 * @date 2025-12-20
 *
 * @copyright Copyright (c) 2025
 *
 *****************************************************************************/
#ifndef _TRACK_INIT_OBJECT_POOL_HPP_
#define _TRACK_INIT_OBJECT_POOL_HPP_

#include <deque>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <stdexcept>
#include <cstdint>
namespace track_project::trackinit
{

    template <typename T>
    class ObjectPool
    {
        // 允许测试代码访问私有成员：友元模板
        template <typename U>
        friend struct ObjectPoolTestFriend;

    private:
        // 核心数据结构
        std::deque<T> memory_pool_;                  // 使用deque确保指针不失效
        std::vector<uint8_t> alloc_flags_;           // 分配标志位（0/1）
        std::vector<size_t> free_slots_;             // 空闲槽位栈，使用 vector 充当 LIFO
        std::unordered_map<T *, size_t> ptr_to_idx_; // 指针到索引的映射

        // 禁用复制
        ObjectPool(const ObjectPool &) = delete;
        ObjectPool &operator=(const ObjectPool &) = delete;

    public:
        /*****************************************************************************
         * @brief 构造对象池
         *
         * @param initial_targets 预分配空间大小（强制）
         *****************************************************************************/
        explicit ObjectPool(size_t initial_targets)
        {
            if (initial_targets == 0)
            {
                initial_targets = 1; // 最小分配1个槽位
            }
            pre_allocate(initial_targets);
        }

        ~ObjectPool() = default;

        /*****************************************************************************
         * @brief 申请一个对象
         * @return 指向新分配对象的指针，扩容失败抛std::bad_alloc
         *****************************************************************************/
        T *acquire_target()
        {
            if (free_slots_.empty())
            {
                // 两倍扩容策略，初始大于1由构造函数保证
                size_t new_size = memory_pool_.size() * 2;
                expand_pool(new_size);
            }

            size_t idx = free_slots_.back();
            free_slots_.pop_back();

            alloc_flags_[idx] = 1;
            T *ptr = &memory_pool_[idx];
            ptr_to_idx_[ptr] = idx;

            ptr->clear(); // 确保对象初始状态
            return ptr;
        }

        /*****************************************************************************
         * @brief 释放一个对象
         *
         * @param target 要释放的对象指针
         * @return true 成功释放
         * @return false 指针无效或不属于本池
         *****************************************************************************/
        bool release_target(T *target)
        {
            auto it = ptr_to_idx_.find(target);
            if (it == ptr_to_idx_.end())
            {
                return false;
            }

            size_t idx = it->second;

            // 安全检查
            if (idx >= alloc_flags_.size() || alloc_flags_[idx] == 0)
            {
                ptr_to_idx_.erase(it);
                return false;
            }

            target->clear();
            alloc_flags_[idx] = 0;
            free_slots_.push_back(idx);
            ptr_to_idx_.erase(it);

            return true;
        }

        /*****************************************************************************
         * @brief 清除所有对象的分配状态以及内部数据
         *****************************************************************************/
        void clear_all()
        {
            free_slots_.clear();
            free_slots_.reserve(memory_pool_.size());
            for (size_t i = 0; i < memory_pool_.size(); ++i)
            {
                if (alloc_flags_[i])
                {
                    memory_pool_[i].clear();
                    alloc_flags_[i] = 0;
                }
                free_slots_.push_back(i);
            }
            ptr_to_idx_.clear();
        }

        /*****************************************************************************
         * @brief 预分配指定数量的槽位
         *
         * @param target_count 预分配的目标数量
         *****************************************************************************/
        void pre_allocate(size_t target_count)
        {
            if (target_count > memory_pool_.size())
            {
                expand_pool(target_count);
            }
        }

        /*****************************************************************************
         * @brief 获取所有已分配对象的指针数组（只读）
         * @return std::vector<T*> 指针数组，按分配顺序排列
         *****************************************************************************/
        std::vector<T *> get_allocated_ptrs() const
        {
            std::vector<T *> result;
            result.reserve(ptr_to_idx_.size());

            // 按索引顺序收集指针，保证输出有序性
            for (size_t i = 0; i < alloc_flags_.size(); ++i)
            {
                if (alloc_flags_[i] == 1)
                {
                    // 注意：这里需要const_cast，因为deque元素在逻辑上可修改
                    result.push_back(const_cast<T *>(&memory_pool_[i]));
                }
            }

            return result;
        }

        // 仅用于DEBUG的函数
        size_t get_allocated_count() const
        {
            size_t count = 0;
            // 遍历对象池或索引，统计有效指针
            for (const auto &flag : alloc_flags_)
            {
                if (flag == 1)
                {
                    ++count;
                }
            }
            return count;
        }

    private:
        // 内部扩容方法
        void expand_pool(size_t new_size)
        {
            size_t old_size = memory_pool_.size();
            if (new_size <= old_size)
                return;

            memory_pool_.resize(new_size);
            alloc_flags_.resize(new_size, 0);
            free_slots_.reserve(new_size);
            ptr_to_idx_.reserve(new_size);

            // 新槽位加入空闲栈
            for (size_t i = old_size; i < new_size; ++i)
            {
                free_slots_.push_back(i);
            }
        }
    };
} // namespace track_project::trackinit

#endif // _TRACK_INIT_OBJECT_POOL_HPP_