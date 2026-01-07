/*****************************************************************************
 * @file LatestKBuffer .hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 滚动更新缓冲区
 * 1、基于环形缓冲区思路设计，仅保留最新的k个数据
 * 2、仅允许追加写入新x数据和定点修改数据，禁止移除数据
 * 3、支持状态管理shutdown()和restart()，支持DEBUG中使用'<<'输出该容器的状态
 * 4、目前设计为禁止拷贝，移动容器
 * @version 0.1
 * @date 2025-10-25
 *
 * @copyright Copyright (c) 2025
 *
 *****************************************************************************/
#ifndef _LATEST_K_BUFFER_HPP_
#define _LATEST_K_BUFFER_HPP_

#include <memory>
// debug相关
#include <cassert>
#include <iostream>
#include <cstring>

namespace track_project::trackmanager
{

    template <typename T>
    class LatestKBuffer
    {

    public:
        // 构造函数
        explicit LatestKBuffer(size_t capacity)
            : capacity_(capacity), head_(0), tail_(0), size_(0), full_(false)
        {
            assert(capacity_ > 0 && "申请了过小的内存！");
            buffer_ = std::make_unique<T[]>(capacity_); // 在构造函数体内初始化
        }

        // 禁止构造、赋值、拷贝，LatestKBuffer 是唯一的，拥有单独的ID号码
        LatestKBuffer(const LatestKBuffer &) = delete;
        LatestKBuffer &operator=(const LatestKBuffer &) = delete;

        // 允许移动
        LatestKBuffer(LatestKBuffer &&other) = default;
        LatestKBuffer &operator=(LatestKBuffer &&other) = default;

        ~LatestKBuffer() = default;

        /*****************************************************************************
         * @brief 数据清空
         *****************************************************************************/
        void clear() noexcept
        {
            head_ = 0;
            tail_ = 0;
            size_ = 0;
            full_ = false;
        }

        /*****************************************************************************
         * @brief 数据追加写入部分实现
         * 实现了追加写入功能，确保数据一定在结尾追加
         *****************************************************************************/
        // 1. 数据写入
        void push(const T &item) noexcept
        {
            buffer_[head_] = item;
            _advance_head();
        }

        // 2. 移动写入
        void push(T &&item) noexcept
        {
            buffer_[head_] = std::move(item);
            _advance_head();
        }

        // 3. 原地构造
        template <typename... Args>
        void emplace(Args &&...args) noexcept
        {
            buffer_[head_] = T(std::forward<Args>(args)...);
            _advance_head();
        }

        /*****************************************************************************
         * @brief 基础参数访问/修改功能
         *****************************************************************************/
        // 定点修改数据
        T &operator[](size_t index) noexcept
        {
            return buffer_[(tail_ + index) % capacity_];
        }

        // 基础参数访问
        const T &operator[](size_t index) const noexcept
        {
            return buffer_[(tail_ + index) % capacity_];
        }

        // 批量拷贝到目标数组
        size_t copy_to(T *dest, size_t maxCount) const noexcept
        {
            if (!dest || maxCount == 0 || empty())
                return 0;

            size_t actualCount = std::min(size_, maxCount);

            if (tail_ + actualCount <= capacity_)
            {
                // 连续内存，直接拷贝
                _memcpy(dest, &buffer_[tail_], actualCount);
            }
            else
            {
                // 需要分两段拷贝
                size_t firstChunk = capacity_ - tail_;
                _memcpy(dest, &buffer_[tail_], firstChunk);
                _memcpy(dest + firstChunk, &buffer_[0], (actualCount - firstChunk));
            }

            return actualCount;
        }

        // 基本信息
        size_t capacity() const noexcept { return capacity_; }
        size_t size() const noexcept { return size_; }
        bool empty() const noexcept { return size_ == 0; }
        bool full() const noexcept { return full_; }

        /*****************************************************************************
         * @brief DEBUG用
         *****************************************************************************/
        friend std::ostream &operator<<(std::ostream &os, const LatestKBuffer &buffer)
        {
            os << "LatestKBuffer [";
            os << "capacity=" << buffer.capacity_ << ", ";
            os << "size=" << buffer.size_ << ", ";
            os << "head=" << buffer.head_ << ", ";
            os << "tail=" << buffer.tail_ << ", ";
            os << "full=" << std::boolalpha << buffer.full_ << ", ";
            os << "]";

            // 有数据时输出数据
            os << " data=[";
            for (size_t i = 0; i < buffer.size_; ++i)
            {
                os << buffer[i];
                if (i < buffer.size_ - 1)
                    os << ", ";
            }
            os << "]";

            return os;
        }

    private:
        std::unique_ptr<T[]> buffer_; // data区域，内存格式为连续的数组
        size_t capacity_;             // 作为计数器不对外输出结果，不适合使用u32或u64
        size_t head_;
        size_t tail_;
        size_t size_;
        bool full_;

        // 智能拷贝控制
        void _memcpy(T *dest, const T *src, size_t count) const noexcept
        {
            if constexpr (std::is_trivially_copyable_v<T>)
            {
                // POD类型：用memcpy获得最高性能
                std::memcpy(dest, src, count * sizeof(T));
            }
            else
            {
                // 非POD类型：用循环保证安全
                for (size_t i = 0; i < count; ++i)
                {
                    dest[i] = src[i];
                }
            }
        }

        // 数据写入位置控制
        void _advance_head() noexcept
        {
            if (full_)
            {
                tail_ = (tail_ + 1) % capacity_;
            }
            else
            {
                ++size_; // 只有非满时才增加size
            }

            head_ = (head_ + 1) % capacity_;
            full_ = (head_ == tail_);
        }
    };
} // namespace track_project::trackmanager
#endif // _UNIQUE_APPEND_RING_BUFFER_H_
