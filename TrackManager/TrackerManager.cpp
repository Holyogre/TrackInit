#include "TrackerManager.hpp"
#include "../utils/Logger.hpp"

namespace track_project::trackmanager
{

    // 定义最大外推次数,当>MAX_EXTRAPOLATION_TIMES时，终结对应航迹
    constexpr std::uint32_t MAX_EXTRAPOLATION_TIMES = 3;

    // 构造函数：预开辟空间，空间上构造目标
    TrackerManager::TrackerManager(std::uint32_t track_size, std::uint32_t track_length) : next_track_id_(1), track_length(track_length)
    {

        // 预分配内存，提高性能
        buffer_pool_.reserve(track_size);
        track_id_to_pool_index_.reserve(track_size); // rehash会给creat任务带来约一倍延时
        free_slots_.reserve(track_size);

        // 初始化内存池
        for (std::uint32_t i = 0; i < track_size; ++i)
        {
            buffer_pool_.push_back(TrackerContainer(track_length));
            free_slots_.push_back(i);
        }
    }

    // 申请新航迹存储器，申请一个最新的空闲内存池，更新ID编号
    std::uint32_t TrackerManager::create_track()
    {

        if (free_slots_.empty())
        {
            LOG_DEBUG << "无法申请新航迹，已不具备新航迹";
            return 0; // 内存池已满
        }

        std::uint32_t pool_index = free_slots_.back();
        free_slots_.pop_back();

        // 修改索引
        std::uint32_t track_id = next_track_id_;
        track_id_to_pool_index_[track_id] = pool_index;

        // 修改container属性
        buffer_pool_[pool_index].header.start(next_track_id_);

        // 修改计数器
        next_track_id_++;

        return track_id;
    }

    // 释放航迹存储器，释放航迹-内存池队，添加空空闲内存池编号到末尾，调用结构体内置clear
    bool TrackerManager::delete_track(std::uint32_t track_id)
    {

        auto it = track_id_to_pool_index_.find(track_id);

        // 异常处理
        if (it == track_id_to_pool_index_.end())
        {
            LOG_DEBUG << "删除航迹失败，该航迹号" << track_id << "不存在";
            return false; // 航迹不存在
        }

        std::uint32_t pool_index = it->second;

        // 清空对应的缓冲区
        buffer_pool_[pool_index].clear();

        // 释放资源，放到空内存区中
        track_id_to_pool_index_.erase(it);
        free_slots_.push_back(pool_index);

        return true;
    }

    // 存放一个数据点，依据是否外推修改航迹状态，并管理航迹，FALSE时候请求所有流水线删除id对应的容器
    bool TrackerManager::push_track_point(std::uint32_t track_id, TrackPoint point)
    {
        // 搜索目标航迹
        auto it = track_id_to_pool_index_.find(track_id);

        // 异常处理
        if (it == track_id_to_pool_index_.end())
        {
            LOG_DEBUG << "添加航迹点失败，该航迹号" << track_id << "不存在";
            return false; // 航迹不存在
        }

        // 获取航迹
        std::uint32_t pool_index = it->second;
        TrackerContainer &track = buffer_pool_[pool_index];

        // 存入数据
        track.data.push(point);

        // 若航迹外推次数过多或是置信度过低，请求删除航迹
        if (track.header.state == 2)
        {
            TrackerManager::delete_track(track.header.track_id);
            return false;
        }

        // 数据处理
        track.header.point_num = static_cast<std::uint32_t>(track.data.size()); // 更新点迹数量
        if (point.is_associated)                                                // 关联点继续
        {
            if (track.header.extrapolation_count > 0)
            {
                track.header.extrapolation_count--;
            }
            track.header.state = 0;
        }
        else if (track.header.extrapolation_count < MAX_EXTRAPOLATION_TIMES) // 未超过最大关联次数
        {
            track.header.extrapolation_count++;
            track.header.state = 1;
        }
        else // 恰好超过最大外推次数
        {
            track.header.state = 2;
        }

        return true;
    }

    // 将source源航迹数据追加到target目标航迹，然后交换ID，最后删除源航迹对应的内存池
    bool TrackerManager::merge_tracks(std::uint32_t source_track_id, std::uint32_t target_track_id)
    {
        // 搜索目标航迹
        auto target_it = track_id_to_pool_index_.find(target_track_id);
        auto source_it = track_id_to_pool_index_.find(source_track_id);

        // 异常处理
        if (target_it == track_id_to_pool_index_.end() || source_it == track_id_to_pool_index_.end())
        {
            LOG_DEBUG << "航迹融合失败，源航迹" << source_track_id << "不存在";
            return false; // 航迹不存在
        }

        // 获取航迹
        std::uint32_t target_pool_index = target_it->second;
        TrackerContainer &target_track = buffer_pool_[target_pool_index];
        std::uint32_t source_pool_index = source_it->second;
        TrackerContainer &source_track = buffer_pool_[source_pool_index];

        // 异常处理
        std::uint32_t target_size = static_cast<std::uint32_t>(target_track.data.size());
        std::uint32_t source_size = static_cast<std::uint32_t>(source_track.data.size());
        if (target_size < MAX_EXTRAPOLATION_TIMES || source_size < MAX_EXTRAPOLATION_TIMES)
        {
            LOG_DEBUG << "航迹融合失败，目标航迹" << source_track_id << "不存在";
            return false;
        }

        // 1.更新外推点为新航迹点
        for (size_t i = 1; i <= MAX_EXTRAPOLATION_TIMES; i++)
        {
            target_track.data[target_size - i] = source_track.data[source_size - i];
        }

        // 2.交换航迹状态，更新内存池索引
        std::swap(target_track.header, source_track.header);
        source_it->second = target_pool_index;

        // 3.删除target_id对应的容器
        delete_track(target_track_id);

        return true;
    }

    // 重置整个缓冲区,所有内存池改为空弦状态，重置内存编号
    void TrackerManager::clear_all()
    {

        for (auto &buffer : buffer_pool_)
        {
            buffer.clear();
        }
        track_id_to_pool_index_.clear();
        free_slots_.clear();

        // 重建空闲槽位列表
        for (std::uint32_t i = 0; i < buffer_pool_.size(); ++i)
        {
            free_slots_.push_back(i);
        }

        next_track_id_ = 1;
    }

    // 获取活跃的航迹号,返回一个包含所有活跃航迹ID的向量
    std::vector<std::uint32_t> TrackerManager::get_active_track_ids() const
    {
        std::vector<std::uint32_t> active_ids;
        active_ids.reserve(track_id_to_pool_index_.size());
        for (const auto &pair : track_id_to_pool_index_)
        {
            active_ids.push_back(pair.first);
        }
        return active_ids;
    }
    // 获取id对应的航迹头部只读引用，若不存在返回nullptr
    const TrackerManager::TrackerHeader *TrackerManager::get_header_ref(std::uint32_t track_id) const
    {
        auto it = track_id_to_pool_index_.find(track_id);
        if (it == track_id_to_pool_index_.end())
            return nullptr;
        const TrackerContainer &track = buffer_pool_[it->second];
        return &track.header;
    }

    // 获取id对应的航迹数据只读引用，若不存在返回nullptr
    const LatestKBuffer<TrackerManager::TrackPoint> *TrackerManager::get_data_ref(std::uint32_t track_id) const
    {
        auto it = track_id_to_pool_index_.find(track_id);
        if (it == track_id_to_pool_index_.end())
            return nullptr;
        const TrackerContainer &track = buffer_pool_[it->second];
        return &track.data;
    }
}
