#include "LogicBased.hpp"
#include "../utils/Logger.hpp"

namespace track_project::trackinit
{
    ProcessStatus LogicBasedTracker::process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_tracks)
    {
        if (points.empty())
        {
            return ProcessStatus::NO_POINT; // 没有点迹，这批次数据不处理，不影响后续关联
        }

        // 移动数据
        shift_batches_and_hypotheses();
        point_batches_[0] = points;            // 存储最新批次数据
        batch_timestamps_[0] = points[0].time; // 存储最新批次时间戳

        // 清空输出航迹
        new_tracks.clear();

        // ==================== STEP1: 扩展假设树 ====================
        for (auto point : points)
        {
            // 查询满足条件的假设节点
            auto candidate_nodes = query_nodes_by_location(point.longitude, point.latitude);

            if (candidate_nodes.empty())
            {
                continue; // 没有满足条件的假设节点，继续处理下一个点迹
            }

            // 没有满足条件的假设节点，以该点迹为基础生成新的假设树
            ProcessStatus status = extend_hypotheses(point, candidate_nodes);
        }
    }

    // 确保[0]索引对应的总是最新的数据
    void LogicBasedTracker::shift_batches_and_hypotheses()
    {
        // 1. 清空最旧批次的数据和假设节点
        point_batches_[3].clear();
        hypothesis_layers_[3].clear();

        // 2. 将批次数据整体后移（索引2→3, 1→2, 0→1）
        for (int i = 3; i > 0; --i)
        {
            // 移动点迹批次
            point_batches_[i] = std::move(point_batches_[i - 1]);
            // 移动假设节点层
            hypothesis_layers_[i] = std::move(hypothesis_layers_[i - 1]);
            // 移动时间戳
            batch_timestamps_[i] = batch_timestamps_[i - 1];
        }

        // 4. 清空最新的批次位置（索引0）
        point_batches_[0].clear();
        hypothesis_layers_[0].clear();
    }

    ProcessStatus LogicBasedTracker::extend_hypotheses(const TrackPoint &point, const std::vector<HypothesisNode *> &parent_nodes)
    {
        // TODO
    }
}