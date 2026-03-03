#include "LogicBasedInitiator.hpp"
#include "../utils/Logger.hpp"
#include <cmath>

namespace track_project::trackinit
{
    ProcessStatus LogicBasedInitiator::process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_tracks)
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
            auto candidate_nodes_ = query_nodes_by_location(point.longitude, point.latitude);

            if (candidate_nodes_.empty())
            {
                continue; // 没有满足条件的假设节点，继续处理下一个点迹
            }

            // 没有满足条件的假设节点，以该点迹为基础生成新的假设树
            ProcessStatus status = extend_hypotheses(point, candidate_nodes_);
        }
    }

    // 使用move语义将旧数据和假设树整体后移，清空最新批次位置。。这步不会涉及到内存分配和释放
    void LogicBasedInitiator::shift_batches_and_hypotheses()
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

    // query_nodes_by_location
    std::vector<LogicBasedInitiator::HypothesisNode *> LogicBasedInitiator::query_nodes_by_location(double x, double y) const
    {
        // 计算分辨率
        double x_bin_size = (2 * LOGIC_BASED_ABS_X_MAX) / LOGIC_BASED_MAX_X_BINS; // X轴每个bin的大小
        double y_bin_size = (2 * LOGIC_BASED_ABS_Y_MAX) / LOGIC_BASED_MAX_Y_BINS; // Y轴每个bin的大小

        // 离散化输入坐标，四舍五入，计算对应的bin索引
        size_t x_index = static_cast<size_t>(std::round((x + LOGIC_BASED_ABS_X_MAX) / x_bin_size));
        size_t y_index = static_cast<size_t>(std::round((y + LOGIC_BASED_ABS_Y_MAX) / y_bin_size));

        // 先查表查看当前点的误差分布情况，根据误差分布情况计算需要查询的bin范围
        size_t x_bin_range = error_distribution_table_[x_index * LOGIC_BASED_MAX_Y_BINS + y_index].first;  // 误差分布表格中存储的x方向误差范围，单位bin数量
        size_t y_bin_range = error_distribution_table_[x_index * LOGIC_BASED_MAX_Y_BINS + y_index].second; // 误差分布表格中存储的y方向误差范围，单位bin数量

        // 依据range，计算所有假设位置
        std::vector<HypothesisNode *> candidate_nodes;
        for (int dx = -x_bin_range; dx <= x_bin_range; ++dx)
        {
            for (int dy = -y_bin_range; dy <= y_bin_range; ++dy)
            {
                // 计算需要查询的bin索引
                size_t query_x_index = x_index + dx;
                size_t query_y_index = y_index + dy;

                // 检查索引是否越界
                if (query_x_index >= LOGIC_BASED_MAX_X_BINS || query_y_index >= LOGIC_BASED_MAX_Y_BINS)
                {
                    continue; // 索引越界，跳过
                }

                // 计算bin索引
                size_t bin_index = query_x_index * LOGIC_BASED_MAX_Y_BINS + query_y_index;
                // 将该bin中的假设节点加入候选列表
                candidate_nodes.insert(candidate_nodes.end(), current_hypothesis_index_[bin_index].begin(), current_hypothesis_index_[bin_index].end());
            }
        }
        return candidate_nodes;
    }

    ProcessStatus LogicBasedInitiator::extend_hypotheses(const TrackPoint &point, const std::vector<HypothesisNode *> &parent_nodes)
    {
        // TODO
    }
}