#include "LogicBasedInitiator.hpp"

#include <cmath>

#include "../include/defsystem.h"
#include "../utils/Logger.hpp"

namespace track_project::trackinit
{
    LogicBasedInitiator::LogicBasedInitiator()
    {

        // 预留空间
        for (size_t i = 0; i < 4; ++i)
        {
            point_batches_[i].reserve(1000);                             // 每批次预留1000个点迹的空间
            hypothesis_layers_[i].reserve(MAX_BINS * MAX_NODE_PER_BINS); // 避免内存碎片一次开到位
        }

        // 初始化索引表，避免内存碎片一次性开到位
        current_hypothesis_index_.resize(MAX_BINS);
        history_hypothesis_index_.resize(MAX_BINS);
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            current_hypothesis_index_[i].reserve(MAX_NODE_PER_BINS);
            history_hypothesis_index_[i].reserve(MAX_NODE_PER_BINS);
        }

        clear_all(); // 初始化数据结构，清空所有数据

        // 误差分布表格，依据先验THETA和RHO的SIGMA计算获得
        error_distribution_table_.resize(MAX_BINS);
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            // 计算BINS对应的x,y坐标
            size_t x_index = i / LOGIC_BASED_NUM_Y_BINS;
            size_t y_index = i % LOGIC_BASED_NUM_Y_BINS;
            double x = (x_index + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_X / LOGIC_BASED_NUM_X_BINS) - LOGIC_BASED_MAX_ABS_X;
            double y = (y_index + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_Y / LOGIC_BASED_NUM_Y_BINS) - LOGIC_BASED_MAX_ABS_Y;

            // 转换为极坐标
            double rho = std::sqrt(x * x + y * y);
            double theta = std::atan2(y, x); // 弧度，范围 -π 到 π

            // 常数
            double sigma_r = track_project::base_r_sigma_km;
            double sigma_theta_rad = track_project::base_theta_sigma_deg * M_PI / 180.0; // 度转弧度

            // 计算x方向上的方差，\sigma^2_x=sigma^2_r \times cos^2(theta) + rho^2 \times sigma^2_r \times sin^2(theta)
            double sigma_x2 = std::pow(sigma_r * std::cos(theta), 2) +
                              std::pow(rho * sigma_theta_rad * std::sin(theta), 2);
            //\sigma^2_y=sigma^2_r \times sin^2(theta) + rho^2 \times sigma^2_r \times cos^2(theta)
            double sigma_y2 = std::pow(sigma_r * std::sin(theta), 2) +
                              std::pow(rho * sigma_theta_rad * std::cos(theta), 2);
            double sigma_x = std::sqrt(sigma_x2);
            double sigma_y = std::sqrt(sigma_y2);

            error_distribution_table_[i] = std::make_pair(sigma_x, sigma_y);

            // 在循环里加个调试
            // 找 rho > 400km 的点
            if (rho > 300)
            {
                std::cout << "rho=" << rho << ", theta=" << theta * 180 / M_PI
                          << "°, sigma_x=" << sqrt(sigma_x2)
                          << ", sigma_y=" << sqrt(sigma_y2) << std::endl;
            }
        }
    }

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

        return ProcessStatus::SUCCESS;
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

        // TODO，依据DOPPLER反推可能的波门大小，从而确定需要查询的bin范围，目前先暂时扩大范围，之后再试试根据误差分布表格计算合理的范围
        // 后面的索引还要改，不是简单的获取一次波们就可以了，得先依据所有可能的波们，再结合误差分布表格计算出所有可能的bin范围，最后查询这些bin中的假设节点
        // 所以逻辑错了，明天改吧

        // 计算分辨率
        double x_bin_size = (2 * LOGIC_BASED_MAX_ABS_X) / LOGIC_BASED_NUM_X_BINS; // X轴每个bin的大小
        double y_bin_size = (2 * LOGIC_BASED_MAX_ABS_Y) / LOGIC_BASED_NUM_Y_BINS; // Y轴每个bin的大小

        // 离散化输入坐标，四舍五入，计算对应的bin索引
        size_t x_index = static_cast<size_t>(std::round((x + LOGIC_BASED_MAX_ABS_X) / x_bin_size));
        size_t y_index = static_cast<size_t>(std::round((y + LOGIC_BASED_MAX_ABS_Y) / y_bin_size));

        // 先查表查看当前点的误差分布情况，根据误差分布情况计算需要查询的bin范围
        size_t x_bin_range = error_distribution_table_[x_index * LOGIC_BASED_NUM_Y_BINS + y_index].first;  // 误差分布表格中存储的x方向误差范围，单位bin数量
        size_t y_bin_range = error_distribution_table_[x_index * LOGIC_BASED_NUM_Y_BINS + y_index].second; // 误差分布表格中存储的y方向误差范围，单位bin数量

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
                if (query_x_index >= LOGIC_BASED_NUM_X_BINS || query_y_index >= LOGIC_BASED_NUM_Y_BINS)
                {
                    continue; // 索引越界，跳过
                }

                // 计算bin索引
                size_t bin_index = query_x_index * LOGIC_BASED_NUM_Y_BINS + query_y_index;
                // 将该bin中的假设节点加入候选列表
                candidate_nodes.insert(candidate_nodes.end(), current_hypothesis_index_[bin_index].begin(), current_hypothesis_index_[bin_index].end());
            }
        }
        return candidate_nodes;
    }

    ProcessStatus LogicBasedInitiator::extend_hypotheses(const TrackPoint &point, const std::vector<HypothesisNode *> &parent_nodes)
    {
        // TODO
        return ProcessStatus::SUCCESS;
    }

    void LogicBasedInitiator::clear_all()
    {
        // 预留空间
        for (size_t i = 0; i < 4; ++i)
        {
            point_batches_[i].clear();
            hypothesis_layers_[i].clear();
            batch_timestamps_[i] = Timestamp();
        }

        // 初始化索引表，避免内存碎片，一次性开到位
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            current_hypothesis_index_[i].clear();
            history_hypothesis_index_[i].clear();
        }

        // 误差分布表格不清空，因为这个不随数据到来二变化，只收到build_error_distribution_table的调用才会更新
    }
}