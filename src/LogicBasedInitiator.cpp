#include "LogicBasedInitiator.hpp"

#include <cmath>
#include "unordered_set"

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
        point_batches_[0] = points;             // 存储最新批次数据
        timestamp_batches_[0] = points[0].time; // 存储最新批次时间戳

        // 清空输出航迹
        new_tracks.clear();

        // ==================== STEP1: 扩展假设树 ====================
        for (auto point : points)
        {
            // 查询满足条件的假设节点
            auto candidate_nodes_ = query_nodes_by_points(point);

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
            timestamp_batches_[i] = timestamp_batches_[i - 1];
        }

        // 4. 清空最新的批次位置（索引0）
        point_batches_[0].clear();
        hypothesis_layers_[0].clear();
    }

    // 只有DOPPLER参数，因此需要外推所有可能的X,Y所在位置，然后依据误差分布函数进一步扩大搜索范围，最后合并重复假设
    // doppler粗筛，位置精筛
    std::vector<LogicBasedInitiator::HypothesisNode *> LogicBasedInitiator::query_nodes_by_points(const TrackPoint &point) const
    {
        // 参数提取
        double x = point.x;             // 单位km
        double y = point.y;             // 单位km
        double doppler = point.doppler; // 单位m/s

        // 计算可能的航向范围
        std::pair<double, double> heading_range = calculate_heading_range(x, y, doppler);
        heading_range.first += M_PI; // 由于是倒推前一个点，所以要反向
        heading_range.second += M_PI;

        // 计算当前点迹的sigma_x和sigma_y，航迹起始阶段不确定性太大，尽量保证检测率，虚警率交给后一步判断
        auto [sigma_x, sigma_y] = error_distribution_table_[location_to_bin_index(x, y)];
        double search_radius_x = 2.58 * sigma_x; // 2.58SIGMA对应99.99%的概率
        double search_radius_y = 2.58 * sigma_y; // 2.58SIGMA对应99.99%的概率

        // 计算时间
        double dt = static_cast<double>(timestamp_batches_[0].milliseconds - timestamp_batches_[1].milliseconds) / 1000.0;

        // 计算依据对应航向所有可能的x_index,y_index列表，多搜索一圈以防万一,步长为离散分辨率,heading的单位为rad,北偏东
        std::unordered_set<size_t> bin_indices;                                            // 使用unordered_set自动去重
        double heading_resolution_rad = LOGIC_BASED_HEADING_RESOLUTION_DEG * M_PI / 180.0; // 航向离散分辨率，单位弧度
        for (double heading = heading_range.first; heading <= heading_range.second + heading_resolution_rad; heading += heading_resolution_rad)
        {
            // 计算对应的x,y位置
            double vx = track_project::velocity_max * std::sin(heading);
            double vy = track_project::velocity_max * std::cos(heading);
            double doppler_calculated = (vx * x + vy * y) / std::sqrt(x * x + y * y); // 计算该航向下的doppler值

            // 依据航向计算前一个假设的基准XY坐标
            double prev_x = x - vx * dt; // 依据航向和速度计算前一个假设的XY坐标
            double prev_y = y - vy * dt;

            // 计算对应的XY的SIGMA范围，扩大搜索范围，实际上标准的做法是用更大的区间反向搜索是否符合SIGMA范围，因为本算法中，误差分布函数与坐标相关
            auto [prev_sigma_x, prev_sigma_y] = error_distribution_table_[location_to_bin_index(prev_x, prev_y)];
            double search_radius_pre_x = 2.58 * prev_sigma_x; // 2.58SIGMA对应99.99%的概率
            double search_radius_pre_y = 2.58 * prev_sigma_y; // 2.58SIGMA对应99.99%的概率

            // 依据当前点迹的不确定度和前一个假设的不确定度扩大搜索范围
            double total_search_radius_x = search_radius_x + search_radius_pre_x + LOGIC_BASED_PROTECTIVE_RADIUS_KM;
            double total_search_radius_y = search_radius_y + search_radius_pre_y + LOGIC_BASED_PROTECTIVE_RADIUS_KM;
            double x_min = prev_x - total_search_radius_x;
            double x_max = prev_x + total_search_radius_x;
            double y_min = prev_y - total_search_radius_y;
            double y_max = prev_y + total_search_radius_y;
            auto [x_index_min, y_index_min] = location_to_xy_index(x_min, y_min);
            auto [x_index_max, y_index_max] = location_to_xy_index(x_max, y_max);

            // 加入满足条件的bin索引，用unordered_set自动去重，确保每个bin索引唯一
            for (size_t x_index = x_index_min; x_index <= x_index_max; ++x_index)
            {
                for (size_t y_index = y_index_min; y_index <= y_index_max; ++y_index)
                {
                    size_t index = x_index * LOGIC_BASED_NUM_Y_BINS + y_index;
                    bin_indices.insert(index);
                }
            }
        }

        // 依据doppler进一步进行筛选，只保留满足条件的假设节点
        std::vector<HypothesisNode *> candidate_nodes_; // 候选节点列表
        for (auto bin_index : bin_indices)
        {
            // 获取该bin索引对应的假设节点列表，加入候选节点列表
            const auto &nodes_in_bin = history_hypothesis_index_[bin_index];

            // 遍历假设节点信息
            for (auto node : nodes_in_bin)
            {
                if (doppler >= node->vr_min && doppler <= node->vr_max)
                {
                    candidate_nodes_.push_back(node);
                }
            }
        }

        return candidate_nodes_;
    }

    // 依据x,y和doppler计算船只的航向范围
    std::pair<double, double> LogicBasedInitiator::calculate_heading_range(double x, double y, double doppler) const
    {
        // 依据doppler和velocity_max计算所有可能的航向角度序列
        double abs_doppler = std::abs(doppler); // 取绝对值，单位m/s
        if (abs_doppler > track_project::velocity_max)
        {
            return {}; // 速度超过最大值，跳过该点迹
        }

        // 计算基准航向，依据doppler的正负确定是朝向还是背向，为极坐标表示方法，便于调用math库函数计算
        double line_of_sight_rad = std::atan2(y, x); // 观测方向（极坐标）
        double base_dir_rad = 0.0;
        if (doppler > 0) // 表示靠近还是原理，以此决定是否加PI
        {
            base_dir_rad = line_of_sight_rad + M_PI;
            if (base_dir_rad >= 2 * M_PI)
            {
                base_dir_rad -= 2 * M_PI;
            }
        }
        else // 不然，值域为(-pi,pi)，需要归一化到(0,2*pi)
        {
            base_dir_rad = line_of_sight_rad;
            if (base_dir_rad < 0)
            {
                base_dir_rad += 2 * M_PI;
            }
        }

        // 计算可能的速度方向与视线方向的夹角,abs_doppler不可能为0，此时检测不出来
        double angle_rad = std::acos(abs_doppler / track_project::velocity_max); // 计算夹角，弧度，范围[0,pi/2)
        angle_rad = std::clamp(angle_rad, 1e-4, M_PI / 2.0 - 1e-4);              // 确保夹角在合理范围内

        double heading1 = base_dir_rad - angle_rad;
        double heading2 = base_dir_rad + angle_rad;

        return std::make_pair(heading1, heading2);
    }

    // 计算x,y坐标对应的离散x,y索引
    std::pair<size_t, size_t> LogicBasedInitiator::location_to_xy_index(double x, double y) const
    {
        //  计算分辨率
        double x_bin_size = (2 * LOGIC_BASED_MAX_ABS_X) / LOGIC_BASED_NUM_X_BINS; // X轴每个bin的大小
        double y_bin_size = (2 * LOGIC_BASED_MAX_ABS_Y) / LOGIC_BASED_NUM_Y_BINS; // Y轴每个bin的大小

        // 离散化输入坐标，四舍五入，计算对应的x,y索引
        size_t x_index = static_cast<size_t>(std::round((x + LOGIC_BASED_MAX_ABS_X) / x_bin_size));
        size_t y_index = static_cast<size_t>(std::round((y + LOGIC_BASED_MAX_ABS_Y) / y_bin_size));
        x_index = std::clamp(x_index, static_cast<size_t>(0), LOGIC_BASED_NUM_X_BINS - 1);
        y_index = std::clamp(y_index, static_cast<size_t>(0), LOGIC_BASED_NUM_Y_BINS - 1);

        return std::make_pair(x_index, y_index);
    }

    // 计算BIN值对应的索引
    size_t LogicBasedInitiator::location_to_bin_index(double x, double y) const
    {
        // 调用location_to_xy_index获取x_index和y_index
        auto [x_index, y_index] = location_to_xy_index(x, y);
        size_t bin_index = x_index * LOGIC_BASED_NUM_Y_BINS + y_index;

        return bin_index;
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
            timestamp_batches_[i] = Timestamp();
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