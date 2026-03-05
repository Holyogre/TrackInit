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

    // 查询假设节点，输入当前点迹的经纬度和DOPPLER，经过反推查询所有满足要求的假设树节点
    std::vector<LogicBasedInitiator::HypothesisNode *> LogicBasedInitiator::query_nodes_by_points(const TrackPoint &point) const
    {
        // 参数提取
        double x = point.x;
        double y = point.y;
        double doppler = point.doppler;                                                                                    // m/s
        double dt = static_cast<double>(timestamp_batches_[0].milliseconds - timestamp_batches_[1].milliseconds) / 1000.0; // s
        double heading_resolution_rad = LOGIC_BASED_HEADING_RESOLUTION_DEG * M_PI / 180.0;                                 // 航向分辨率，单位弧度
        // sigma_x,sigma_y参数提取
        auto [sigma_x, sigma_y] = error_distribution_table_[location_to_bin_index(x, y)];
        double x_protected = LOGIC_BASED_PROTECTIVE_RADIUS_KM + sigma_x; // 保护半径，单位km
        double y_protected = LOGIC_BASED_PROTECTIVE_RADIUS_KM + sigma_y; // 保护半径，单位km

        // ===== 第一步：从多普勒反推可能的航向范围 =====
        std::pair<double, double> heading_center_and_range = calculate_heading_range(x, y, doppler);

        // ===== 第二步：收集所有航向产生的prev点，并记录边界值 =====
        struct xy_range
        {
            double x_min = std::numeric_limits<double>::max();    // 8字节
            double x_max = std::numeric_limits<double>::lowest(); // 8字节
            double y_min = std::numeric_limits<double>::max();    // 8字节
            double y_max = std::numeric_limits<double>::lowest(); // 8字节
        };
        std::unordered_map<size_t, xy_range> bin_indeces; // 记录涉及的bin索引，避免后续遍历全部bin
        for (double heading = heading_center_and_range.first - heading_center_and_range.second;
             heading <= heading_center_and_range.first + heading_center_and_range.second + heading_resolution_rad;
             heading += heading_resolution_rad) // 每10度一个航向
        {
            // 计算航向对应的速度分量
            double vx = track_project::velocity_max * std::sin(heading); // 北偏东角度，sin对应x分量
            double vy = track_project::velocity_max * std::cos(heading); // 北偏东角度，cos对应y分量

            // 反推prev点位置
            double prev_x = x - vx * dt / 1000.0; // km
            double prev_y = y - vy * dt / 1000.0; // km

            // 计算prev点对应的bin索引
            size_t bin_index = location_to_bin_index(prev_x, prev_y);
            bin_indeces[bin_index].x_min = std::min(bin_indeces[bin_index].x_min, prev_x);
            bin_indeces[bin_index].x_max = std::max(bin_indeces[bin_index].x_max, prev_x);
            bin_indeces[bin_index].y_min = std::min(bin_indeces[bin_index].y_min, prev_y);
            bin_indeces[bin_index].y_max = std::max(bin_indeces[bin_index].y_max, prev_y);
        }

        // ===== 第三步：依据sigma值，去重，遍历假设区域，将假设输出 =====
        // 第一个区域采用遍历的方式
        auto prev_bin_pos = bin_indeces.begin();
        assert(prev_bin_pos != bin_indeces.end() && "不存在满足条件的bin索引");
        std::array<size_t, 4> prev_index_range = calculate_bin_index_range(prev_bin_pos->second.x_min, prev_bin_pos->second.x_max,
                                                                           prev_bin_pos->second.y_min, prev_bin_pos->second.y_max,
                                                                           x_protected, y_protected);
        std::vector<HypothesisNode *> candidate_nodes = extrack_hypothesis_from_bins(prev_index_range); // 构建第一批返回点

        // 从第二个区域开始，仅计算未重复的bin索引，避免重复计算
        auto current_bin_pos = std::next(prev_bin_pos);
        for (current_bin_pos; current_bin_pos != bin_indeces.end(); ++current_bin_pos)
        {
            // 依据sigma_x,sigma_y，保护半径计算一个std::array<size_t,4>，来表示满足条件的bin范围
            std::array<size_t, 4> current_index_range = calculate_bin_index_range(current_bin_pos->second.x_min, current_bin_pos->second.x_max,
                                                                                  current_bin_pos->second.y_min, current_bin_pos->second.y_max,
                                                                                  x_protected, y_protected);

            // 计算新增部分
            std::vector<std::array<size_t, 4>> rect_differences = get_expansion_regions(prev_index_range, current_index_range);

            for (const auto recv_item : rect_differences)
            {
                std::vector<HypothesisNode *> nodes = extrack_hypothesis_from_bins(recv_item);
                candidate_nodes.insert(candidate_nodes.end(), nodes.begin(), nodes.end());
            }

            // 迭代，保留上一个迭代的结果
            prev_index_range = current_index_range;
            prev_bin_pos = current_bin_pos;
        }

        return candidate_nodes;
    }

    // 输出航向中心，和上下波动范围
    std::pair<double, double> LogicBasedInitiator::calculate_heading_range(double x, double y, double doppler) const
    {
        // 计算距离和视线方向
        double r = std::sqrt(x * x + y * y);
        assert(r > 1 && "存在异常点迹，检测输入端");

        // 船只相对于雷达站的北偏东角度（弧度）
        double los_angle = std::atan2(x, y); // atan2(x, y)得到相对于北的角度

        // 中心航向
        double heading_center = (doppler < 0) ? los_angle : los_angle + M_PI; // 以视线方向作为中心

        // 获取abs_doppler
        double abs_doppler_m_s = std::abs(doppler) / 1000.0; // 转换为km/s，单位统一为km和s，避免数值不稳定

        // 如果Doppler足够小直接给半圆范围，避免数值不稳定
        if (abs_doppler_m_s < 1e-6)
        {
            // 可能的航向范围是个半圆
            double heading_range = M_PI / 2; // 单一方向

            return std::make_pair(heading_center, heading_range);
        }

        // 对于过大多普勒只给一个方向
        if (abs_doppler_m_s > track_project::velocity_max - 1e-3)
        {
            double heading_range = 0.0; // 单一方向

            return std::make_pair(heading_center, heading_range);
        }

        // 最小可能的夹角（对应最大船速）
        double cos_aspect_angle = abs_doppler_m_s / track_project::velocity_max;
        cos_aspect_angle = std::clamp(cos_aspect_angle, 0.0, 1.0); // 确保在有效范围内
        double aspect_angle = std::acos(cos_aspect_angle);         // 航向与径向的夹角

        return std::make_pair(heading_center, aspect_angle);
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

    // 保护范围实际上应该是另一个点的sigma值，加上宏定义的边界值，计算该波门的真实边界值
    std::array<size_t, 4> LogicBasedInitiator::calculate_bin_index_range(double x_min, double x_max, double y_min, double y_max,
                                                                         double x_protected, double y_protected) const
    {
        // 搜索sigma值,min,max都在一个波们内，用哪个都无所谓
        auto [sigma_x, sigma_y] = error_distribution_table_[location_to_bin_index(x_min, y_min)];

        // 累计保护半径，计算边界值
        double x_min_protected = x_min - sigma_x - x_protected;
        double x_max_protected = x_max + sigma_x + x_protected;
        double y_min_protected = y_min - sigma_y - y_protected;
        double y_max_protected = y_max + sigma_y + y_protected;

        // 计算离散索引，边界条件由location_to_xy_index函数内部处理
        auto [x_min_index, y_min_index] = location_to_xy_index(x_min_protected, y_min_protected);
        auto [x_max_index, y_max_index] = location_to_xy_index(x_max_protected, y_max_protected);

        return {x_min_index, x_max_index, y_min_index, y_max_index};
    }

    // 添加假设点

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