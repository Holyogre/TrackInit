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
            point_batches_[i].reserve(1000);                                         // 每批次预留1000个点迹的空间
            hypothesis_layers_[i].reserve(MAX_BINS * LOGIC_BASED_MAX_NODE_PER_BINS); // 避免内存碎片一次开到位
        }

        // 初始化索引表，避免内存碎片一次性开到位
        current_hypothesis_index_.resize(MAX_BINS);
        history_hypothesis_index_.resize(MAX_BINS);
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            current_hypothesis_index_[i].reserve(LOGIC_BASED_MAX_NODE_PER_BINS);
            history_hypothesis_index_[i].reserve(LOGIC_BASED_MAX_NODE_PER_BINS);
        }

        clear_all(); // 初始化数据结构，清空所有数据

        // 误差分布表格，依据先验THETA和RHO的SIGMA计算获得
        error_distribution_table_.resize(MAX_BINS);
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            // 计算BINS对应的x,y坐标
            auto [x_index, y_index] = bin_index_to_xy_index(i);
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
        shift_batches_and_hypotheses(points);
        point_batches_[0] = points;             // 存储最新批次数据
        timestamp_batches_[0] = points[0].time; // 存储最新批次时间戳

        // 清空输出航迹
        new_tracks.clear();

        for (const auto &point : point_batches_[0])
        {
            // 查询满足条件的假设节点
            auto candidate_nodes_ = make_nodes_by_points(point);

            if (candidate_nodes_.empty())
            {
                continue; // 没有满足条件的假设节点，继续处理下一个点迹
            }

            // 预处理假设节点，对于无可置疑的假设节点直接拓展假设，或直接输出为航迹
            ProcessStatus status = preprocess_nodes(candidate_nodes_, new_tracks);
            if (status != ProcessStatus::SUCCESS)
            {
                return status; // 如果预处理失败，直接返回错误码
            }
        }

        return ProcessStatus::SUCCESS;
    }

    // 使用swap交换位置，清空最新批次位置。。这步不会涉及到内存分配和释放
    void LogicBasedInitiator::shift_batches_and_hypotheses(const std::vector<TrackPoint> &new_points)
    {
        // 1. 所有数据后移动，最后一批数据放到最前面来
        for (int i = 3; i > 0; --i)
        {
            std::swap(point_batches_[i], point_batches_[i - 1]);
            std::swap(hypothesis_layers_[i], hypothesis_layers_[i - 1]);
            std::swap(timestamp_batches_[i], timestamp_batches_[i - 1]);
        }
        std::swap(current_hypothesis_index_, history_hypothesis_index_);

        // 2. 清空最新的批次位置（索引0）
        point_batches_[0].clear();
        hypothesis_layers_[0].clear();
        for (auto &bin : current_hypothesis_index_)
        {
            bin.clear();
        }

        // 放置数据
        point_batches_[0] = new_points; // 存储最新批次数据
    }

    // 查询假设节点，输入当前点迹的经纬度和DOPPLER，经过反推查询所有满足要求的假设树节点
    std::vector<LogicBasedInitiator::HypothesisNode> LogicBasedInitiator::make_nodes_by_points(const TrackPoint &point) const
    {
        // 参数提取
        double x = point.x;
        double y = point.y;
        size_t bin_index = location_to_bin_index(x, y);
        double doppler = point.doppler;                                                                                    // m/s
        double dt = static_cast<double>(timestamp_batches_[0].milliseconds - timestamp_batches_[1].milliseconds) / 1000.0; // s
        double heading_resolution_rad = LOGIC_BASED_HEADING_RESOLUTION_DEG * M_PI / 180.0;                                 // 航向分辨率，单位弧度
        // sigma_x,sigma_y参数提取
        auto [sigma_x, sigma_y] = error_distribution_table_[bin_index];
        double x_protected = LOGIC_BASED_PROTECTIVE_RADIUS_KM + 1.96 * sigma_x; // 保护半径，单位km
        double y_protected = LOGIC_BASED_PROTECTIVE_RADIUS_KM + 1.96 * sigma_y; // 保护半径，单位km

        // ===== 第一步：从多普勒反推可能的航向范围 =====
        auto [heading_center, heading_range] = calculate_heading_range(x, y, doppler);
        double heading_start = heading_center - heading_range;
        double heading_end = heading_center + heading_range;

        // ===== 第二步：收集所有航向产生的prev点，并记录边界值 =====
        // 用于记录外推点迹的边界情况
        struct prev_info
        {
            size_t prev_bin_index;
            double x, y;
            double heading_start, heading_end;
            double bias_x_max, bias_y_max;
            double doppler_min, doppler_max;
        };
        double prev_x_min = std::numeric_limits<double>::max();
        double prev_x_max = std::numeric_limits<double>::lowest();
        double prev_y_min = std::numeric_limits<double>::max();
        double prev_y_max = std::numeric_limits<double>::lowest();
        std::vector<prev_info> prev_info_list; // 记录涉及的prev_infonfo，包含bin索引和doppler范围
        size_t prev_bin = MAX_BINS + 1;        // 用于记录前一个bin的参数
        // 依据heading，计算doppler和边界情况
        for (double heading = heading_start; heading <= heading_end + heading_resolution_rad; heading += heading_resolution_rad)
        {
            // 计算航向对应的速度分量
            double vx = track_project::velocity_max * std::sin(heading); // 北偏东角度，sin对应x分量
            double vy = track_project::velocity_max * std::cos(heading); // 北偏东角度，cos对应y分量

            // 反推prev点位置
            double prev_x = x - vx * dt / 1000.0; // km
            double prev_y = y - vy * dt / 1000.0; // km

            // 计算doppler，指向雷达站为正方向
            double prev_doppler = -(vx * prev_x + vy * prev_y) / hypot(prev_x, prev_y); // m/s

            // 计算prev点对应的bin索引
            auto prev_bin_index = location_to_bin_index(prev_x, prev_y);

            // 若无新点，存放新点，同时更新边界值
            if (prev_bin_index != prev_bin)
            {
                prev_bin = prev_bin_index;
                auto bin_error = error_distribution_table_[prev_bin_index]; // 该bin对应的误差分布参数
                double temp_bias_x_max = x_protected + 1.96 * bin_error.first;
                double temp_bias_y_max = y_protected + 1.96 * bin_error.second;
                double heading_start = heading;
                double heading_end = heading;
                prev_info_list.push_back({prev_bin_index, prev_x, prev_y, heading_start, heading_end, temp_bias_x_max, temp_bias_y_max, prev_doppler, prev_doppler});
            }
            else // 该点已经存在，更新边界值
            {
                auto &temp_prev_info = prev_info_list.back();
                temp_prev_info.x = (prev_x + temp_prev_info.x) / 2.0; // 近似取中心点，一般来说一个波们差不多也就遍历个1-3次，无所谓了
                temp_prev_info.y = (prev_y + temp_prev_info.y) / 2.0;
                temp_prev_info.heading_start = std::min(temp_prev_info.heading_start, heading);
                temp_prev_info.heading_end = std::max(temp_prev_info.heading_end, heading);
                temp_prev_info.bias_x_max = std::max(temp_prev_info.bias_x_max, x_protected + 1.96 * error_distribution_table_[prev_bin_index].first);
                temp_prev_info.bias_y_max = std::max(temp_prev_info.bias_y_max, y_protected + 1.96 * error_distribution_table_[prev_bin_index].second);
                temp_prev_info.doppler_min = std::min(temp_prev_info.doppler_min, prev_doppler);
                temp_prev_info.doppler_max = std::max(temp_prev_info.doppler_max, prev_doppler);
            }

            // 更新搜索区域边界
            prev_x_min = std::min(prev_x_min, prev_x);
            prev_x_max = std::max(prev_x_max, prev_x);
            prev_y_min = std::min(prev_y_min, prev_y);
            prev_y_max = std::max(prev_y_max, prev_y);
        }

        // ===== 第三步：构造搜索区域，遍历搜索区域，提取符合要求的点迹 =====
        std::array<size_t, 4> bin_index_range = calculate_bin_index_range(prev_x_min, prev_x_max, prev_y_min, prev_y_max, x_protected, y_protected);
        std::vector<HypothesisNode> candidate_nodes;                                        // 存储满足条件的假设节点
        for (size_t x_index = bin_index_range[0]; x_index <= bin_index_range[1]; ++x_index) // 为了确保假设只被遍历一次
        {
            for (size_t y_index = bin_index_range[2]; y_index <= bin_index_range[3]; ++y_index)
            {
                size_t hypothesis_bin_index = xy_index_to_bin_index(x_index, y_index);

                // 遍历该bin中的假设节点
                for (HypothesisNode *node : history_hypothesis_index_[hypothesis_bin_index]) // 节点是稀疏的，绝大多数情况这个for循环不会执行
                {
                    TrackPoint hypothesis_point = *(node->associated_point); // 假设节点关联的点迹
                    double hypothesis_heading_start = node->heading_start;
                    double hypothesis_heading_end = node->heading_end;
                    std::pair<double, double> temp_heading_range(1e10, -1e10);

                    for (const auto &prev_info_cell : prev_info_list) // 至多30次循环
                    {
                        // 1. 判断doppler是否在范围内
                        if (hypothesis_point.doppler < prev_info_cell.doppler_min - LOGIC_BASED_PROTECTIVE_DOPPLER_M_S ||
                            hypothesis_point.doppler > prev_info_cell.doppler_max + LOGIC_BASED_PROTECTIVE_DOPPLER_M_S)
                        {
                            continue; // 不满足doppler条件，跳过该节点
                        }

                        // 2. 判断位置是否在保护半径范围内
                        if (std::fabs(hypothesis_point.x - prev_info_cell.x) > prev_info_cell.bias_x_max ||
                            std::fabs(hypothesis_point.y - prev_info_cell.y) > prev_info_cell.bias_y_max)
                        {
                            continue; // 不满足位置条件，跳过该节点
                        }

                        // 计算当前航向可能范围的并集
                        temp_heading_range.first = std::min(temp_heading_range.first, prev_info_cell.heading_start);
                        temp_heading_range.second = std::max(temp_heading_range.second, prev_info_cell.heading_end);
                    }

                    if (temp_heading_range.first > 1e5) // 快速退出
                    {
                        continue;
                    }

                    // 和历史假设取并集
                    hypothesis_heading_start = std::min(temp_heading_range.first, hypothesis_heading_start);
                    hypothesis_heading_end = std::max(temp_heading_range.second, hypothesis_heading_end);

                    if (hypothesis_heading_end < hypothesis_heading_start) // 不存在重合航向，舍去该点
                    {
                        continue; // 不满足航向条件，跳过该节点
                    }
                    // 满足条件，加入候选节点列表
                    size_t temp_depth = node->depth + 1;
                    candidate_nodes.emplace_back(temp_depth, bin_index, &point, node, hypothesis_heading_start, hypothesis_heading_end);
                }
            }
        }

        // 如果没有匹配假设，生成新假设
        if (candidate_nodes.empty())
        {
            candidate_nodes.emplace_back(0, bin_index, &point, nullptr, heading_start, heading_end);
        }

        return candidate_nodes;
    }

    // 输出航向中心，和上下波动范围,单位弧度,北偏东
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
        double abs_doppler_m_s = std::abs(doppler); // 为km和s，避免数值不稳定

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
        size_t bin_index = xy_index_to_bin_index(x_index, y_index);

        return bin_index;
    }

    // 保护范围实际上应该是另一个点的sigma值，加上宏定义的边界值，计算该波门的真实边界值
    std::array<size_t, 4> LogicBasedInitiator::calculate_bin_index_range(double x_min, double x_max, double y_min, double y_max,
                                                                         double x_protected, double y_protected) const
    {
        // 搜索sigma值,min,max都在一个波们内，用哪个都无所谓
        auto [sigma_x_1, sigma_y_1] = error_distribution_table_[location_to_bin_index(x_min, y_min)];
        auto [sigma_x_2, sigma_y_2] = error_distribution_table_[location_to_bin_index(x_max, y_max)];
        auto [sigma_x_3, sigma_y_3] = error_distribution_table_[location_to_bin_index(x_min, y_max)];
        auto [sigma_x_4, sigma_y_4] = error_distribution_table_[location_to_bin_index(x_max, y_min)];
        double sigma_x = std::max({sigma_x_1, sigma_x_2, sigma_x_3, sigma_x_4});
        double sigma_y = std::max({sigma_y_1, sigma_y_2, sigma_y_3, sigma_y_4});

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

    ProcessStatus LogicBasedInitiator::extend_hypotheses(const std::vector<HypothesisNode> &node)
    {
        if (hypothesis_layers_[0].size() + node.size() >= MAX_BINS * LOGIC_BASED_NUM_Y_BINS * LOGIC_BASED_MAX_NODE_PER_BINS ||
            node.size() >= LOGIC_BASED_MAX_NODE_PER_BINS)
        {
            // 假设节点总数超过上限，返回错误状态
            return ProcessStatus::TOO_MANY_HYPOTHSIS;
        }

        size_t node_size = node.size();
        if (node_size == 0)
        {
            return ProcessStatus::SUCCESS;
        }

        // 遍历node列表，往前搜索一个假设，若前置假设来源于同一个点，存放到冲突列表中
        struct temp_conflict
        {
            TrackPoint *prev_point_pos;
            std::vector<size_t> node_index;
        };
        std::vector<temp_conflict> conflict_node_list;
        for (size_t index = 0; index < node_size; index++)
        {
            const HypothesisNode &possible_node = node[index];
            switch (possible_node.depth)
            {
            case 0: // 第零层假设必然不冲突，直接放入对应bin就是了
            case 1: // 第一层假设必然不冲突，因为前置点迹必然不同
            {
                hypothesis_layers_[0].push_back(possible_node);
                size_t bin_index = possible_node.bin_index;
                current_hypothesis_index_[bin_index].push_back(&hypothesis_layers_[0].back()); // 存放索引
                break;
            }
            case 2: // 二、三层假设可能冲突，放入冲突区域
            case 3:
            {
                bool is_conflict = false;
                const TrackPoint *possible_node_prev_point_pos = possible_node.parent_node->associated_point;
                for (auto &conflict_node : conflict_node_list)
                {
                    if (conflict_node.prev_point_pos == possible_node_prev_point_pos)
                    {
                        is_conflict = true;
                        conflict_node.node_index.push_back(index);
                        break;
                    }
                }
                if (!is_conflict)
                {
                    conflict_node_list.emplace_back(possible_node_prev_point_pos, index);
                }
                break;
            }
            default:
            {
                return ProcessStatus::WRONG_PREV_NODE;
                break;
            }
            }
        }

        //遍历冲突列表，对于无多个假设项，直接存入，对于多个假设项，存入冲突假设列表
        //todo

        return ProcessStatus::SUCCESS;
    }

    void LogicBasedInitiator::clear_all()
    {
        // 清空存储信息
        for (size_t i = 0; i < 4; ++i)
        {
            point_batches_[i].clear();
            hypothesis_layers_[i].clear();
            timestamp_batches_[i] = Timestamp();
        }

        // 名为vector实际作为array使用的区域，需要逐个清空
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            current_hypothesis_index_[i].clear();
            history_hypothesis_index_[i].clear();
        }

        // 误差分布表格不清空，因为这个不随数据到来二变化，只收到build_error_distribution_table的调用才会更新
    }
}