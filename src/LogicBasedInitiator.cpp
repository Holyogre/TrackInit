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

        // 初始化冲突假设列表
        for (size_t i = 0; i < MAX_INPUT_POINTS; ++i)
        {
            conflicts_hypothesis_[i].reserve(LOGIC_BASED_MAX_NODE_PER_BINS * 10); // 预留假设空间
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
            double sigma_r = track_project::BASE_R_SIGMA_KM;
            double sigma_theta_rad = track_project::BASE_THETA_SIGMA_DEG * M_PI / 180.0; // 度转弧度

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
        if (points.empty()) // 没数据也不能把之前的结果全删了。。反正靠time更新的
        {
            return ProcessStatus::NO_POINT; // 没有点迹，这批次数据不处理，不影响后续关联
        }

        // 移动所有数据，准备进行下一批次航迹起始
        rotate_to_new_batch(points, new_tracks);

        // 逐个点迹处理假设
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

        // 统一处理冲突假设
        auto status = reprocess_conflict_nodes(new_tracks);
        if (status != ProcessStatus::SUCCESS)
        {
            return status; // 如果预处理失败，直接返回错误码
        }

        // 调用回调函数，输出结果
        assert(trackCallback_ != nullptr && "HoughSlice类调用的时候没有绑定回调函数");
        trackCallback_(new_tracks);

        return ProcessStatus::SUCCESS;
    }

    // 使用swap交换位置，清空最新批次位置。。这步不会涉及到内存分配和释放
    void LogicBasedInitiator::rotate_to_new_batch(const std::vector<TrackPoint> &new_points, std::vector<std::array<TrackPoint, 4>> &new_tracks)
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

        // 3. 清空冲突假设列表
        for (auto &conflict_ : conflicts_hypothesis_)
        {
            conflict_.clear();
        }

        // 4.清空输出航迹
        new_tracks.clear();

        // 4. 放置数据
        point_batches_[0] = new_points;             // 存储最新批次数据
        timestamp_batches_[0] = new_points[0].time; // 存储最新批次时间戳
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
            double sog = std::fabs(doppler / std::cos(heading - heading_center)); // 反推SOG，单位km/s

            // 计算航向对应的速度分量
            double vx = sog * std::sin(heading); // 北偏东角度，sin对应x分量
            double vy = sog * std::cos(heading); // 北偏东角度，cos对应y分量

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
                        constexpr double doppler_tolerance = 1.96 * track_project::BASE_DOPPLER_TOLERANCE_M_S + LOGIC_BASED_PROTECTIVE_DOPPLER_M_S;
                        if (hypothesis_point.doppler < prev_info_cell.doppler_min - doppler_tolerance ||
                            hypothesis_point.doppler > prev_info_cell.doppler_max + doppler_tolerance)
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
        if (abs_doppler_m_s > track_project::VELOCITY_MAX - 1e-3)
        {
            double heading_range = 0.0; // 单一方向

            return std::make_pair(heading_center, heading_range);
        }

        // 最小可能的夹角（对应最大船速）
        double cos_aspect_angle = abs_doppler_m_s / track_project::VELOCITY_MAX;
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

    ProcessStatus LogicBasedInitiator::preprocess_nodes(std::vector<HypothesisNode> &node, std::vector<std::array<TrackPoint, 4>> &new_tracks)
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

        // ===== 第一步：构建冲突组结构 =====
        struct ConflictGroup
        {
            const TrackPoint *prev_point;
            std::vector<size_t> node_indices;
        };
        std::vector<ConflictGroup> conflict_groups_depth2; // 存储冲突组，每组包含相同前置点迹的假设节点索引
        std::vector<ConflictGroup> conflict_groups_depth3; // 存储冲突组，每组包含相同前置点迹的假设节点索引
        conflict_groups_depth2.reserve(node_size);
        conflict_groups_depth3.reserve(node_size);

        // ===== 第二步：遍历node列表，往前搜索一个假设，若前置假设来源于同一个点，存放到冲突列表中 =====
        // 冲突点迹存放逻辑
        auto addToConflictGroup = [](std::vector<ConflictGroup> &groups, const HypothesisNode &possible_node, size_t index)
        {
            const TrackPoint *prev_point = possible_node.parent_node->associated_point;
            auto it = std::find_if(groups.begin(), groups.end(),
                                   [prev_point](const ConflictGroup &g)
                                   { return g.prev_point == prev_point; });
            if (it != groups.end())
            {
                it->node_indices.push_back(index);
            }
            else
            {
                groups.push_back({prev_point, {index}});
            }
        };
        for (size_t index = 0; index < node_size; ++index)
        {
            const HypothesisNode &possible_node = node[index];
            switch (possible_node.depth)
            {
            case 0:
            case 1:
                hypothesis_layers_[0].push_back(possible_node);
                current_hypothesis_index_[possible_node.bin_index].push_back(&hypothesis_layers_[0].back());
                break;
            case 2:
                addToConflictGroup(conflict_groups_depth2, possible_node, index);
                break;
            case 3:
                addToConflictGroup(conflict_groups_depth3, possible_node, index);
                break;
            default:
                return ProcessStatus::WRONG_PREV_NODE;
            }
        }

        // ===== 第三步：遍历冲突列表，对于无唯一假设项，直接存入，对于多个假设项，存入冲突假设列表 =====
        auto addConflictsToHypothesis = [this](const ConflictGroup &group, const std::vector<HypothesisNode> &node)
        {
            const TrackPoint *prev_point = group.prev_point;
            size_t prev_point_index = static_cast<size_t>(prev_point - point_batches_[1].data());
            for (size_t index : group.node_indices)
            {
                conflicts_hypothesis_[prev_point_index].push_back(node[index]);
            }
        };
        // 处理depth=2的冲突假设，对符合标准的假设直接存入假设层，对不符合标准的假设存入冲突假设列表
        for (const auto &group : conflict_groups_depth2)
        {
            if (group.node_indices.size() == 1)
            {
                const HypothesisNode &node_to_add = node[group.node_indices[0]];
                hypothesis_layers_[0].push_back(node_to_add);
                size_t bin_index = node_to_add.bin_index;
                current_hypothesis_index_[bin_index].push_back(&hypothesis_layers_[0].back());
            }
            else
            {
                addConflictsToHypothesis(group, node);
            }
        }

        // 处理depth=3的冲突假设，对于符合标准的假设直接存入new_track，对不符合标准的假设存入冲突假设列表，后续可以添加depth=3特有的处理逻辑
        for (const auto &group : conflict_groups_depth3)
        {
            if (group.node_indices.size() == 1)
            {
                // 前向索引四批次数据，输出为航迹
                const HypothesisNode &node_to_add = node[group.node_indices[0]];
                std::array<TrackPoint, 4> new_track;
                new_track[0] = *(node_to_add.parent_node->parent_node->associated_point); // depth=3的前两个节点关联的点迹
                new_track[1] = *(node_to_add.parent_node->associated_point);
                new_track[2] = *(node_to_add.associated_point);
                new_track[3] = *(node_to_add.associated_point); // depth=3的最后一个节点关联的点迹
                if (filter_init_track_func_(new_track))
                {
                    new_track[3].confidence = 1.0; // 独占航迹置信度默认为1
                    new_tracks.push_back(new_track);
                }
            }
            else
            {
                addConflictsToHypothesis(group, node);
            }
        }

        return ProcessStatus::SUCCESS;
    }

    // 若有多个假设，满足1、当前假设关联点迹相同；2、此前假设关联点迹相同；则只准留下一个假设
    // 若有多个假设，满足1、当前假设关联点迹不同；2、此前假设关联点迹相同；则计算一个置信度，分配前一个点迹的所有权
    ProcessStatus LogicBasedInitiator::reprocess_conflict_nodes(std::vector<std::array<TrackPoint, 4>> &new_tracks)
    {
        struct Conflict_track
        {
            size_t prev_point_index; // 冲突点迹索引
            std::vector<TrackPoint *> current_point_index;
            std::vector<double> score1; // 最重要的评分
            std::vector<double> score2; // 当最重要的评分近似的时候，选用此评分作为依据

            // 同步存放参数的函数，保证每个点迹对应的评分是一一对应的
            void add_candidate(TrackPoint *point, double s1, double s2)
            {
                current_point_index.push_back(point);
                score1.push_back(s1);
                score2.push_back(s2);
            }
        };

        std::vector<Conflict_track> conflict_tracks; // 存储冲突组合，记录分数
        conflict_tracks.reserve(100);                // 最大1000条航迹，存储100条必然够了

        // 冲突列表和点迹列表的索引是一一对应的，故能如此做
        for (size_t point_index = 0; point_index < conflicts_hypothesis_.size(); ++point_index)
        {
            auto &conflict_nodes = conflicts_hypothesis_[point_index];
            if (conflict_nodes.size() <= 1)
            {
                continue;
            }

            // 从前置航迹索引中，存入一个冲突航迹
            Conflict_track temp_conflict_track{
                point_index,
                std::vector<TrackPoint *>(),
                std::vector<double>(),
                std::vector<double>()};
            conflict_tracks.push_back(std::move(temp_conflict_track));
            auto &current_conflict = conflict_tracks.back();

            // 利用associated_point连续的特性进行分组处理
            const TrackPoint *current_point = nullptr;
            double best_score1 = 1e10;
            double best_score2 = 1e10;
            TrackPoint *best_point = nullptr;
            int group_count = 0;
            bool has_processed_group = false; // 标记是否处理过任何组

            // 遍历所有冲突节点
            for (size_t i = 0; i < conflict_nodes.size(); ++i)
            {
                const auto &node = conflict_nodes[i];
                const TrackPoint *node_point = node.associated_point;

                // 计算当前节点的评分
                double score1 = evaluate_motion_consistency(node);
                double score2 = evaluate_heading_consistency(node);

                // 利用数据特性！ associated_point相同的节点一定是连续存储的，且已经按照associated_point排序好了
                if (node_point != current_point)
                {
                    // 处理上一组的结果
                    if (current_point != nullptr)
                    {
                        has_processed_group = true;
                        // 无论组内有多少节点，都把该组对应的最优节点加入
                        current_conflict.add_candidate(best_point, best_score1, best_score2);
                    }

                    // 重置为新的一组
                    current_point = node_point;
                    best_point = const_cast<TrackPoint *>(node_point);
                    best_score1 = score1;
                    best_score2 = score2;
                    group_count = 1;
                }
                else // 如果是重复的点迹，只保留最好的结果
                {
                    group_count++;

                    // 评分规则：先比score1，如果近似再比score2
                    const double EPSILON = 0.1; // 如果距离小于一个检测范围的话，就只保留一次最好的结果
                    if (score1 < best_score1 - EPSILON)
                    {
                        best_score1 = score1;
                        best_score2 = score2;
                        best_point = const_cast<TrackPoint *>(node_point);
                    }
                    else if (std::abs(score1 - best_score1) < EPSILON && score2 < best_score2)
                    {
                        best_score2 = score2;
                        best_point = const_cast<TrackPoint *>(node_point);
                    }
                }
            }

            // 处理最后一组
            if (current_point != nullptr)
            {
                current_conflict.add_candidate(best_point, best_score1, best_score2);
            }

            // 验证：确保有候选点被添加
            assert(!current_conflict.current_point_index.empty() && "冲突点迹没有生成任何候选点");
        }

        // 依据score1和score2计算depth=2时候的置信度，
        // 依据score1和score2和depth=2时候的置信度计算depth=3的时候的置信度，
        for (auto &conflict : conflict_tracks)
        {
            // TODO: 在这里实现置信度计算和筛选逻辑
            // 可以用 conflict.score1, conflict.score2 和 new_tracks 中的信息
        }

        return ProcessStatus::SUCCESS;
    }

    // 评估航向一致性，输入一个假设节点，回溯获取该节点和父节点的航向中心，计算航向变化量，变化量越小越好
    double LogicBasedInitiator::evaluate_heading_consistency(const HypothesisNode &node)
    {
        if (node.depth < 1)
            return 0.0; // 至少需要前一时刻的航向

        // 回溯获取航向序列
        std::vector<double> heading_centers;
        const HypothesisNode *current = &node;
        while (current != nullptr && heading_centers.size() <= node.depth)
        {
            double heading_center = (current->heading_start + current->heading_end) / 2;
            heading_centers.push_back(heading_center);
            current = current->parent_node;
        }

        if (heading_centers.size() < 2)
            return 0.0;

        // 计算航向变化：最新航向 / 旧航向（考虑角度周期性）
        double latest_heading = heading_centers[0];   // 当前节点
        double previous_heading = heading_centers[1]; // 父节点

        // 处理角度周期性（例如从355度变到5度，实际变化只有10度）
        double heading_diff = std::abs(latest_heading - previous_heading);
        heading_diff = std::min(heading_diff, 2 * M_PI - heading_diff);

        // 航向变化越小越好
        return heading_diff;
    }

    // 评估运动一致性
    double LogicBasedInitiator::evaluate_motion_consistency(const HypothesisNode &node)
    {
        if (node.depth < 2)
            return 0.0;

        // 回溯获取点迹序列
        std::vector<const TrackPoint *> point_chain;
        const HypothesisNode *current = &node;
        while (current != nullptr && point_chain.size() <= node.depth)
        {
            if (current->associated_point != nullptr)
            {
                point_chain.push_back(current->associated_point);
            }
            current = current->parent_node;
        }

        if (point_chain.size() < 2)
            return 0.0;

        if (node.depth == 2 && point_chain.size() >= 3)
        {
            double dx1 = point_chain[1]->x - point_chain[0]->x;
            double dy1 = point_chain[1]->y - point_chain[0]->y;
            double dx2 = point_chain[2]->x - point_chain[1]->x;
            double dy2 = point_chain[2]->y - point_chain[1]->y;
            return (dx1 - dx2) * (dx1 - dx2) + (dy1 - dy2) * (dy1 - dy2);
        }
        else if (node.depth >= 3 && point_chain.size() >= 4)
        {
            std::vector<double> dx, dy;
            for (size_t i = 0; i < point_chain.size() - 1; i++)
            {
                dx.push_back(point_chain[i + 1]->x - point_chain[i]->x);
                dy.push_back(point_chain[i + 1]->y - point_chain[i]->y);
            }

            double mean_dx = 0.0, mean_dy = 0.0;
            for (size_t i = 0; i < dx.size(); i++)
            {
                mean_dx += dx[i];
                mean_dy += dy[i];
            }
            mean_dx /= dx.size();
            mean_dy /= dy.size();

            double var_dx = 0.0, var_dy = 0.0;
            for (size_t i = 0; i < dx.size(); i++)
            {
                var_dx += (dx[i] - mean_dx) * (dx[i] - mean_dx);
                var_dy += (dy[i] - mean_dy) * (dy[i] - mean_dy);
            }
            var_dx /= dx.size();
            var_dy /= dy.size();

            return var_dx + var_dy;
        }

        return 0.0;
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

        // 清空冲突列表
        for (auto &conflict_ : conflicts_hypothesis_)
        {
            conflict_.clear();
        }

        // 误差分布表格不清空，因为这个不随数据到来二变化，只收到build_error_distribution_table的调用才会更新
    }
}