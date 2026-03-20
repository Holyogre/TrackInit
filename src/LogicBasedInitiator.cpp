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

        clear_all(); // 初始化数据结构，清空所有数据

        // 误差分布表格，依据先验THETA和RHO的SIGMA计算获得
        error_distribution_table_.resize(MAX_BINS);
        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            // 计算BINS对应的x,y坐标
            auto [x_index, y_index] = bin_index_to_xy_index(i);
            double sbcpp_x = static_cast<double>(x_index);
            double sbcpp_y = static_cast<double>(y_index);
            double x = (sbcpp_x + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_X / LOGIC_BASED_NUM_X_BINS) - LOGIC_BASED_MAX_ABS_X;
            double y = (sbcpp_y + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_Y / LOGIC_BASED_NUM_Y_BINS) - LOGIC_BASED_MAX_ABS_Y;

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
            ProcessStatus status = process_nodes(candidate_nodes_, new_tracks);
            if (status != ProcessStatus::SUCCESS)
            {
                return status; // 如果预处理失败，直接返回错误码
            }
        }

        // 统一处理冲突假设
        // auto status = reprocess_conflict_nodes(new_tracks);
        // if (status != ProcessStatus::SUCCESS)
        // {
        //     return status; // 如果预处理失败，直接返回错误码
        // }

        // 调用回调函数，输出结果
        assert(trackCallback_ != nullptr && "HoughSlice类调用的时候没有绑定回调函数");
        trackCallback_(new_tracks);

        return ProcessStatus::SUCCESS;
    }

    // 使用swap交换位置，清空最新批次位置。。这步不会涉及到内存分配和释放
    void LogicBasedInitiator::rotate_to_new_batch(const std::vector<TrackPoint> &new_points, std::vector<std::array<TrackPoint, 4>> &new_tracks)
    {
        // 1. 所有数据后移动，最后一批数据放到最前面来
        for (size_t i = 3; i > 0; --i)
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

        // 3.清空输出航迹
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
                    hypothesis_heading_start = std::max(temp_heading_range.first, hypothesis_heading_start);
                    hypothesis_heading_end = std::min(temp_heading_range.second, hypothesis_heading_end);

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

    ProcessStatus LogicBasedInitiator::process_nodes(std::vector<HypothesisNode> &node, std::vector<std::array<TrackPoint, 4>> &new_tracks)
    {
        if (hypothesis_layers_[0].size() + node.size() >= MAX_BINS * LOGIC_BASED_NUM_Y_BINS * LOGIC_BASED_MAX_NODE_PER_BINS ||
            node.size() >= LOGIC_BASED_MAX_NODE_PER_BINS)
        {
            // 假设节点总数超过上限，返回错误状态
            return ProcessStatus::TOO_MANY_HYPOTHSIS;
        }

        const size_t node_size = node.size();
        assert(node_size > 0 && "process_nodes函数输入的假设节点列表不能为空");
        if (node[0].depth == 0) // depth=0必然只有一个假设，且不需要进行距离和航向一致性评估，直接加入即可
        {
            // 新假设直接加入，不需要处理
            for (size_t index = 0; index < node_size; ++index)
            {
                hypothesis_layers_[0].push_back(node[index]);
                current_hypothesis_index_[node[index].bin_index].push_back(&hypothesis_layers_[0].back()); // 存放下索引
            }
            return ProcessStatus::SUCCESS;
        }

        // 获取时间差
        double dt = static_cast<double>(timestamp_batches_[0].milliseconds - timestamp_batches_[1].milliseconds) / 1000.0; // s

        // ===== 第一步：计算参考值、统计信息 =====
        std::vector<double> distance_ref_values(node_size); // 距离参考值数组
        std::vector<double> angle_ref_values(node_size);    // 角度参考值数组
        size_t depth_max = 0;                               // 提取最大深度
        double softmax_angle_sum = 0.0;                     // 两个参数softmax的分母，和一个参数的softmax分子
        double softmax_distance_sum = 0.0;                  // 两个参数softmax的分母，和一个参数的softmax分子
        for (size_t index = 0; index < node_size; ++index)
        {
            // 计算偏移量
            double angle_ref = evaluate_heading(node[index]);
            double distance_ref = evaluate_distance(node[index], dt);

            // 计算参考值
            angle_ref_values[index] = std::exp(-angle_ref);
            distance_ref_values[index] = std::exp(-distance_ref);
            softmax_angle_sum += angle_ref_values[index]; // 计算分母
            softmax_distance_sum += distance_ref_values[index];

            depth_max = std::max(depth_max, node[index].depth);
        }

        // ===== 第二步：计算置信度 =====
        std::vector<double> confidence_values(node_size);
        for (size_t index = 0; index < node_size; ++index)
        {
            // 计算原始匹配质量（距离+角度）
            double raw_match = LOGIC_BASED_HEADING_CONFIDENCE_WEIGHT * (angle_ref_values[index] / softmax_angle_sum) +
                               LOGIC_BASED_DISTANCE_CONFIDENCE_WEIGHT * (distance_ref_values[index] / softmax_distance_sum);

            double parent_conf = node[index].parent_node->confidence; // 父假设置信度

            // 按深度衰减（重复乘以原始匹配质量）
            double decayed_match = raw_match;
            for (size_t d = 0; d < (depth_max - node[index].depth); ++d)
            {
                decayed_match *= raw_match; // 每层衰减一次
            }

            // 最终置信度 = 衰减后的匹配质量 × 父置信度
            confidence_values[index] = decayed_match * parent_conf;
        }

        // ===== 第三步：K-BEST剪枝 =====
        std::array<size_t, LOGIC_BASED_K_BEST> best_indices;
        size_t k_num = k_best_selection(confidence_values, best_indices);

        // ===== 第四步：重新分配权重，分配假设节点归属 =====
        double confidence_sum = 0.0;
        for (size_t i = 0; i < k_num; ++i)
        {
            confidence_sum += confidence_values[best_indices[i]];
        }
        // 防止除零（若总和为0，则均匀分配）
        if (confidence_sum < 1e-12)
            confidence_sum = 1.0;
        size_t wait_to_track_index = SIZE_MAX;
        double hypothesis_to_track_best_confidence = 0.0;
        for (size_t i = 0; i < k_num; ++i)
        {
            size_t best_index = best_indices[i];
            node[best_index].confidence = confidence_values[best_index] / confidence_sum; // 归一化

            if (node[best_index].depth == 3)
            {
                if (node[best_index].confidence > hypothesis_to_track_best_confidence)
                {
                    hypothesis_to_track_best_confidence = node[best_index].confidence;
                    wait_to_track_index = best_index;
                }
                continue; // depth=3 节点不加入假设层
            }

            // 加入对应深度的假设层（注意：原代码用 hypothesis_layers_[0]，应改为 node[best_index].depth）
            hypothesis_layers_[node[best_index].depth].push_back(node[best_index]);
            current_hypothesis_index_[node[best_index].bin_index].push_back(
                &hypothesis_layers_[node[best_index].depth].back());
        }

        // 第五步：处理depth=3的节点，直接输出//TODO理论上应该置信度不够的时候再等一批，但是框架不允许了，面向论文CODE，暂时不动了
        if (wait_to_track_index != SIZE_MAX)
        {
            const auto &best_node = node[wait_to_track_index];

            // debug
            assert(best_node.associated_point != nullptr && best_node.parent_node != nullptr &&
                   best_node.parent_node->associated_point != nullptr && best_node.parent_node->parent_node != nullptr &&
                   best_node.parent_node->parent_node->associated_point != nullptr && best_node.parent_node->parent_node->parent_node != nullptr &&
                   best_node.parent_node->parent_node->parent_node->associated_point != nullptr &&
                   "depth=3的假设节点及其父节点必须都关联有效点迹");

            // debug
            LOG_DEBUG << "heading范围: [" << best_node.heading_start << ", " << best_node.heading_end << "], confidence: " << best_node.confidence;

            new_tracks.push_back({*(best_node.parent_node->parent_node->parent_node->associated_point),
                                  *(best_node.parent_node->parent_node->associated_point),
                                  *(best_node.parent_node->associated_point),
                                  *(best_node.associated_point)});
        }

        return ProcessStatus::SUCCESS;
    }

    // 计算与父节点外推点迹最小距离误差
    double LogicBasedInitiator::evaluate_distance(const HypothesisNode &node, double dt) const
    {
        // 获取当前节点的观测点迹，必须存在
        const TrackPoint *current_point = node.associated_point;
        if (!current_point)
        {
            // 如果点迹为空，返回一个极大值表示距离不一致
            return std::numeric_limits<double>::max();
        }

        // 获取父节点及其点迹（depth>0时保证存在）
        const HypothesisNode *parent = node.parent_node;
        if (!parent || !parent->associated_point)
        {
            return std::numeric_limits<double>::max();
        }
        const TrackPoint *parent_point = parent->associated_point;

        // 提取坐标（单位：km）
        double xn = current_point->x;
        double yn = current_point->y;
        double xf = parent_point->x;
        double yf = parent_point->y;

        // 径向速度转换：从 m/s 转换为 km/s
        double v_r_km_per_s = current_point->doppler / 1000.0;

        // 计算当前点距离原点的距离 sqrt(xn^2 + yn^2)
        double rn = std::sqrt(xn * xn + yn * yn);
        if (rn < 1e-9)
        {
            // 如果当前点在原点附近，无法计算有效距离，返回极大值
            return std::numeric_limits<double>::max();
        }

        // 计算 C 值（所有量统一为 km 和 km/s）
        double C = -v_r_km_per_s * dt * rn + (xn * xn + yn * yn);

        // 计算分子绝对值
        double numerator = std::fabs(xn * xf + yn * yf + C);

        // 最终距离误差 d_i
        double d_i = numerator / rn;

        return d_i;
    }

    // 计算与父节点外推点迹航向一致性，输入当前节点和时间差，获取当前节点和父节点的航向范围，计算当前节点航向中心与父节点航向范围的关系，变化量越小越好
    double LogicBasedInitiator::evaluate_heading(const HypothesisNode &node) const
    {
        // 获取当前节点的航向范围
        double h_ns = node.heading_start;
        double h_ne = node.heading_end;

        // 获取父节点（必须存在）
        const HypothesisNode *parent = node.parent_node;
        if (!parent)
        {
            // 若无父节点，返回一个极大值表示航向不一致
            return std::numeric_limits<double>::max();
        }
        double h_fs = parent->heading_start;
        double h_fe = parent->heading_end;

        // 计算当前节点航向中心 h_nc
        double h_nc = (h_ns + h_ne) / 2.0;

        // 定义角度差规范化函数：返回 [0, π] 内的最小弧度差
        auto angle_diff = [](double a, double b) -> double
        {
            double diff = std::fabs(a - b);
            // 处理角度周期性（假设输入在 [0, 2π)）
            diff = std::fmod(diff, 2.0 * M_PI);
            if (diff > M_PI)
            {
                diff = 2.0 * M_PI - diff;
            }
            return diff;
        };

        // 计算三项
        double term1 = angle_diff(h_nc, h_fs);
        double term2 = angle_diff(h_nc, h_fe);
        double term3 = angle_diff(h_fs, h_fe); // 父节点范围宽度

        // 最终 h_out
        double h_out = term1 + term2 - term3;

        return h_out;
    }

    size_t LogicBasedInitiator::k_best_selection(const std::vector<double> &confidence, std::array<size_t, LOGIC_BASED_K_BEST> &best_indices) const
    {
        const size_t N = confidence.size();
        const size_t K = LOGIC_BASED_K_BEST;
        const size_t selected = std::min(N, K); // 实际选中的数量

        // 创建索引数组 [0, 1, ..., N-1]
        std::vector<size_t> idx(N);
        for (size_t i = 0; i < N; ++i)
        {
            idx[i] = i;
        }

        // 部分排序：只排序前 selected 个最大的元素（降序）
        std::partial_sort(idx.begin(), idx.begin() + (long int)selected, idx.end(),
                          [&confidence](size_t a, size_t b)
                          {
                              return confidence[a] > confidence[b]; // 按置信度从大到小
                          });

        // 将前 selected 个索引复制到输出数组
        for (size_t i = 0; i < selected; ++i)
        {
            best_indices[i] = idx[i];
        }

        // 将剩余位置标记为无效（可选，便于调试）
        for (size_t i = selected; i < K; ++i)
        {
            best_indices[i] = SIZE_MAX;
        }

        return selected;
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