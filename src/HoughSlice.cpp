#include <cassert> //静态断言检查避免忘了用回调函数
#include <algorithm>
#include <cmath>

#include "HoughSlice.hpp"
#include "Func_dbscan.hpp"
#include "../include/defsystem.h"
#include "../utils/Logger.hpp"

namespace track_project::trackinit
{
    ProcessStatus SliceHough::process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        // 清空输出航迹
        new_track.clear();

        // STEP1：对于新来的点迹生成聚类
        process_cluster_generation(points);
        // STEP1 DEBUG宏，用于查看当前聚类信息
        LOG_DEBUG << "=== 聚类信息可视化 ===";
        LOG_DEBUG << "当前聚类数量：" << ClustArea.get_allocated_count(); // ⚠仅建议DEBUG下使用这个函数，因为真的很浪费时间
        for (auto &clust : ClustArea.get_allocated_ptrs())
        {
            LOG_DEBUG << "当前聚类中心：(" << clust->center_x << "," << clust->center_y << ")，共" << clust->current_batch_index << "批次："
                      << "第一批次点迹数量：" << clust->point_list[0].size()
                      << ", 第二批次点迹数量：" << clust->point_list[1].size()
                      << ", 第三批次点迹数量：" << clust->point_list[2].size()
                      << ", 第四批次点迹数量：" << clust->point_list[3].size();
        }

        // STEP2：对于当前聚类中的所有点迹进行霍夫变换投票
        LOG_DEBUG << "=== 霍夫投票过程可视化 ===";
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index != 3)
            {
                LOG_DEBUG << "聚类(" << it_clust->center_x << "," << it_clust->center_y << "): 批次数不足，跳过投票";
                continue; // 批次数不足，跳过
            }

            for (size_t batch = 0; batch <= 3; ++batch)
            {
                for (const auto &point : it_clust->point_list[batch]) // 对于每个点执行官投票
                {
                    process_point_for_hough_vote(batch, point, *it_clust);
                }
            }

            // STEP2 DEBUG宏，用于查看霍夫变换空间投票结果
#ifndef NDEBUG
            {
                // 生成文件名
                std::stringstream filename;
                filename << "/home/holyogre/TrackInit/hough_debug_"
                         << std::fixed << std::setprecision(1) << it_clust->center_x << "_"
                         << std::setprecision(1) << it_clust->center_y << ".dat";

                // 打开文件
                std::ofstream file(filename.str(), std::ios::binary);
                if (file.is_open())
                {
                    // 写入维度
                    uint32_t dims[2] = {static_cast<uint32_t>(HOUGH_RHO_DIM),
                                        static_cast<uint32_t>(HOUGH_THETA_DIM)};
                    file.write(reinterpret_cast<const char *>(dims), 2 * sizeof(uint32_t));

                    // 统计信息

                    // 写入数据并统计
                    for (size_t theta = 0; theta < HOUGH_THETA_DIM; ++theta)
                    {
                        for (size_t rho = 0; rho < HOUGH_RHO_DIM; ++rho)
                        {
                            BitArray<256> votes = it_clust->vote_area[theta][rho];
                            file.write(reinterpret_cast<const char *>(&votes), sizeof(BitArray<256>));
                        }
                    }

                    file.close();
                }
            }
#endif
        }

        // 峰值检测和航迹生成
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index != 3)
            {
                LOG_DEBUG << "聚类(" << it_clust->center_x << "," << it_clust->center_y << "): 批次数不足，跳过峰值检测";
                continue; // 批次数不足，跳过
            }

            // STEP3:峰值检测
            std::vector<std::array<double, 3>> detected_lines = process_extract_peak_from_hough_space(*it_clust);
            LOG_DEBUG << "聚类中心(" << it_clust->center_x << "," << it_clust->center_y << ")检测到直线数量：" << detected_lines.size();

            // STEP4:回溯点迹，生成航迹
            process_backtrack_points(detected_lines, *it_clust, new_track);
        }

        // 调用回调函数发送航迹
        assert(trackCallback_ != nullptr && "HoughSlice类调用的时候没有绑定回调函数"); // debug模式下检查
        trackCallback_(new_track);

        // 删除已经执行过航迹生成的霍夫变换切片
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index == 3)
            {
                ClustArea.release_target(it_clust);
            }
        }

        return ProcessStatus::SUCCESS;
    }

    // 先剔除在聚类中心的点，再剔除离群点，最后对剩余点进行聚类
    void SliceHough::process_cluster_generation(const std::vector<TrackPoint> &points)
    {
        // 用于记录当前点迹是否被聚类
        std::vector<bool> visited(points.size(), false);

        // 循环访问聚类，看能不能存放到之前的聚类里面去
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (const auto &it_clust : clustAera_list)
        {
            it_clust->current_batch_index++; // 批次索引加1

            //  检测点迹是否在中心点附近
            for (auto it_point = points.begin(); it_point != points.end(); ++it_point)
            {
                const double dx = it_point->x - it_clust->center_x;
                const double dy = it_point->y - it_clust->center_y;
                const double dist_square = dx * dx + dy * dy;
                if (dist_square <= SLICEHOUGH_CLUSTER_RADIUS_KM * SLICEHOUGH_CLUSTER_RADIUS_KM)
                {
                    // 点迹在聚类范围内，加入该聚类
                    size_t batch = it_clust->current_batch_index; // 当前批次数据存放位置
                    if (batch >= it_clust->point_list.size())
                    {
                        continue; // 超出批次数量，忽略
                    }
                    it_clust->point_list[batch].push_back(*it_point);

                    // 标记该点迹已被聚类
                    size_t point_idx = std::distance(points.begin(), it_point);
                    visited[point_idx] = true;
                }
            }
        }

        // 将未被聚类的点迹全部存放到新的向量里面
        std::vector<TrackPoint> unvisited_points;
        for (size_t i = 0; i < points.size(); ++i)
        {
            if (!visited[i])
            {
                unvisited_points.push_back(points[i]);
            }
        }

        // 新聚类
        std::vector<std::vector<size_t>> _clust = dbscan(unvisited_points, SLICEHOUGH_CORE_POINT_RADIUS_KM, 2);

        // 计算相关参数，申请新的聚类空间
        for (const auto &cluster : _clust)
        {
            if (cluster.empty()) // 忽略空集
            {
                continue;
            }

            // 计算聚类中心
            double sum_x = 0.0;
            double sum_y = 0.0;
            for (size_t idx : cluster)
            {
                sum_x += unvisited_points[idx].x;
                sum_y += unvisited_points[idx].y;
            }
            double center_x = sum_x / cluster.size();
            double center_y = sum_y / cluster.size();

            // 申请新的聚类空间
            Slice *new_slice = ClustArea.acquire_target();

            new_slice->center_x = center_x;
            new_slice->center_y = center_y;
            new_slice->current_batch_index = 0;

            // 将点迹加入新的聚类
            for (size_t idx : cluster)
            {
                new_slice->point_list[0].push_back(unvisited_points[idx]);
            }
        }
    }

    // 对于单个点迹进行处理，在霍夫空间中进行投票,减少计算量，内部使用极坐标
    void SliceHough::process_point_for_hough_vote(size_t batch, const TrackPoint &point, Slice &it_clust)
    {
        // 雷达站位于0,0位置
        double rel_x = point.x - it_clust.center_x;
        double rel_y = point.y - it_clust.center_y;

        // 依据doppler和velocity_max计算所有可能的航向角度序列
        double speed = std::abs(point.doppler); // 取绝对值，单位m/s
        if (speed > track_project::velocity_max)
        {
            return; // 速度超过最大值，跳过该点迹
        }

        // 计算基准航向，依据doppler的正负确定是朝向还是背向，为极坐标表示方法，便于调用math库函数计算
        double line_of_sight_rad = std::atan2(point.y, point.x); // 观测航向
        double base_dir_rad = 0.0;
        if (point.doppler > 0) // 表示靠近还是原理，以此决定是否加PI
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

        // 计算可能的速度方向与视线方向的夹角
        double angle_rad = std::acos(speed / track_project::velocity_max); // 计算夹角，弧度，范围[0,pi/2]

        double heading1 = base_dir_rad - angle_rad;
        double heading2 = base_dir_rad + angle_rad;

        //  区间为(heading1,heading2)∪(heading3,heading4)，且heading1<=heading2<=heading3<=heading4
        if (heading1 < 0) // 仅有可能heading1小于0
        {
            double heading3 = 0.0, heading4 = M_PI;
            // 射线方向性可以由doppler确定，故可以将区间调整为(0,heading2)∪(heading3+π,π)
            heading3 = heading1 + M_PI;
            heading4 = M_PI;
            heading1 = 0.0;
            heading2 = heading2;
            vote_in_hough_space(heading1, heading2, rel_x, rel_y, batch, point.doppler, it_clust.vote_area);
            vote_in_hough_space(heading3, heading4, rel_x, rel_y, batch, point.doppler, it_clust.vote_area);
        }
        else if (heading2 > 2 * M_PI) // 仅有可能heading2大于2π
        {
            double heading3 = 0.0, heading4 = M_PI;
            // 射线方向性可以由doppler确定，由angle_rad保证heading3必然位于(π,2π)，故区间为(0,heading2-2π)∪(heading1-π,π)
            heading3 = heading1 - M_PI;
            heading4 = M_PI;
            heading1 = 0.0;
            heading2 = heading2 - 2 * M_PI;
            vote_in_hough_space(heading1, heading2, point.doppler, rel_x, rel_y, batch, it_clust.vote_area);
            vote_in_hough_space(heading3, heading4, point.doppler, rel_x, rel_y, batch, it_clust.vote_area);
        }
        else
        {
            vote_in_hough_space(heading1, heading2, point.doppler, rel_x, rel_y, batch, it_clust.vote_area);
        }
    }

    // 对于霍夫变换空间中的指定区间进行特殊投票
    void SliceHough::vote_in_hough_space(const double heading_start, const double heading_end, const double doppler,
                                         const double rel_x, const double rel_y, const size_t batch,
                                         std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area)
    {
        // 计算起始和结束的角度索引
        std::uint32_t angle_idx_start = static_cast<std::uint32_t>(heading_start * 180.0 / M_PI / SLICEHOUGH_THETA_RESOLUTION_DEG);
        std::uint32_t angle_idx_end = static_cast<std::uint32_t>(heading_end * 180.0 / M_PI / SLICEHOUGH_THETA_RESOLUTION_DEG);

        // 计算速度索引 - 由宏控制位数
        double ratio = doppler / track_project::velocity_max;
        ratio = std::clamp(ratio, -1.0, 1.0);

        // 根据SLICEHOUGH_DOPPLER_BIT_NUM计算每侧级别数
        constexpr size_t LEVELS_PER_SIDE = SLICEHOUGH_DOPPLER_BIT_NUM / 2; // 正负各一半
        constexpr size_t BYTES_PER_BATCH = SLICEHOUGH_DOPPLER_BIT_NUM / 8; // 每批次的字节数

        // 创建掩码缓冲区（按字节存储，小端序）
        std::vector<uint8_t> mask_bytes(BYTES_PER_BATCH, 0);

        // 速度位的存储布局（连续均匀分布）：
        // 低半区 (bits 0 到 SLICEHOUGH_DOPPLER_BIT_NUM/2 - 1)：存储负速度（从最大负到最小负）
        //    [bit0: 最大负速度 (-1.0), bit1: 次大负速度 (-0.9), ..., bit31: 最小负速度 (-0.1)]
        // 高半区 (bits SLICEHOUGH_DOPPLER_BIT_NUM/2 到 SLICEHOUGH_DOPPLER_BIT_NUM-1)：存储正速度（从最小正到最大正）
        //    [bit32: 最小正速度 (+0.1), bit33: 次小正速度 (+0.2), ..., bit63: 最大正速度 (+1.0)]
        //
        // 这样整个速度范围是连续的：-1.0, -0.9, ..., -0.1, +0.1, +0.2, ..., +1.0
        // ratio == 0 的情况：不投票（不可能存在这种点，除非出BUG了）
        if (ratio == 0)
        {
            return;
        }
        else if (ratio < 0)
        {
            // 负数：使用低半区 (bits 0 到 SLICEHOUGH_DOPPLER_BIT_NUM/2 - 1)
            double abs_ratio = -ratio; // 0.0 ~ 1.0
            int level = static_cast<int>(std::floor(abs_ratio * LEVELS_PER_SIDE));
            level = std::clamp(level, 0, static_cast<int>(LEVELS_PER_SIDE) - 1);

            // 负速度：从最大负（bit0）到最小负（bit31）
            // level=0 (|速度|最大, -1.0) → bit_pos=0
            // level=31 (|速度|最小, -0.1) → bit_pos=31
            int bit_pos = level; // 注意：这里直接用level，因为负速度在低半区

            // 设置对应位
            size_t byte_idx = bit_pos / 8;
            size_t bit_in_byte = bit_pos % 8;
            mask_bytes[byte_idx] |= (1 << bit_in_byte);
        }
        else if (ratio > 0)
        {
            // 正数：使用高半区 (bits SLICEHOUGH_DOPPLER_BIT_NUM/2 到 SLICEHOUGH_DOPPLER_BIT_NUM-1)
            int level = static_cast<int>(std::floor(ratio * LEVELS_PER_SIDE));
            level = std::clamp(level, 0, static_cast<int>(LEVELS_PER_SIDE) - 1);

            // 正速度：从最小正（bit32）到最大正（bit63）
            // level=0 (+0.1) → bit_pos=32
            // level=31 (+1.0) → bit_pos=63
            int bit_pos = SLICEHOUGH_DOPPLER_BIT_NUM / 2 + level;

            // 设置对应位
            size_t byte_idx = bit_pos / 8;
            size_t bit_in_byte = bit_pos % 8;
            mask_bytes[byte_idx] |= (1 << bit_in_byte);
        }

        // 遍历所有角度索引进行投票
        for (std::uint32_t angle_idx = angle_idx_start; angle_idx < angle_idx_end; ++angle_idx)
        {
            double theta = angle_idx * SLICEHOUGH_THETA_RESOLUTION_DEG * M_PI / 180.0;
            double distance = rel_x * std::cos(theta) + rel_y * std::sin(theta);

            // 计算距离索引
            int distance_idx = static_cast<int>((distance + 2 * SLICEHOUGH_CLUSTER_RADIUS_KM) / SLICEHOUGH_RHO_RESOLUTION_KM);
            if (distance_idx < 0 || distance_idx >= static_cast<int>(HOUGH_RHO_DIM))
            {
                continue;
            }

            // 获取对应的投票单元
            BitArray<HOUGH_VOTE_BIT_NUM> &cell = (angle_idx >= HOUGH_THETA_DIM)
                                                     ? vote_area[angle_idx - HOUGH_THETA_DIM][distance_idx]
                                                     : vote_area[angle_idx][distance_idx];

            // 计算这个批次在总位数中的起始字节偏移
            size_t byte_offset = batch * BYTES_PER_BATCH;

            // 使用or_bytes进行投票
            cell.or_bytes(byte_offset, mask_bytes.data(), BYTES_PER_BATCH);
        }
    }

    // 峰值检测,从霍夫变换空间中提取检测到的直线参数列表
    std::vector<std::array<double, 3>> SliceHough::process_extract_peak_from_hough_space(Slice &cluster) const
    {
        std::vector<std::array<double, 3>> detected_lines; // 存储检测到的直线参数（theta, rho, doppler）

        for (size_t angle_idx = 0; angle_idx < HOUGH_THETA_DIM; ++angle_idx)
        {
            for (size_t distance_idx = 0; distance_idx < HOUGH_RHO_DIM; ++distance_idx)
            {
                // 提取4个批次的共同投票位（返回SLICEHOUGH_DOPPLER_BIT_NUM位的BitArray）
                auto common_votes = peak_filter(angle_idx, distance_idx, cluster.vote_area);

                // 如果没有共同投票，跳过
                bool has_votes = false;
                for (size_t i = 0; i < SLICEHOUGH_DOPPLER_BIT_NUM; ++i)
                {
                    if (common_votes.get_bit(i))
                    {
                        has_votes = true;
                        break;
                    }
                }
                if (!has_votes)
                    continue;

                // 遍历所有速度位
                for (size_t speed_bit = 0; speed_bit < SLICEHOUGH_DOPPLER_BIT_NUM; ++speed_bit)
                {
                    if (common_votes.get_bit(speed_bit))
                    {
                        double ratio = 0.0;

                        // 根据布局计算ratio：
                        // 低半区 (bits 0 到 SLICEHOUGH_DOPPLER_BIT_NUM/2 - 1)：负速度，从最大负(-1.0)到最小负
                        // 高半区 (bits SLICEHOUGH_DOPPLER_BIT_NUM/2 到 SLICEHOUGH_DOPPLER_BIT_NUM-1)：正速度，从最小正到最大正(+1.0)

                        constexpr size_t HALF_BITS = SLICEHOUGH_DOPPLER_BIT_NUM / 2;
                        constexpr size_t MAX_LEVEL = HALF_BITS - 1;
                        constexpr double SPEED_STEP = 1.0 / HALF_BITS; // 每级的速度增量

                        if (speed_bit < HALF_BITS) // 负速度区域
                        {
                            // speed_bit = 0        → -1.0
                            // speed_bit = MAX_LEVEL → -SPEED_STEP
                            // 公式：从 -1.0 逐渐增加到接近0
                            ratio = -1.0 + speed_bit * SPEED_STEP;
                        }
                        else // 正速度区域
                        {
                            // speed_bit = HALF_BITS       → +SPEED_STEP
                            // speed_bit = SLICEHOUGH_DOPPLER_BIT_NUM-1 → +1.0
                            size_t level = speed_bit - HALF_BITS;
                            ratio = (level + 1) * SPEED_STEP; // +1 是为了跳过0
                        }

                        double doppler = ratio * track_project::velocity_max;

                        // 计算theta和rho值
                        double theta = angle_idx * SLICEHOUGH_THETA_RESOLUTION_DEG * M_PI / 180.0; // 弧度
                        double rho = distance_idx * SLICEHOUGH_RHO_RESOLUTION_KM - 2 * SLICEHOUGH_CLUSTER_RADIUS_KM;

                        // 保存检测到的直线参数
                        detected_lines.push_back({theta, rho, doppler});
                    }
                }
            }
        }

        return detected_lines;
    }

    // 提取目标区域峰值
    BitArray<SLICEHOUGH_DOPPLER_BIT_NUM> SliceHough::peak_filter(
        const size_t angle_idx, const size_t distance_idx,
        const std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area) const
    {
        const auto &cell = vote_area[angle_idx][distance_idx];

        constexpr size_t BYTES_PER_BATCH = SLICEHOUGH_DOPPLER_BIT_NUM / 8;
        constexpr size_t BATCH_NUM = HOUGH_VOTE_BIT_NUM / SLICEHOUGH_DOPPLER_BIT_NUM;

        // 分配两个缓冲区
        std::vector<uint8_t> buf1(BYTES_PER_BATCH);
        std::vector<uint8_t> buf2(BYTES_PER_BATCH);

        // 读取第一个批次
        cell.read_bytes(buf1.data(), 0, BYTES_PER_BATCH);

        // 创建结果对象并写入
        BitArray<SLICEHOUGH_DOPPLER_BIT_NUM> result;
        result.write_bytes(buf1.data(), 0, BYTES_PER_BATCH);

        // 依次处理后续批次
        for (size_t batch = 1; batch < BATCH_NUM; ++batch)
        {
            cell.read_bytes(buf2.data(), batch * BYTES_PER_BATCH, BYTES_PER_BATCH);

            BitArray<SLICEHOUGH_DOPPLER_BIT_NUM> current;
            current.write_bytes(buf2.data(), 0, BYTES_PER_BATCH);

            result &= current;
        }

        return result;
    }

    // 回溯点迹
    void SliceHough::process_backtrack_points(const std::vector<std::array<double, 3>> &detected_lines, const Slice &cluster,
                                              std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        const double RHO_TOL = SLICEHOUGH_RHO_RESOLUTION_KM / 2.0;
        const double DOPPLER_TOL = 1.0;
        const double CENTER_X = cluster.center_x;
        const double CENTER_Y = cluster.center_y;

        for (const auto &line : detected_lines)
        {
            double theta = line[0];
            double rho = line[1];
            double doppler = line[2];

            double cos_theta = std::cos(theta);
            double sin_theta = std::sin(theta);

            std::array<TrackPoint, 4> track;
            bool valid = true;

            // 检查4个批次
            for (size_t batch = 0; batch < 4 && valid; ++batch)
            {
                const auto &points = cluster.point_list[batch];
                bool found = false;

                for (const auto &point : points)
                {
                    double rel_x = point.x - CENTER_X;
                    double rel_y = point.y - CENTER_Y;
                    double point_rho = rel_x * cos_theta + rel_y * sin_theta;

                    if (std::abs(point_rho - rho) < RHO_TOL &&
                        std::abs(point.doppler - doppler) < DOPPLER_TOL)
                    {
                        track[batch] = point;
                        found = true;
                        break;
                    }
                }

                if (!found)
                    valid = false;
            }

            if (valid)
            {
                new_track.push_back(track);
            }
        }
    }

    void SliceHough::clear_all()
    {
        ClustArea.clear_all();
    }
} // namespace track_project::trackinit
