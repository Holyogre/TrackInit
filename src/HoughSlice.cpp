#include <cassert> //静态断言检查避免忘了用回调函数
#include <algorithm>
#include <cmath>

#include "HoughSlice.hpp"
#include "Func_dbscan.hpp"
#include "../include/defsystem.h"
#include "../utils/Logger.hpp"

namespace track_project::trackinit
{
    // 主程序
    ProcessStatus HoughSlice::process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        // 清空输出航迹
        new_track.clear();

        // ==================== STEP1: 聚类(HOUGHSLICE)生成 ====================
        process_cluster_generation(points);
// STEP1 DEBUG宏，用于查看当前聚类信息
#ifndef NDEBUG
        LOG_DEBUG << "=== [STEP1] 聚类信息可视化 ===";
        LOG_DEBUG << "当前聚类数量：" << ClustArea.get_allocated_count(); // ⚠仅建议DEBUG下使用这个函数，因为真的很浪费时间
        for (auto &clust : ClustArea.get_allocated_ptrs())
        {
            LOG_DEBUG << "当前聚类中心：(" << clust->center_x << "," << clust->center_y << ")，共" << clust->current_batch_index << "批次：";
            for (int batch = 0; batch < HOUGHSLICE_BATCH_NUM; ++batch)
            {
                LOG_DEBUG << "第" << batch + 1 << "批次点迹数量：" << clust->point_list[batch].size();
            }
        }
#endif

        // ==================== STEP2: 霍夫变换投票 ====================
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index != HOUGHSLICE_BATCH_NUM - 1)
            {
                LOG_DEBUG << "聚类(" << it_clust->center_x << "," << it_clust->center_y << "): 批次数不足，跳过投票";
                continue; // 批次数不足，跳过
            }

            for (size_t batch = 0; batch < HOUGHSLICE_BATCH_NUM; ++batch)
            {
                for (const auto &point : it_clust->point_list[batch]) // 对于每个点执行霍夫变换投票
                {
                    process_point_for_hough_vote(batch, point, *it_clust);
                }
            }
        }

// STEP2 DEBUG宏，用于查看霍夫变换空间投票结果
#ifndef NDEBUG
        LOG_DEBUG << "=== [STEP2] 霍夫投票过程可视化 ===";
        for (auto &it_clust : clustAera_list)
        {
            // 生成文件名
            std::stringstream filename;
            filename << "../hough_debug_"
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

                // 写入数据并统计
                for (size_t theta = 0; theta < HOUGH_THETA_DIM; ++theta)
                {
                    for (size_t rho = 0; rho < HOUGH_RHO_DIM; ++rho)
                    {
                        BitArray<HOUGHSLICE_DOPPLER_BIT_NUM * HOUGHSLICE_BATCH_NUM> votes = it_clust->vote_area[theta][rho];
                        file.write(reinterpret_cast<const char *>(&votes), HOUGHSLICE_DOPPLER_BIT_NUM * HOUGHSLICE_BATCH_NUM / 8);
                    }
                }

                file.close();
            }
            LOG_DEBUG << "结果存储在" << filename.str();
        }
#endif

        // ========================================STEP3~5:峰值检测、点迹凝聚、航迹生成========================================
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index != HOUGHSLICE_BATCH_NUM - 1)
            {
                LOG_DEBUG << "聚类(" << it_clust->center_x << "," << it_clust->center_y << "): 批次数不足，跳过峰值检测";
                continue; // 批次数不足，跳过
            }

            // ==================== STEP3: 峰值检测 ====================
            std::vector<std::vector<std::array<size_t, 2>>> detected_lines = process_extract_peak_from_hough_space(*it_clust);

#ifndef NDEBUG
            LOG_DEBUG << "===== [STEP3] 峰值检测 =====";
            LOG_DEBUG << "聚类中心(" << it_clust->center_x << "," << it_clust->center_y << ")";
            size_t total_peaks = 0;
            std::vector<size_t> doppler_counts(detected_lines.size(), 0);
            for (size_t d = 0; d < detected_lines.size(); ++d)
            {
                doppler_counts[d] = detected_lines[d].size();
                total_peaks += doppler_counts[d];
            }
            LOG_DEBUG << "检测到峰值总数：" << total_peaks;

            // ==================== 保存峰值点迹到CSV文件 ====================
            // 生成文件名，包含时间戳或聚类中心信息
            char filename[256];
            snprintf(filename, sizeof(filename), "../peaks_cluster_%.0f_%.0f_%ld.csv",
                     it_clust->center_x, it_clust->center_y, time(nullptr));

            FILE *fp = fopen(filename, "w");
            if (fp)
            {
                // 写入CSV表头（可选，如果不需要可以注释掉）
                fprintf(fp, "THETA,RHO,DOPPLER\n");

                size_t written_count = 0;
                for (size_t d = 0; d < detected_lines.size(); ++d)
                {
                    for (const auto &peak : detected_lines[d])
                    {
                        // peak[0] = theta索引, peak[1] = rho索引
                        fprintf(fp, "%zu,%zu,%zu\n", peak[0], peak[1], d);
                        written_count++;
                    }
                }

                fclose(fp);
                LOG_DEBUG << "已保存 " << written_count << " 个峰值点到CSV文件: " << filename;
            }
            else
            {
                LOG_ERROR << "无法创建CSV文件: " << filename;
            }

#endif

            // ==================== STEP4: 点迹凝聚 ====================
            std::vector<std::array<double, 3>> condensed_lines = process_condense_detected_lines(detected_lines);
#ifndef NDEBUG
            LOG_DEBUG << "===== [STEP4] 点迹凝聚 =====";
            LOG_DEBUG << "凝聚后直线数量：" << condensed_lines.size();
            // 输出每条凝聚直线的详细参数
            for (size_t i = 0; i < condensed_lines.size(); ++i)
            {
                const auto &line = condensed_lines[i];
                LOG_DEBUG << "直线[" << i << "]: theta=" << line[0]
                          << " rad (" << line[0] * 180.0 / M_PI << "°)"
                          << ", rho=" << line[1] << " km"
                          << ", doppler=" << line[2] << " m/s";
            }

#endif

            // ==================== STEP5: 点迹回溯 ====================
            // STEP5不需要DEBUG日志，可以从显控查看结果
            process_backtrack_points(condensed_lines, *it_clust, new_track);
        }

        // 调用回调函数发送航迹
        assert(trackCallback_ != nullptr && "HoughSlice类调用的时候没有绑定回调函数"); // debug模式下检查
        trackCallback_(new_track);

        // 删除已经执行过航迹生成的霍夫变换切片
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index == HOUGHSLICE_BATCH_NUM - 1)
            {
                ClustArea.release_target(it_clust);
            }
        }

        return ProcessStatus::SUCCESS;
    }

    // 先剔除在聚类中心的点，再剔除离群点，最后对剩余点进行聚类
    void HoughSlice::process_cluster_generation(const std::vector<TrackPoint> &points)
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
                if (dist_square <= HOUGHSLICE_CLUSTER_RADIUS_KM * HOUGHSLICE_CLUSTER_RADIUS_KM)
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
        std::vector<std::vector<size_t>> _clust = dbscan(unvisited_points, HOUGHSLICE_CORE_POINT_RADIUS_KM, 2);

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
    void HoughSlice::process_point_for_hough_vote(size_t batch, const TrackPoint &point, Slice &it_clust)
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

        // ==================== 计算 doppler_tolerance_bits ====================
        // 计算时间间隔 dt（秒）：使用 point_list 的第一个和最后一个批次的时间戳
        int doppler_tolerance_bits = 0;
        if (!it_clust.point_list[0].empty() && !it_clust.point_list[HOUGHSLICE_BATCH_NUM - 1].empty())
        {
            int64_t t_start = it_clust.point_list[0][0].time.milliseconds;
            int64_t t_end = it_clust.point_list[HOUGHSLICE_BATCH_NUM - 1][0].time.milliseconds;
            double dt = static_cast<double>(t_end - t_start) / 1000.0 / (HOUGHSLICE_BATCH_NUM - 1); // 转换为秒

            if (dt > 0)
            {
                // 计算点迹到雷达的距离（km转m）
                double distance_m = std::sqrt(point.x * point.x + point.y * point.y) * 1000.0;

                // 计算航迹在 (BATCH_NUM-1) 个时间间隔内可能移动的最大距离
                double max_displacement = track_project::velocity_max * dt * (HOUGHSLICE_BATCH_NUM - 1);

                // 计算视线角变化 α = arctan(max_displacement / distance)
                double alpha = std::atan2(max_displacement, distance_m);

                // 计算多普勒速度的极限变化量 Δv = max_velocity * (1 - cos(α))
                double delta_v = track_project::velocity_max * (1.0 - std::cos(alpha));

                // 转换为位数：doppler_tolerance_bits = round((Δv / velocity_max) * DOPPLER_BIT_NUM)
                doppler_tolerance_bits = static_cast<int>(std::round((delta_v / track_project::velocity_max) * HOUGHSLICE_DOPPLER_BIT_NUM));
                doppler_tolerance_bits = std::clamp(doppler_tolerance_bits, 0, static_cast<int>(HOUGHSLICE_DOPPLER_BIT_NUM) - 1);
            }
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
            vote_in_hough_space(heading1, heading2, point.doppler, rel_x, rel_y, batch, doppler_tolerance_bits, it_clust.vote_area);
            vote_in_hough_space(heading3, heading4, point.doppler, rel_x, rel_y, batch, doppler_tolerance_bits, it_clust.vote_area);
        }
        else if (heading2 > 2 * M_PI) // 仅有可能heading2大于2π
        {
            double heading3 = 0.0, heading4 = M_PI;
            // 射线方向性可以由doppler确定，由angle_rad保证heading3必然位于(π,2π)，故区间为(0,heading2-2π)∪(heading1-π,π)
            heading3 = heading1 - M_PI;
            heading4 = M_PI;
            heading1 = 0.0;
            heading2 = heading2 - 2 * M_PI;
            vote_in_hough_space(heading1, heading2, point.doppler, rel_x, rel_y, batch, doppler_tolerance_bits, it_clust.vote_area);
            vote_in_hough_space(heading3, heading4, point.doppler, rel_x, rel_y, batch, doppler_tolerance_bits, it_clust.vote_area);
        }
        else
        {
            vote_in_hough_space(heading1, heading2, point.doppler, rel_x, rel_y, batch, doppler_tolerance_bits, it_clust.vote_area);
        }
    }

    // 对于霍夫变换空间中的指定区间进行特殊投票
    void HoughSlice::vote_in_hough_space(const double heading_start, const double heading_end, const double doppler,
                                         const double rel_x, const double rel_y, const size_t batch,
                                         const int doppler_tolerance_bits,
                                         std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area)
    {
        // 计算起始和结束的角度索引
        std::uint32_t angle_idx_start = static_cast<std::uint32_t>(heading_start * 180.0 / M_PI / HOUGHSLICE_THETA_RESOLUTION_DEG);
        std::uint32_t angle_idx_end = static_cast<std::uint32_t>(heading_end * 180.0 / M_PI / HOUGHSLICE_THETA_RESOLUTION_DEG);

        // 计算速度索引 - 由宏控制位数
        double ratio = doppler / track_project::velocity_max;
        ratio = std::clamp(ratio, -1.0, 1.0);

        // 根据HOUGHSLICE_DOPPLER_BIT_NUM计算每侧级别数
        constexpr size_t LEVELS_PER_SIDE = HOUGHSLICE_DOPPLER_BIT_NUM / 2; // 正负各一半
        constexpr size_t BYTES_PER_BATCH = HOUGHSLICE_DOPPLER_BIT_NUM / 8; // 每批次的字节数

        // 创建掩码缓冲区（按字节存储，小端序）
        std::vector<uint8_t> mask_bytes(BYTES_PER_BATCH, 0);

        // 构造一个批次的速度位的存储布局（连续均匀分布）：
        // 低半区 (bits 0 到 HOUGHSLICE_DOPPLER_BIT_NUM/2 - 1)：存储负速度（从最大负到最小负）
        //    [bit0: 最大负速度 (-1.0), bit1: 次大负速度 (-0.9), ..., bit31: 最小负速度 (-0.1)]
        // 高半区 (bits HOUGHSLICE_DOPPLER_BIT_NUM/2 到 HOUGHSLICE_DOPPLER_BIT_NUM-1)：存储正速度（从最小正到最大正）
        //    [bit32: 最小正速度 (+0.1), bit33: 次小正速度 (+0.2), ..., bit63: 最大正速度 (+1.0)]
        //
        // 这样整个速度范围是连续的：-1.0, -0.9, ..., -0.1, +0.1, +0.2, ..., +1.0
        // ratio == 0 的情况：不投票（不可能存在这种点，除非出BUG了）
        if (ratio == 0)
        {
            return;
        }

        // 计算中心位位置
        int center_bit_pos = 0;
        if (ratio < 0)
        {
            double abs_ratio = -ratio; // 0.0 ~ 1.0
            // 希望 abs_ratio=0.1 → bit=31, abs_ratio=1.0 → bit=0
            int level = static_cast<int>(std::floor((1.0 - abs_ratio) * LEVELS_PER_SIDE));
            level = std::clamp(level, 0, static_cast<int>(LEVELS_PER_SIDE) - 1);
            center_bit_pos = level; // 现在 level=31 for 0.1, level=0 for 1.0
        }
        else // ratio > 0
        {
            // 希望 ratio=0.1 → bit=32, ratio=1.0 → bit=63
            int level = static_cast<int>(std::floor((1.0 - ratio) * LEVELS_PER_SIDE));
            level = std::clamp(level, 0, static_cast<int>(LEVELS_PER_SIDE) - 1);
            center_bit_pos = HOUGHSLICE_DOPPLER_BIT_NUM / 2 + (LEVELS_PER_SIDE - 1 - level);
        }

        // 设置中心位及其相邻的 ±doppler_tolerance_bits 范围内的位
        int bit_start = std::max(0, center_bit_pos - doppler_tolerance_bits);
        int bit_end = std::min(static_cast<int>(HOUGHSLICE_DOPPLER_BIT_NUM) - 1, center_bit_pos + doppler_tolerance_bits);

        for (int bit_pos = bit_start; bit_pos <= bit_end; ++bit_pos)
        {
            size_t byte_idx = bit_pos / 8;
            size_t bit_in_byte = bit_pos % 8;
            mask_bytes[byte_idx] |= (1 << bit_in_byte);
        }

        // 遍历所有角度索引进行投票
        for (std::uint32_t angle_idx = angle_idx_start; angle_idx < angle_idx_end; ++angle_idx)
        {
            double theta = angle_idx * HOUGHSLICE_THETA_RESOLUTION_DEG * M_PI / 180.0;
            double distance = rel_x * std::cos(theta) + rel_y * std::sin(theta);

            // 计算距离索引
            int distance_idx = static_cast<int>((distance + 2 * HOUGHSLICE_CLUSTER_RADIUS_KM) / HOUGHSLICE_RHO_RESOLUTION_KM);
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
    std::vector<std::vector<std::array<size_t, 2>>> HoughSlice::process_extract_peak_from_hough_space(Slice &cluster) const
    {
        // 外层vector按doppler索引分组，内层vector存储该doppler下的(angle_idx, distance_idx)索引对
        std::vector<std::vector<std::array<size_t, 2>>> detected_lines_by_doppler(HOUGHSLICE_DOPPLER_BIT_NUM);

        for (size_t angle_idx = 0; angle_idx < HOUGH_THETA_DIM; ++angle_idx)
        {
            for (size_t distance_idx = 0; distance_idx < HOUGH_RHO_DIM; ++distance_idx)
            {
                // 提取HOUGHSLICE_BATCH_NUM个批次的共同投票位
                auto common_votes = peak_filter(angle_idx, distance_idx, cluster.vote_area);

                // 如果没有共同投票，跳过
                if (common_votes.none())
                {
                    continue;
                }

                LOG_DEBUG << "检测到峰值：angle_idx=" << angle_idx << ", distance_idx=" << distance_idx;

                // 遍历所有速度位
                for (size_t speed_bit = 0; speed_bit < HOUGHSLICE_DOPPLER_BIT_NUM; ++speed_bit)
                {
                    if (common_votes.get_bit(speed_bit))
                    {
                        // 直接按doppler索引分组存储(angle_idx, distance_idx)
                        detected_lines_by_doppler[speed_bit].push_back({angle_idx, distance_idx});
                    }
                }
            }
        }

        // 打印 (35,590) 周围一圈的峰值
        LOG_DEBUG << "===== (35,590) 周围一圈 =====";
        for (int di = -1; di <= 1; ++di)
        {
            for (int dj = -1; dj <= 1; ++dj)
            {
                int ni = 35 + di;
                int nj = 590 + dj;
                LOG_DEBUG << "(" << ni << "," << nj << "): " << cluster.vote_area[ni][nj];
            }
        }

        // 打印 (36,589) 周围一圈的峰值
        LOG_DEBUG << "===== (36,589) 周围一圈 =====";
        for (int di = -1; di <= 1; ++di)
        {
            for (int dj = -1; dj <= 1; ++dj)
            {
                int ni = 36 + di;
                int nj = 589 + dj;
                LOG_DEBUG << "(" << ni << "," << nj << "): " << cluster.vote_area[ni][nj];
            }
        }

        return detected_lines_by_doppler;
    }

    // 提取目标区域峰值
    BitArray<HOUGHSLICE_DOPPLER_BIT_NUM> HoughSlice::peak_filter(
        const size_t angle_idx, const size_t distance_idx,
        const std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area) const
    {
        const auto &cell = vote_area[angle_idx][distance_idx];

        constexpr size_t BYTES_PER_BATCH = HOUGHSLICE_DOPPLER_BIT_NUM / 8;

        // 分配两个缓冲区
        std::vector<uint8_t> buf1(BYTES_PER_BATCH);
        std::vector<uint8_t> buf2(BYTES_PER_BATCH);

        // 读取倒数第4批次
        cell.read_bytes(buf1.data(), 0, BYTES_PER_BATCH);

        // 创建结果对象并写入
        BitArray<HOUGHSLICE_DOPPLER_BIT_NUM> result;
        result.write_bytes(buf1.data(), 0, BYTES_PER_BATCH);

        // 依次处理后续批次
        for (size_t batch = 1; batch < HOUGHSLICE_BATCH_NUM; ++batch)
        {
            cell.read_bytes(buf2.data(), batch * BYTES_PER_BATCH, BYTES_PER_BATCH);

            BitArray<HOUGHSLICE_DOPPLER_BIT_NUM> current;
            current.write_bytes(buf2.data(), 0, BYTES_PER_BATCH);

            result &= current;
        }

        return result;
    }

    // 霍夫空间参数凝聚
    std::vector<std::array<double, 3>> HoughSlice::process_condense_detected_lines(
        const std::vector<std::vector<std::array<size_t, 2>>> &detected_lines) const
    {
        std::vector<std::array<double, 3>> condensed_lines;

        // 凝聚阈值（使用整数索引）
        const size_t THETA_CLUSTER_TOL = HOUGHSLICE_THETA_CLUSTER_TOL_DEG / HOUGHSLICE_THETA_RESOLUTION_DEG; // 角度索引差
        const size_t RHO_CLUSTER_TOL = HOUGHSLICE_RHO_CLUSTER_TOL_KM / HOUGHSLICE_RHO_RESOLUTION_KM;         // 距离索引差

        // 辅助函数：将doppler位索引转换为实际速度值
        auto doppler_bit_to_velocity = [](size_t bit) -> double
        {
            constexpr size_t HALF_BITS = HOUGHSLICE_DOPPLER_BIT_NUM / 2;
            constexpr double SPEED_STEP = 1.0 / HALF_BITS;

            if (bit < HALF_BITS)
            {
                return (-1.0 + bit * SPEED_STEP) * track_project::velocity_max;
            }
            else
            {
                size_t level = bit - HALF_BITS;
                return ((level + 1) * SPEED_STEP) * track_project::velocity_max;
            }
        };

        // 辅助函数：将索引转换为实际值
        auto idx_to_theta = [](size_t angle_idx) -> double
        {
            return angle_idx * HOUGHSLICE_THETA_RESOLUTION_DEG * M_PI / 180.0;
        };

        auto idx_to_rho = [](size_t distance_idx) -> double
        {
            return distance_idx * HOUGHSLICE_RHO_RESOLUTION_KM - 2 * HOUGHSLICE_CLUSTER_RADIUS_KM;
        };

        // 遍历每个多普勒速度组
        for (size_t doppler_bit = 0; doppler_bit < detected_lines.size(); ++doppler_bit)
        {
            const auto &lines_at_doppler = detected_lines[doppler_bit];
            if (lines_at_doppler.empty())
                continue;

            double doppler = doppler_bit_to_velocity(doppler_bit);

            // 对同一多普勒速度下的(angle_idx, distance_idx)进行聚类
            std::vector<bool> clustered(lines_at_doppler.size(), false);

            for (size_t i = 0; i < lines_at_doppler.size(); ++i)
            {
                if (clustered[i])
                    continue;

                const auto &idx_i = lines_at_doppler[i];
                size_t angle_i = idx_i[0];
                size_t rho_i = idx_i[1];

                // 找到所有属于同一聚类的点
                std::vector<size_t> cluster_indices = {i};
                clustered[i] = true;

                for (size_t j = i + 1; j < lines_at_doppler.size(); ++j)
                {
                    if (clustered[j])
                        continue;

                    const auto &idx_j = lines_at_doppler[j];
                    size_t angle_j = idx_j[0];
                    size_t rho_j = idx_j[1];

                    // 计算索引差（使用整数运算）
                    size_t angle_diff = angle_i > angle_j ? angle_i - angle_j : angle_j - angle_i;
                    // 角度索引需要考虑周期性,实际上m度和180-n度，当m和n非常小的时候任然是非常接近的，不能只看数值
                    if (angle_diff > HOUGH_THETA_DIM / 2)
                    {
                        angle_diff = HOUGH_THETA_DIM - angle_diff;
                    }
                    size_t rho_diff = rho_i > rho_j ? rho_i - rho_j : rho_j - rho_i;

                    if (angle_diff <= THETA_CLUSTER_TOL && rho_diff <= RHO_CLUSTER_TOL)
                    {
                        cluster_indices.push_back(j);
                        clustered[j] = true;
                    }
                }

                // 凝聚：取聚类中所有点的平均索引作为代表点
                if (cluster_indices.size() > 1)
                {
                    size_t angle_sum = 0, rho_sum = 0;
                    for (auto idx : cluster_indices)
                    {
                        angle_sum += lines_at_doppler[idx][0];
                        rho_sum += lines_at_doppler[idx][1];
                    }
                    // 整数平均，四舍五入
                    size_t angle_avg = (angle_sum + cluster_indices.size() / 2) / cluster_indices.size();
                    size_t rho_avg = (rho_sum + cluster_indices.size() / 2) / cluster_indices.size();

                    // 转换为实际值
                    double theta = idx_to_theta(angle_avg);
                    double rho = idx_to_rho(rho_avg);
                    condensed_lines.push_back({theta, rho, doppler});
                }
                else
                {
                    // 单个点，直接转换
                    double theta = idx_to_theta(angle_i);
                    double rho = idx_to_rho(rho_i);
                    condensed_lines.push_back({theta, rho, doppler});
                }
            }
        }

        return condensed_lines;
    }

    // 回溯点迹，只输出最新的4批次中符合条件的点迹，组成航迹（接口设定如此）
    void HoughSlice::process_backtrack_points(const std::vector<std::array<double, 3>> &detected_lines, const Slice &cluster,
                                              std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        const double DOPPLER_TOL = track_project::velocity_max / HOUGHSLICE_DOPPLER_BIT_NUM; // +-速度分辨率的一半
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

            // 对于HOUGHSLICE_BATCH_NUM批次数据，只检查最新的4个批次，依据多普勒速度和距离约束去寻找符合条件的点迹
            for (size_t batch = HOUGHSLICE_BATCH_NUM - 4; batch < HOUGHSLICE_BATCH_NUM && valid; ++batch)
            {
                const auto &points = cluster.point_list[batch];
                bool found = false;

                for (const auto &point : points)
                {
                    double rel_x = point.x - CENTER_X;
                    double rel_y = point.y - CENTER_Y;
                    double point_rho = rel_x * cos_theta + rel_y * sin_theta;

                    if (std::abs(point_rho - rho) < HOUGHSLICE_RHO_CLUSTER_TOL_KM && std::abs(point.doppler - doppler) < DOPPLER_TOL)
                    {
                        track[batch] = point;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                new_track.push_back(track);
            }
        }
    }

    void HoughSlice::clear_all()
    {
        ClustArea.clear_all();
    }
} // namespace track_project::trackinit
