#include "HoughSlice.hpp"
#include "Func_dbscan.hpp"
#include "../include/defsystem.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace track_project::trackinit
{
    ProcessStatus SliceHough::process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        // 清空输出航迹
        new_track.clear();

        // STEP1：对于新来的点迹生成聚类
        clust_gen(points);

        // 读取各个聚类区域，进行霍夫变换投票
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (auto &it_clust : clustAera_list)
        {
            if (it_clust->current_batch_index != 3)
            {
                continue; // 批次数不足，跳过
            }

            // 对该聚类区域内的点迹进行霍夫变换投票
            for (size_t batch = 0; batch <= 3; ++batch)
            {
                for (const auto &point : it_clust->point_list[batch])
                {
                    // 计算点迹相对于聚类中心的坐标
                    double rel_x = point.x - it_clust->center_x;
                    double rel_y = point.y - it_clust->center_y;

                    // 依据doppler和velocity_max计算所有可能的航向角度序列
                    double speed = std::abs(point.doppler); // 取绝对值，单位m/s
                    if (speed > track_project::velocity_max)
                    {
                        continue; // 速度超过最大值，跳过该点迹
                    }
                    double angle_rad = std::acos(speed / track_project::velocity_max); // 计算夹角，弧度，必定为正数
                    // 计算两个可能的航向角度（南偏东）
                    double beta = std::atan2(rel_y, rel_x);
                    if (point.doppler > 0) // 速度是正值时是靠近雷达站的，此时需要对beta加π；操作后值域为(0,2*pi)，无需归一化
                    {
                        beta += M_PI;
                        if (beta >= 2 * M_PI)
                        {
                            beta -= 2 * M_PI;
                        }
                    }
                    else // 不然，值域为(-pi,pi)，需要归一化到(0,2*pi)
                    {
                        if (beta < 0)
                        {
                            beta += 2 * M_PI;
                        }
                    }
                    double heading1 = beta - angle_rad;
                    double heading2 = beta + angle_rad;
                    //  区间为(heading1,heading2)∪(heading3,heading4)，且heading1<=heading2<=heading3<=heading4
                    if (heading1 < 0) // 仅有可能heading1小于0
                    {
                        double heading3 = 0.0, heading4 = 2 * M_PI;
                        heading3 = heading1 + 2 * M_PI; // 调整到正区间
                        heading4 = 2 * M_PI;
                        heading1 = 0.0;
                        heading2 = heading2;
                        voteInHoughSpace(heading1, heading2, rel_x, rel_y, batch, *it_clust);
                        voteInHoughSpace(heading3, heading4, rel_x, rel_y, batch, *it_clust);
                    }
                    else if (heading2 > 2 * M_PI) // 仅有可能heading2大于2π
                    {
                        double heading3 = 0.0, heading4 = 2 * M_PI;
                        heading3 = heading1;
                        heading4 = 2 * M_PI;
                        heading1 = 0.0;
                        heading2 = heading2 - 2 * M_PI;
                        voteInHoughSpace(heading1, heading2, rel_x, rel_y, batch, *it_clust);
                        voteInHoughSpace(heading3, heading4, rel_x, rel_y, batch, *it_clust);
                    }
                    else
                    {
                        voteInHoughSpace(heading1, heading2, rel_x, rel_y, batch, *it_clust);
                    }

                    // TODO 做成一个投票函数，解耦合，遍历所有角度，计算对应的距离截距
                }
            }
        }

        // TODO 准备明天写峰值检测和航迹生成部分

        return ProcessStatus::SUCCESS;
    }

    void SliceHough::voteInHoughSpace(double heading_start, double heading_end, double rel_x, double rel_y, size_t batch, Slice &it_clust)
    {
        // 计算起始和结束的角度索引
        std::uint32_t angle_idx_start = static_cast<std::uint32_t>(heading_start * 180.0 / M_PI / SLICEHOUGH_ANGLE_RESOLUTION_DEG);
        std::uint32_t angle_idx_end = static_cast<std::uint32_t>(heading_end * 180.0 / M_PI / SLICEHOUGH_ANGLE_RESOLUTION_DEG);

        // 遍历所有聚类区域，进行投票
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (std::uint32_t angle_idx = angle_idx_start; angle_idx < angle_idx_end; ++angle_idx)
        {

            double theta = angle_idx * SLICEHOUGH_ANGLE_RESOLUTION_DEG * M_PI / 180.0; // 转为弧度
            double distance = rel_x * std::cos(theta) + rel_y * std::sin(theta);

            // 计算距离索引
            int distance_idx = static_cast<int>((distance + SLICEHOUGH_CLUSTER_RADIUS_KM) / SLICEHOUGH_DIST_RESOLUTION_KM);
            if (distance_idx < 0 || distance_idx >= static_cast<int>(DISTANCE_BINS))
            {
                continue; // 距离索引越界，跳过
            }

            // 在对应位置投票，使用8位存储不同批次的信息
            it_clust.vote_area[angle_idx][distance_idx] += (1 << (batch * 8));
        }
    }

    // 先剔除在聚类中心的点，再剔除离群点，最后对剩余点进行聚类
    void SliceHough::clust_gen(const std::vector<TrackPoint> &points)
    {
        // 用于记录当前点迹是否被聚类
        std::vector<bool> visited(points.size(), false);

        // 循环访问聚类，看能不能存放到之前的聚类里面去
        std::vector<Slice *> clustAera_list = ClustArea.get_allocated_ptrs();
        for (const auto &it_clust : clustAera_list)
        {
            //  检测点迹是否在中心点附近
            for (auto it_point = points.begin(); it_point != points.end(); ++it_point)
            {
                const double dx = it_point->x - it_clust->center_x;
                const double dy = it_point->y - it_clust->center_y;
                const double dist_square = dx * dx + dy * dy;
                if (dist_square <= SLICEHOUGH_CLUSTER_RADIUS_KM * SLICEHOUGH_CLUSTER_RADIUS_KM)
                {
                    // 点迹在聚类范围内，加入该聚类
                    size_t batch = it_clust->current_batch_index + 1; // 下一批次数据存放位置
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

        // 拆分聚类
        std::vector<std::vector<size_t>> _clust = dbscan(unvisited_points, SLICEHOUGH_CLUSTER_RADIUS_KM / 2.0, 3);

        // 计算相关参数，申请新的聚类空间
        for (const auto &cluster : _clust)
        {
            if (cluster.size() < 0) // 忽略空集
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

    std::uint32_t SliceHough::peak_filter(const Slice &slice, size_t angle_idx, size_t distance_idx) const
    {
        int current_value = slice.vote_area[angle_idx][distance_idx];
        std::uint8_t peak_count = 0; // 记录每8位数据的最小值

        // 检测是否每个八位都是非零的
        for (int batch = 0; batch < 4; ++batch)
        {
            std::uint8_t batch_peak_value = (current_value >> (batch * 8)) & 0xFF;
            // 若有一个批次的值为零，直接返回不存在航迹
            if (batch_peak_value == 0)
            {
                return 0;
            }

            // 第一批次的时候，直接赋值
            if (batch == 0)
            {
                peak_count = batch_peak_value;
                continue;
            }

            // 此后计算每个批次的值，并更新peak_count
            peak_count = std::min(peak_count, batch_peak_value);
        }

        return static_cast<std::uint32_t>(peak_count);
    }

    void SliceHough::clear_all()
    {
        ClustArea.clear_all();
    }

}
