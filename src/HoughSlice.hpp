
/*****************************************************************************
 * @file hough.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 特点：
 * 1、对于每个点迹聚集区域，创建霍夫变换切面，每个切面只处理该区域内的点迹
 * 2、对于时间序列点迹，不同时刻、不同多普勒速度的信息保留在不同位中，每批次数据占据16bits，分别表示是否存在对应的多普勒速度
 * 3、对于角度约束，使用doppler和最大速度进行约束，同时基于多普勒速度一致性去寻找点迹
 *时间复杂度    分析：
 * 聚类生成：O(N)

霍夫投票：O(N × 角度分辨率)

峰值检测：O(θ维度 × ρ维度 × 16)

点迹回溯：O(峰值数 × 4 × 平均点迹数)
 * @version 0.2
 * @date 2025-12-17
 *
 * @copyright Copyright (c) 2025
 *
 *****************************************************************************/
#ifndef _HOUGH_SLICE_HPP_
#define _HOUGH_SLICE_HPP_

#include "TrackInitBase.hpp"
#include "ObjectPool.hpp"
#include <cstring>

namespace track_project::trackinit
{
    class SliceHough : public TrackInitBase
    {
        // 基准测试访问器友元，用于测试
        friend struct test_HoughSlice;

    private:
        // 航向角离散化参数，由于doppler位的引入，射线方向无意义，仅需"直线与y轴夹角":\alpha、"多普勒速度":doppler即可确定航向
        static constexpr size_t HOUGH_THETA_DIM = 180 / SLICEHOUGH_THETA_RESOLUTION_DEG;
        // 截距离散化参数，依据聚类判定直径2*1.41*RADIUS近似计算得到
        static constexpr size_t HOUGH_RHO_DIM = 4 * SLICEHOUGH_CLUSTER_RADIUS_KM / SLICEHOUGH_RHO_RESOLUTION_KM;

        struct Slice // 单个切面的内容
        {
            int current_batch_index;                           // 当前批次索引，从0开始
            std::array<std::vector<TrackPoint>, 4> point_list; // 历史点迹检索
            double center_x, center_y;                         // 聚类中心点坐标
            // 角度索引依据北偏东做分割，正北索引为0，是射线而非直线；截距索引依据负到正做分割，0点为 -2R，2R点为 +2R
            std::array<std::array<std::uint64_t, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> vote_area;

            // 为ObjectPool添加clear方法
            void clear()
            {
                current_batch_index = 0;

                // 清空point_list
                for (auto &vec : point_list)
                {
                    vec.clear();
                }

                // 聚类中心点清空
                center_x = 0.0;
                center_y = 0.0;

                // 清零vote_area
                vote_area = {};
            }
        };

    public:
        // 初始状态，预留20个聚类区域
        SliceHough() : ClustArea(20) {};
        virtual ~SliceHough() noexcept = default;

        /*****************************************************************************
         * @brief 检测直线，在cpp文件中
         * 当点迹过于分散，开辟了过多的空间，就反馈点迹过多
         *****************************************************************************/
        ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track) override;

        /*****************************************************************************
         * @brief 获取 name 对象
         *****************************************************************************/
        std::string get_name() const override { return "SliceHough"; }

        /*****************************************************************************
         * @brief 重置算法状态，清空所有数据
         *****************************************************************************/
        void clear_all() override;

    private:
        /*****************************************************************************
         * @brief 聚类算法，依据当前批次数据的情况，构建聚类，然后申请空间
         *
         * @param points 传进来的点迹
         *****************************************************************************/
        void process_cluster_generation(const std::vector<TrackPoint> &points);

        /*****************************************************************************
         * @brief 对于单个点迹进行霍夫变换投票
         *
         * @param batch 当前批次索引
         * @param point 当前点迹
         * @param it_clust 当前聚类区域
         *****************************************************************************/
        void process_point_for_hough_vote(size_t batch, const TrackPoint &point, Slice &it_clust);

        /*****************************************************************************
         * @brief 对于霍夫变换中的特定区域执行投票
         *
         * @param heading_start 起始航向角度，单位弧度
         * @param heading_end 结束航向角度，单位弧度
         * @param rel_x 点迹相对于聚类中心的x坐标，单位km
         * @param rel_y 点迹相对于聚类中心的y坐标，单位km
         * @param batch 当前批次索引
         * @param vote_area 投票区域，位于Slice结构体内
         *****************************************************************************/
        void vote_in_hough_space(const double heading_start, const double heading_end, const double doppler,
                                 const double rel_x, const double rel_y, const size_t batch,
                                 std::array<std::array<std::uint64_t, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area);

        /*****************************************************************************
         * @brief 从霍夫变换空间中提取峰值，返回检测到的直线参数列表
         *
         * @param cluster 单个霍夫变换切面
         * @return std::vector<std::array<double, 3>> 检测到的直线参数列表，每个元素为(theta, rho, doppler)
         *****************************************************************************/
        std::vector<std::array<double, 3>> process_extract_peak_from_hough_space(Slice &cluster) const;

        /*****************************************************************************
         * @brief 峰值过滤器，从霍夫变换空间中检测峰值，并推算峰值数量
         *
         * @param slice 单个切面
         * @param angle_idx 角度索引
         * @param distance_idx 距离索引
         * @param vote_area 投票区域，位于Slice结构体内
         * @return 返回一个峰值投票数，将速度均匀分成16份，每一位表示该位置上是否有对应速度点存在
         *****************************************************************************/
        std::uint16_t peak_filter(const size_t angle_idx, const size_t distance_idx,
                                  std::array<std::array<std::uint64_t, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area) const;

        /*****************************************************************************
         * @brief 回溯点迹，依据检测到的直线参数，从聚类中提取符合条件的点迹，组成航迹
         *
         * @param detected_lines 检测到的直线参数列表
         * @param cluster 单个霍夫变换切面
         * @param new_track 输出航迹列表
         *****************************************************************************/
        void process_backtrack_points(const std::vector<std::array<double, 3>> &detected_lines, const Slice &cluster,
                                      std::vector<std::array<TrackPoint, 4>> &new_track);

    private:
        ObjectPool<Slice> ClustArea;
    };
}

#endif
