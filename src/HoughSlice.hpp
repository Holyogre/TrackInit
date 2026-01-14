
/*****************************************************************************
 * @file hough.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 特点：
 * 1、对于每个点迹聚集区域，创建霍夫变换切面，每个切面只处理该区域内的点迹
 * 2、对于时间序列点迹，不同时刻的信息保留在不同位中，各占据8bit
 * 3、对于角度约束，使用doppler和最大速度进行约束
 *
 * @version 0.1
 * @date 2025-12-10
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
        friend struct SliceHoughBenchAccessor;

    private:
        // 角度和距离离散化参数
        static constexpr size_t ANGLE_BINS = 360 / SLICEHOUGH_ANGLE_RESOLUTION_DEG;
        static constexpr size_t DISTANCE_BINS = 2 * SLICEHOUGH_CLUSTER_RADIUS_KM / SLICEHOUGH_DIST_RESOLUTION_KM;

        struct Slice // 单个切面的内容
        {
            int current_batch_index;                           // 当前批次索引，从0开始
            std::array<std::vector<TrackPoint>, 4> point_list; // 历史点迹检索
            double center_x, center_y;                         // 聚类中心点坐标
            // 角度索引依据南偏东做分割，正南索引为0，是射线而非直线；截距索引依据负到正做分割，0点为 -R，2R点为 +R
            std::uint32_t vote_area[static_cast<std::uint32_t>(ANGLE_BINS)][static_cast<std::uint32_t>(DISTANCE_BINS)];

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
                std::memset(vote_area, 0, sizeof(vote_area));
            }
        };

    public:
        // 共计四批次，每个批次允许新增5个聚类
        SliceHough() : ClustArea(5 * 4) {};
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
        void clust_gen(const std::vector<TrackPoint> &points);

        /*****************************************************************************
         * @brief 峰值过滤器，从霍夫变换空间中检测峰值，并推算峰值数量
         *
         * @param slice 单个切面
         * @param angle_idx 角度索引
         * @param distance_idx 距离索引
         * @return
         *****************************************************************************/
        std::uint32_t peak_filter(const Slice &slice, size_t angle_idx, size_t distance_idx) const;

        /*****************************************************************************
         * @brief 在霍夫空间中投票
         *
         * @param heading_start 起始航向角度，单位弧度
         * @param heading_end 结束航向角度，单位弧度
         *****************************************************************************/
        void voteInHoughSpace(double heading_start, double heading_end, double rel_x, double rel_y, size_t batch, Slice &it_clust);

    private:
        ObjectPool<Slice> ClustArea;
    };
}

#endif
