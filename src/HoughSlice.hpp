
/*****************************************************************************
 * @file hough.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 特点：
 * 1、对于每个点迹聚集区域，创建霍夫变换切面，每个切面只处理该区域内的点迹
 * 2、对于时间序列点迹，不同时刻、不同多普勒速度的信息保留在不同位中，每批次数据占据16bits，分别表示是否存在对应的多普勒速度
 * 3、对于角度约束，使用doppler和最大速度进行约束，同时基于多普勒速度一致性去寻找点迹
 * 4、“使用K批次点迹生成4个航迹”，其中K可以利用def_init头文件中的HOUGHSLICE_BATCH_NUM修改，
 *时间复杂度    分析：
 * 聚类生成：O(N)

霍夫投票：O(N × 角度分辨率)

峰值检测：O(θ维度 × ρ维度 × 16)

点迹回溯：O(峰值数 × 4 × 平均点迹数)
 * @version 0.3 我确信霍夫变换在航迹起始根本不好用，别用了，原因：
 1. 有SOG的话效果可能还行，DOPPLER根本没法用，近距离处误差太大了
 2. 而且启航时间太短了，THETA的分辨率糊成一坨，引入误差之后几乎没法计算，我们不可能容忍这么多关联时间的
 3. 航迹回溯的时间复杂度已经接近逻辑法了，还是别用霍夫变换了
 * @date 2025-12-17
 *
 * @copyright Copyright (c) 2025
 *
 *****************************************************************************/
#ifndef _HOUGH_SLICE_HPP_
#define _HOUGH_SLICE_HPP_

#include <cstring>

#include "TrackInitBase.hpp"
#include "ObjectPool.hpp"
#include "BitArray.hpp"
#include "LatestKBuffer.hpp"

namespace track_project::trackinit
{
    class HoughSlice : public TrackInitBase
    {
        // 基准测试访问器友元，用于测试
        friend struct test_HoughSlice;

    private:
        // 航向角离散化参数，由于doppler位的引入，射线方向无意义，仅需"直线与y轴夹角":\alpha、"多普勒速度":doppler即可确定航向
        static constexpr size_t HOUGH_THETA_DIM = 180 / HOUGHSLICE_THETA_RESOLUTION_DEG;
        // 截距离散化参数，搜索空间是圆，所以RHO的极限也是R
        static constexpr size_t HOUGH_RHO_DIM = 2 * HOUGHSLICE_CLUSTER_RADIUS_KM / HOUGHSLICE_RHO_RESOLUTION_KM;
        // 霍夫变换空间的位数，等于批次数*每批次的速度位数
        static constexpr size_t HOUGH_VOTE_BIT_NUM = HOUGHSLICE_BATCH_NUM * HOUGHSLICE_DOPPLER_BIT_NUM;

        struct Slice // 单个切面的内容，雷达站位于x=0,y=0处，坐标单位为km，霍夫变换的角度是东偏北（标准极坐标，不要和cog搞混）
        {
            int current_batch_index;                                              // 当前批次索引，从0开始
            std::array<std::vector<TrackPoint>, HOUGHSLICE_BATCH_NUM> point_list; // 历史点迹检索
            double center_x, center_y;                                            // 聚类中心点坐标
            // 截距索引依据负到正做分割，0点为 -R，2R点为 +R
            // 角度索引为法线方向，范围为(-90,90)，0点为y轴正方向，90点为x轴正方向，-90点为x轴负方向，按索引恰好为HEADING方向
            std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> vote_area;

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
                for (auto &theta_array : vote_area)
                {
                    for (auto &bit_array : theta_array)
                    {
                        bit_array.clear();
                    }
                }
            }
        };

    public:
        // 初始状态，预留20个聚类区域
        HoughSlice() : time_buffer(HOUGHSLICE_BATCH_NUM), ClustArea(20) {};
        virtual ~HoughSlice() noexcept = default;

        /*****************************************************************************
         * @brief 检测直线，在cpp文件中
         * 当点迹过于分散，开辟了过多的空间，就反馈点迹过多
         *****************************************************************************/
        ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track) override;

        /*****************************************************************************
         * @brief 获取 name 对象
         *****************************************************************************/
        std::string get_name() const override { return "HoughSlice"; }

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
         * @param doppler 多普勒速度，单位m/s
         * @param rel_x 点迹相对于聚类中心的x坐标，单位km
         * @param rel_y 点迹相对于聚类中心的y坐标，单位km
         * @param batch 当前批次索引
         * @param doppler_tolerance_bits 多普勒容差位数，在中心位两侧各设置此数量的位
         * @param vote_area 投票区域，位于Slice结构体内
         *****************************************************************************/
        void vote_in_hough_space(const double heading_start, const double heading_end, const double doppler,
                                 const double rel_x, const double rel_y, const size_t batch,
                                 const int doppler_tolerance_bits,
                                 std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area);

        /*****************************************************************************
         * @brief 从霍夫变换空间中提取峰值，返回按多普勒速度分组的检测直线参数
         *
         * @param cluster 单个霍夫变换切面
         * @return std::vector<std::vector<std::array<size_t, 2>>>
         *         返回值是一个二维向量：
         *         - 第一维索引对应多普勒速度位 (0 ~ HOUGHSLICE_DOPPLER_BIT_NUM-1)
         *         - 第二维是该多普勒速度下的所有检测直线，每个元素为(theta, rho)的索引对
         *
         * @note 多普勒速度位的具体含义：
         *       - 低半区 (0 ~ HOUGHSLICE_DOPPLER_BIT_NUM/2 - 1): 负速度，从最大负到最小负
         *       - 高半区 (HOUGHSLICE_DOPPLER_BIT_NUM/2 ~ HOUGHSLICE_DOPPLER_BIT_NUM-1): 正速度，从最小正到最大正
         *
         * @version 0.2 修改返回值类型，将多普勒维度提取到外层索引，便于后续按速度分组处理
         *****************************************************************************/
        std::vector<std::vector<std::array<size_t, 2>>> process_extract_peak_from_hough_space(Slice &cluster) const;

        /*****************************************************************************
         * @brief 峰值过滤器，从霍夫变换空间中检测峰值，并推算峰值数量
         *
         * @param slice 单个切面
         * @param angle_idx 角度索引
         * @param distance_idx 距离索引
         * @param vote_area 投票区域，位于Slice结构体内
         * @return 返回一个峰值投票数，将速度均匀分成16份，每一位表示该位置上是否有对应速度点存在
         *****************************************************************************/
        BitArray<HOUGHSLICE_DOPPLER_BIT_NUM> peak_filter(
            const size_t angle_idx, const size_t distance_idx,
            const std::array<std::array<BitArray<HOUGH_VOTE_BIT_NUM>, HOUGH_RHO_DIM>, HOUGH_THETA_DIM> &vote_area) const;

        /*****************************************************************************
         * @brief 点迹凝聚，对霍夫空间中的点迹拥簇问题，进行适当的凝聚（主要是RHO上的），减少重复航迹生成
         *
         * @param detected_lines 按多普勒速度分组的检测直线参数
         *                       第一维：多普勒速度位索引
         *                       第二维：该速度下的(theta, rho)参数对列表
         * @return std::vector<std::array<double, 3>> 凝聚后的直线参数列表 (theta, rho, doppler)
         *****************************************************************************/
        std::vector<std::array<double, 3>> process_condense_detected_lines(const std::vector<std::vector<std::array<size_t, 2>>> &detected_lines) const;

        /*****************************************************************************
         * @brief 回溯点迹，依据检测到的直线参数，从聚类中提取符合条件的点迹，组成航迹
         *
         * @param detected_lines 检测到的直线参数列表，分别是(theta, rho, doppler)
         * @param cluster 单个霍夫变换切面
         * @param new_track 输出航迹列表
         *****************************************************************************/
        void process_backtrack_points(const std::vector<std::array<double, 3>> &detected_lines, const Slice &cluster,
                                      std::vector<std::array<TrackPoint, 4>> &new_track);

        /*****************************************************************************
         * @brief 依据不准确的DX,DY，直线与X轴的某个夹角，来计算航向
         *
         * @param line_angle 直线与x轴的夹角，单位弧度，范围[0, π)，由霍夫变换空间的theta索引转换得到
         * @param points 组成航迹的四个点迹
         * @return std::pair<bool, double> 是否成功计算航向及真实航向
         *****************************************************************************/
        inline std::pair<bool, double> unwrapHeadingFromDisplacement(double line_angle, const std::array<TrackPoint, 4> &points) const
        {
            // 用位移向量估算真实航向
            double dx = points[3].x - points[0].x;
            double dy = points[3].y - points[0].y;

            // 处理位移为零的情况
            const double EPS = 1e-10;
            if (std::abs(dx) < EPS && std::abs(dy) < EPS)
            {
                return {false, 0.0};
            }

            // line_angle 是直线方向 [0, π)，真实航向要么是它，要么是它 + π
            double candidate1 = line_angle;        // [0, π)
            double candidate2 = line_angle + M_PI; // [π, 2π)
            if (candidate2 >= 2 * M_PI)
                candidate2 -= 2 * M_PI;

            // 用位移的象限做粗粒度判断（抗噪声）
            bool moving_right = dx > 0;
            bool moving_up = dy > 0;

            // 根据运动方向选择候选
            double heading_true;

            if (moving_right && moving_up)
            { // 第一象限：目标朝右上
                // 真实航向应该在 0~90° 之间
                heading_true = (candidate1 <= M_PI / 2) ? candidate1 : candidate2;
            }
            else if (!moving_right && moving_up)
            { // 第二象限：左上
                // 真实航向应该在 90~180° 之间
                heading_true = (candidate1 > M_PI / 2 && candidate1 < M_PI) ? candidate1 : candidate2;
                // 确保结果在 90~180°
                if (heading_true < M_PI / 2)
                    heading_true += M_PI;
            }
            else if (!moving_right && !moving_up)
            { // 第三象限：左下
                // 真实航向应该在 180~270° 之间
                heading_true = (candidate1 > M_PI) ? candidate1 : candidate2;
                // 确保结果在 180~270°
                if (heading_true < M_PI)
                    heading_true += M_PI;
            }
            else
            { // 第四象限：右下 (moving_right && !moving_up)
                // 真实航向应该在 270~360° 之间
                heading_true = (candidate1 > 3 * M_PI / 2) ? candidate1 : candidate2;
                // 确保结果在 270~360°
                if (heading_true < 3 * M_PI / 2)
                    heading_true += M_PI;
            }

            // 归一化到 [0, 2π)
            heading_true = std::fmod(heading_true, 2 * M_PI);
            if (heading_true < 0)
                heading_true += 2 * M_PI;

            return {true, heading_true};
        }

    private:
        ObjectPool<Slice> ClustArea;

        track_project::LatestKBuffer<track_project::Timestamp> time_buffer; // 时间戳缓冲区，存储每批次数据的时间戳
    };
}

#endif
