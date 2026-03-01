/*****************************************************************************
 * @file LogicBased.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 特点，逻辑法，我计划引入一点MHT的思路，这就要求初始航迹的结构可能需要调整
 * @version 0.1
 * @date 2026-03-01
 *
 * @copyright Copyright (c) 2026
 *
 *****************************************************************************/

#ifndef _LOGIC_BASED_HPP_
#define _LOGIC_BASED_HPP_

#include "TrackInitBase.hpp"

namespace track_project::trackinit
{
    class LogicBasedTracker : public TrackInitBase
    {
    private:
        // 宏定义，结构体定义
        // 假设节点（Hypothesis Node）
        struct HypothesisNode
        {
            double heading; // 航向假设
            double vr_min;  // 径向速度下限
            double vr_max;  // 径向速度上限

            int depth; // 节点深度 (0-3)

            TrackPoint *obs;          // 关联的观测点迹
            HypothesisNode *ancestor; // 上一深度节点

            // 置信度（可选）
            double confidence;
        };

    public:
        LogicBasedTracker();
        virtual ~LogicBasedTracker() noexcept = default;

        /*****************************************************************************
         * @brief 检测直线，在cpp文件中
         * 当点迹过于分散，开辟了过多的空间，就反馈点迹过多
         *****************************************************************************/
        ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track) override;

        /*****************************************************************************
         * @brief 重置算法状态，清空所有数据
         *****************************************************************************/
        void clear_all() override;

        /*****************************************************************************
         * @brief 获取 name 对象
         *****************************************************************************/
        std::string get_name() const override { return "LogicBasedTracker"; }

    private:
        ProcessStatus extend_hypotheses(const std::vector<TrackPoint> &points);


        //**********************************工具函数**************************************** */
        /*****************************************************************************
         * @brief svd拟合直线，输入4个点迹，输出直线参数a,b,c，满足ax+by+c=0
         *
         * @param points 输入点迹，包含4个TrackPoint
         * @return std::array<double, 3>  对应变量a,b,c
         *****************************************************************************/
        std::array<double, 3> fit_line_svd(const std::vector<const TrackPoint *> &points) const;

        // TODO

    private:
        // 存储四批传入的点迹
        std::array<std::vector<TrackPoint>, 4> point_batches_;         // 追溯点迹区域，存储四批点迹
        std::array<std::vector<HypothesisNode *>, 4> hypothesis_tree_; // 追溯假设树，存储四批假设节点指针
        std::array<Timestamp, 4> batch_times_;                         // 每批点迹的时间戳，单位秒
    };
}
#endif //_LOGIC_BASED_HPP_