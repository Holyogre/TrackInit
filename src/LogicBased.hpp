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
#include <unordered_map>

namespace track_project::trackinit
{
    class LogicBasedTracker : public TrackInitBase
    {
    private:
        static constexpr size_t MAX_BINS = LOGIC_BASED_MAX_R_BINS * LOGIC_BASED_MAX_THETA_BINS; // 最大允许的距离门数量乘以角度门数量

        // 假设节点（Hypothesis Node）
        struct HypothesisNode
        {
            double heading; // 航向假设
            double vr_min;  // 径向速度下限
            double vr_max;  // 径向速度上限

            size_t depth; // 节点深度(0-3)

            TrackPoint *associated_obs;  // 关联的观测点迹
            HypothesisNode *parent_node; // 上一深度节点

            // 置信度，我计划从两个方面计算，一个时通过雷达视界的置信度参考图，另一个通过DOPPLER的分布情况
            double confidence;
        };

    public:
        LogicBasedTracker();
        virtual ~LogicBasedTracker() noexcept = default;

        /*****************************************************************************
         * @brief 检测直线，在cpp文件中
         * 当点迹过于分散，开辟了过多的空间，就反馈点迹过多
         *****************************************************************************/
        ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_tracks) override;

        /*****************************************************************************
         * @brief 重置算法状态，清空所有数据
         *****************************************************************************/
        void clear_all() override;

        /*****************************************************************************
         * @brief 获取 name 对象
         *****************************************************************************/
        std::string get_name() const override { return "LogicBasedTracker"; }

    private:
        /*****************************************************************************
         * @brief 移动批次数据和假设树，清空最旧的批次数据和对应的假设树，
         * 更新索引表,确保[0]索引对应的总是最新的数据
         *****************************************************************************/
        void shift_batches_and_hypotheses();

        /*****************************************************************************
         * @brief 扩展假设树，若有点迹未被存放到任何假设树中，就以这些点迹为基础生成新的假设树
         *
         * @param points
         * @return ProcessStatus 错误码，SUCCESS表示成功，其他值表示具体错误
         *****************************************************************************/
        ProcessStatus extend_hypotheses(const TrackPoint &point, const std::vector<HypothesisNode *> &parent_nodes);

        //**********************************工具函数**************************************** */
        /*****************************************************************************
         * @brief svd拟合直线，输入4个点迹，输出直线参数a,b,c，满足ax+by+c=0
         *
         * @param points 输入点迹，包含4个TrackPoint
         * @return std::array<double, 3>  对应变量a,b,c
         *****************************************************************************/
        std::array<double, 3> fit_line_svd(const std::vector<const TrackPoint *> &points) const;

        /*****************************************************************************
         * @brief 输入当前点迹的经纬度，查询所有满足要求的假设树节点
         *
         * @param longitude 输入点迹的经度，单位度
         * @param latitude 输入点迹的纬度，单位度
         * @return std::vector<HypothesisNode *> 所有可能的假设列表
         *****************************************************************************/
        std::vector<HypothesisNode *> query_nodes_by_location(double lon, double lat) const;

    private:
        std::array<std::vector<TrackPoint>, 4> point_batches_;         // 追溯点迹区域，存储四批点迹
        std::array<std::vector<HypothesisNode>, 4> hypothesis_layers_; // 各个假设节点存储区域
        std::array<Timestamp, 4> batch_timestamps_;                    // 每批点迹的时间戳，单位秒

        std::array<std::vector<HypothesisNode *>, MAX_BINS> current_hypothesis_index_; // 当前假设索引表
        std::array<std::vector<HypothesisNode *>, MAX_BINS> history_hypothesis_index_; // 历史假设索引表
    };
}
#endif //_LOGIC_BASED_HPP_