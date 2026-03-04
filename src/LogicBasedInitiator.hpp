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

#ifndef _LOGIC_BASED_INITIATOR_HPP_
#define _LOGIC_BASED_INITIATOR_HPP_

#include "TrackInitBase.hpp"
#include <unordered_set>

namespace track_project::trackinit
{
    class LogicBasedInitiator : public TrackInitBase
    {
    public:
        // 友元测试类
        friend class test_LogicBasedInitiator;
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

    private:
        static constexpr size_t MAX_BINS = LOGIC_BASED_NUM_X_BINS * LOGIC_BASED_NUM_Y_BINS; // 最大允许的距离门数量乘以角度门数量
        // 每个bin中假设节点的最大数量，因为第三批次的结果直接输出出去了，因此每个节点中最多拥有1+子节点数量+孙节点数量的假设节点
        static constexpr size_t MAX_NODE_PER_BINS = LOGIC_BASED_MAX_CHILDREN_PER_PARENT_NODE * LOGIC_BASED_MAX_CHILDREN_PER_PARENT_NODE +
                                                    LOGIC_BASED_MAX_CHILDREN_PER_PARENT_NODE + 1;

    public:
        LogicBasedInitiator();
        virtual ~LogicBasedInitiator() noexcept = default;

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
        std::string get_name() const override { return "LogicBasedInitiator"; }

        /*****************************************************************************
         * @brief ⚠️特有：改装版逻辑法特有的函数，用于构建误差分布表格
         *
         * @param error_table 输入参数，存储每个bin的误差分布数据(sigma_x,sigma_y)
         *****************************************************************************/
        inline void build_error_distribution_table(std::vector<std::pair<double, double>> &error_table)
        {
            // 1. 检查大小是否完全匹配
            const size_t expected_size = MAX_BINS;
            const size_t actual_size = error_table.size();

            assert(actual_size == expected_size && "错误：误差分布表格大小不匹配！");

            // 3. 执行拷贝（深拷贝）
            error_distribution_table_ = error_table;
        }

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
        ProcessStatus extend_hypotheses(const TrackPoint &point, const std::unordered_set<HypothesisNode *> &parent_nodes);

        //**********************************工具函数**************************************** */
        /*****************************************************************************
         * @brief svd拟合直线，输入4个点迹，输出直线参数a,b,c，满足ax+by+c=0
         *
         * @param points 输入点迹，包含4个TrackPoint
         * @return std::array<double, 3>  对应变量a,b,c
         *****************************************************************************/
        std::array<double, 3> fit_line_svd(const std::vector<const TrackPoint *> &points) const;

        /*****************************************************************************
         * @brief 输入当前点迹的经纬度和DOPPLER，经过反推查询所有满足要求的假设树节点
         * 为了节省算力，假定相关点迹之间的误差分布函数不会有过大的变化，使用3.29SIGMA的门限对周围进行搜索
         * 未来如果要改进的话，应该是：
         * 1、成对于外推点所在位置，遍历周围一个巨大的空间进行遍历
         * 2. 对于每个BIN，反推其误差分布函数笼罩的范围是否包括当前点迹的外推点迹
         * 3. 如果包括，就将该BIN中的假设节点加入候选列表，并且对于该候选添加一个置信度，
         * 4. 置信度的结果为误差分布函数（离散化后）在该点的值。其中，我的专利提供了一种实时离散化误差分布函数的方法
         *
         * @param x 输入点迹的x坐标，单位米，东方向为x正方向
         * @param y 输入点迹的y坐标，单位米，北方向为y正方向
         * @param doppler 输入点迹的DOPPLER值
         * @return std::unordered_set<HypothesisNode *>
         * @version 0.3 xjl，修改了函数的返回值，遍历轮次提高的时候存在性能瓶颈,优化成这个了,以通用性为代价。。
         * @copyright Copyright (c) 2026
         * @author xjl (xjl20011009@126.com)
         *****************************************************************************/
        std::unordered_set<HypothesisNode *> query_nodes_by_points(double x, double y, double doppler) const;

        /*****************************************************************************
         * @brief 依据doppler和doppler计算船只的航向范围
         * 返回航向范围，北偏东，单位弧度，范围[-0.5PI,2.5PI)以确保航向连续性
         *
         * @param x 输入点迹的x坐标，单位米，东方向为x正方向
         * @param y 输入点迹的y坐标，单位米，北方向为y正方向
         * @param doppler 输入点迹的DOPPLER值，单位m/s，朝着雷达站是正方向
         * @return std::pair<double, double> 航向范围，北偏东，单位弧度，范围[-0.5PI,2.5PI)
         *****************************************************************************/
        std::pair<double, double> calculate_heading_range(double x, double y, double doppler) const;

        /*****************************************************************************
         * @brief 输入经纬度，输出对应的bin索引
         *
         * @param x 输入点迹的x坐标，单位米，东方向为x正方向
         * @param y 输入点迹的y坐标，单位米，北方向为y正方向
         * @return size_t bin索引，范围[0,MAX_BINS)，如果超出范围则返回MAX_BINS表示无效索引
         *****************************************************************************/
        size_t location_to_node_index(double x, double y) const;

    private:
        std::array<std::vector<TrackPoint>, 4> point_batches_;         // 追溯点迹区域，存储四批点迹
        std::array<std::vector<HypothesisNode>, 4> hypothesis_layers_; // 各个假设节点存储区域
        std::array<Timestamp, 4> batch_timestamps_;                    // 每批点迹的时间戳，单位秒

        // 假设列表，大小默认为MAX_BINS，不用ARRAY是因为容易栈溢出
        std::vector<std::vector<HypothesisNode *>> current_hypothesis_index_; // 当前假设索引表，
        std::vector<std::vector<HypothesisNode *>> history_hypothesis_index_; // 历史假设索引表

        // 误差分布表格，存储每个bin的误差分布数据(sigma_x,sigma_y)，大小默认为MAX_BINS，不用ARRAY是因为容易栈溢出
        std::vector<std::pair<double, double>> error_distribution_table_;
    };
}
#endif //_LOGIC_BASED_INITIATOR_HPP_