#include <catch2/catch_all.hpp>
#include "../src/LogicBasedInitiator.hpp"
#include "../include/ManagementService.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <memory>
#include <algorithm>
#include <cstring>
#include <set>
#include <thread>
#include <condition_variable>
#include <csignal>
#include <chrono>
#include <numeric>
#include <iomanip>

#include "../utils/Logger.hpp"

using namespace track_project;
using track_project::trackinit::LogicBasedInitiator;

// 为了匹配 LogicBasedInitiator.hpp 中的 "friend class test_LogicBasedInitiator;"
// 将测试辅助类放到同一个命名空间 track_project::trackinit 下
namespace track_project::trackinit
{

    class test_LogicBasedInitiator
    {
    public:
        // 构造函数持有引用
        explicit test_LogicBasedInitiator(LogicBasedInitiator &init) : initiator_(init) {}

        // 用于在单元测试或 gdb 中快速查看假设簇结构的摘要
        struct HypothesisSummary
        {
            double heading{};      // 航向假设
            double vr_min{};       // 径向速度下限
            double vr_max{};       // 径向速度上限
            std::size_t depth{};   // 节点深度
            double confidence{};   // 置信度
            bool has_associated{}; // 是否有关联观测
            bool has_parent{};     // 是否有父节点
        };

        //======================== 只读访问内部状态 ========================//

        // 获取某一批次点迹数量（0 最新，3 最旧）
        std::size_t get_point_batch_size(std::size_t batch_index) const
        {
            if (batch_index >= initiator_.point_batches_.size())
            {
                return 0;
            }
            return initiator_.point_batches_[batch_index].size();
        }

        // 获取某一层假设节点数量（0~3）
        std::size_t get_layer_node_count(std::size_t layer_index) const
        {
            if (layer_index >= initiator_.hypothesis_layers_.size())
            {
                return 0;
            }
            return initiator_.hypothesis_layers_[layer_index].size();
        }

        // 汇总当前索引表中非空 bin 数量（可以理解为有多少个“簇”被激活）
        std::size_t get_nonempty_current_bin_count() const
        {
            std::size_t cnt = 0;
            for (const auto &bucket : initiator_.current_hypothesis_index_)
            {
                if (!bucket.empty())
                {
                    ++cnt;
                }
            }
            return cnt;
        }

        // 汇总历史索引表中非空 bin 数量
        std::size_t get_nonempty_history_bin_count() const
        {
            std::size_t cnt = 0;
            for (const auto &bucket : initiator_.history_hypothesis_index_)
            {
                if (!bucket.empty())
                {
                    ++cnt;
                }
            }
            return cnt;
        }

        // 获取指定层、指定下标的假设摘要信息
        HypothesisSummary get_node_summary(std::size_t layer_index, std::size_t node_index) const
        {
            HypothesisSummary summary;

            if (layer_index >= initiator_.hypothesis_layers_.size())
            {
                return summary;
            }

            const auto &layer = initiator_.hypothesis_layers_[layer_index];
            if (node_index >= layer.size())
            {
                return summary;
            }

            const auto &node = layer[node_index];
            summary.heading = node.heading;
            summary.vr_min = node.vr_min;
            summary.vr_max = node.vr_max;
            summary.depth = node.depth;
            summary.confidence = node.confidence;
            summary.has_associated = (node.associated_obs != nullptr);
            summary.has_parent = (node.parent_node != nullptr);

            return summary;
        }

        //======================== 调试打印辅助函数 ========================//

        // 打印所有层次的假设树结构，方便在调试终端查看
        void print_hypothesis_tree(std::ostream &os = std::cout) const
        {
            os << "===== LogicBasedInitiator Hypothesis Tree =====\n";
            for (std::size_t layer_idx = 0; layer_idx < initiator_.hypothesis_layers_.size(); ++layer_idx)
            {
                const auto &layer = initiator_.hypothesis_layers_[layer_idx];
                os << "[Layer " << layer_idx << "] node_count=" << layer.size() << "\n";
                for (std::size_t i = 0; i < layer.size(); ++i)
                {
                    const auto &node = layer[i];
                    os << "  (#" << i << ") depth=" << node.depth
                       << ", heading=" << node.heading
                       << ", vr=[" << node.vr_min << ", " << node.vr_max << "]"
                       << ", conf=" << node.confidence;

                    if (node.associated_obs)
                    {
                        os << ", obs_lon=" << node.associated_obs->longitude
                           << ", obs_lat=" << node.associated_obs->latitude;
                    }

                    if (node.parent_node)
                    {
                        os << ", parent_depth=" << node.parent_node->depth;
                    }

                    os << "\n";
                }
            }
            os << "===============================================\n";
        }

        // 打印索引表中每个非空 bin 的大小（可以直观看到“簇”的分布情况）
        void print_index_summary(std::ostream &os = std::cout) const
        {
            os << "===== Current Hypothesis Index =====\n";
            for (std::size_t i = 0; i < initiator_.current_hypothesis_index_.size(); ++i)
            {
                const auto &bucket = initiator_.current_hypothesis_index_[i];
                if (!bucket.empty())
                {
                    os << "  bin=" << i << ", size=" << bucket.size() << "\n";
                }
            }

            os << "===== History Hypothesis Index =====\n";
            for (std::size_t i = 0; i < initiator_.history_hypothesis_index_.size(); ++i)
            {
                const auto &bucket = initiator_.history_hypothesis_index_[i];
                if (!bucket.empty())
                {
                    os << "  bin=" << i << ", size=" << bucket.size() << "\n";
                }
            }
            os << "===================================\n";
        }

    private:
        LogicBasedInitiator &initiator_;
    };

} // namespace track_project::trackinit

// 方便在测试代码中直接使用不带命名空间前缀的名称
using track_project::trackinit::test_LogicBasedInitiator;
