#include <catch2/catch_all.hpp>
#include "../src/LogicBasedInitiator.hpp"
#include "../include/ManagementService.hpp"
#include "../include/defsystem.h"
#include "Func_pointgen.hpp"
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
    static constexpr size_t MAX_BINS = LOGIC_BASED_NUM_X_BINS * LOGIC_BASED_NUM_Y_BINS; // 最大允许的距离门数量乘以角度门数量

    // test_LogicBasedInitiator.h
    class test_LogicBasedInitiator
    {
    public:
        explicit test_LogicBasedInitiator(LogicBasedInitiator &init) : initiator_(init) {}

        // based:
        bool saveErrorDistributionToDat(const std::string &filename) const
        {
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open())
            {
                return false;
            }

            // 直接写入所有sigma_x, sigma_y交替
            for (const auto &item : initiator_.error_distribution_table_)
            {
                double sigma_x = item.first;
                double sigma_y = item.second;
                file.write(reinterpret_cast<const char *>(&sigma_x), sizeof(double));
                file.write(reinterpret_cast<const char *>(&sigma_y), sizeof(double));
            }

            file.close();
            return true;
        }

        // 功能1：打印所有有内容的假设节点分布（包括内容）
        void printHypothesisDistribution() const
        {
            const auto &index = initiator_.current_hypothesis_index_;

            LOG_INFO << "========== 当前假设节点分布 ==========";

            for (size_t x = 0; x < LOGIC_BASED_NUM_X_BINS; ++x)
            {
                for (size_t y = 0; y < LOGIC_BASED_NUM_Y_BINS; ++y)
                {
                    size_t bin_index = x * LOGIC_BASED_NUM_Y_BINS + y; // 修正：x是列索引，y是行索引
                    const auto &nodes = index[bin_index];

                    if (nodes.empty())
                    {
                        continue;
                    }

                    LOG_INFO << "位置 [" << x << "," << y << "] (Bin Index: " << bin_index << ")，节点数量：" << nodes.size();
                    for (const auto *node : nodes)
                    {
                        LOG_INFO << "    节点详细信息：" << *node;
                    }
                }
            }

            LOG_INFO << "=====================================";
        }

        // 功能2：保存节点数量分布到DAT文件
        bool saveHypothesisDistributionToDat(const std::string &filename) const
        {
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open())
            {
                return false;
            }

            const auto &index = initiator_.current_hypothesis_index_;

            // 写入行数(Y方向)和列数(X方向)
            uint32_t rows = LOGIC_BASED_NUM_Y_BINS; // Y方向是行
            uint32_t cols = LOGIC_BASED_NUM_X_BINS; // X方向是列
            file.write(reinterpret_cast<const char *>(&rows), sizeof(uint32_t));
            file.write(reinterpret_cast<const char *>(&cols), sizeof(uint32_t));

            // 按顺序写入每个BIN的节点数量
            for (size_t x = 0; x < LOGIC_BASED_NUM_X_BINS; ++x)
            {
                for (size_t y = 0; y < LOGIC_BASED_NUM_Y_BINS; ++y)
                {
                    size_t bin_index = x * LOGIC_BASED_NUM_Y_BINS + y; // 修正的索引计算
                    uint32_t node_count = static_cast<uint32_t>(index[bin_index].size());
                    file.write(reinterpret_cast<const char *>(&node_count), sizeof(uint32_t));
                }
            }

            file.close();
            return true;
        }

        // 功能3：查询并打印指定区域内的所有假设节点详细信息
        void printHypothesisInRegion(size_t center_x, size_t center_y, int dx, int dy) const
        {
            const auto &index = initiator_.current_hypothesis_index_;

            LOG_INFO << "========== 查询区域 [" << center_x << "," << center_y << "] ± (" << dx << "," << dy << ") 的假设节点 ==========";

            // 计算区域范围
            int x_start = std::max(0, static_cast<int>(center_x) - dx);
            int x_end = std::min(static_cast<int>(LOGIC_BASED_NUM_X_BINS - 1),
                                 static_cast<int>(center_x) + dx);
            int y_start = std::max(0, static_cast<int>(center_y) - dy);
            int y_end = std::min(static_cast<int>(LOGIC_BASED_NUM_Y_BINS - 1),
                                 static_cast<int>(center_y) + dy);

            bool found_any = false;

            for (int x = x_start; x <= x_end; ++x)
            {
                for (int y = y_start; y <= y_end; ++y)
                {
                    size_t bin_index = x * LOGIC_BASED_NUM_Y_BINS + y;
                    const auto &nodes = index[bin_index];

                    if (nodes.empty())
                    {
                        LOG_INFO << "位置 [" << x << "," << y << "] (Bin Index: " << bin_index << ")：无节点";
                        continue;
                    }

                    found_any = true;
                    LOG_INFO << "位置 [" << x << "," << y << "] (Bin Index: " << bin_index << ")，节点数量：" << nodes.size();

                    // 打印每个节点的详细信息
                    for (size_t i = 0; i < nodes.size(); ++i)
                    {
                        const auto *node = nodes[i];
                        if (!node)
                            continue;

                        // 获取关联点迹的坐标（如果有的话）
                        double point_x = 0.0, point_y = 0.0;
                        if (node->associated_point)
                        {
                            point_x = node->associated_point->x;
                            point_y = node->associated_point->y;
                        }

                        LOG_INFO << "  节点[" << i << "]: depth=" << node->depth << ", confidence=" << node->confidence
                                 << ", heading=[" << node->heading_start << ", " << node->heading_end << "] rad";

                        LOG_INFO << "          parent=0x" << (void *)node->parent_node << ", point=(" << point_x << ", " << point_y << ")";

                        // 如果有父节点，打印父节点信息
                        if (node->parent_node)
                        {
                            LOG_INFO << "          parent depth=" << node->parent_node->depth << ", parent conf=" << node->parent_node->confidence;
                        }
                    }
                }
            }

            if (!found_any)
            {
                LOG_INFO << "区域内未找到任何假设节点";
            }

            LOG_INFO << "==========================================================";
        }

        // 辅助函数：打印单个位置的节点
        void printHypothesisAt(size_t x, size_t y) const
        {
            printHypothesisInRegion(x, y, 0, 0);
        }

    private:
        LogicBasedInitiator &initiator_;
    };

} // namespace track_project::trackinit

// 方便在测试代码中直接使用不带命名空间前缀的名称
using track_project::trackinit::test_LogicBasedInitiator;
using namespace track_project::trackinit;

// 用来看图的结果对不对的，也不算验证
TEST_CASE("单目标测试", "[FunctionalityCheck][single_track]")
{
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);

    // 误差分布表格保存
    // REQUIRE(tester.saveErrorDistributionToDat("../error_distribution.dat") == true);

    // 生成目标航迹
    std::vector<std::array<double, 4>> params = {{10.0, 5.0, 100, 50}}; // 初始化
    std::vector<TrackPoint> track_points = generate_target_points_xyv(0, params);

    // 处理点迹
    std::vector<std::array<TrackPoint, 4>> new_tracks;
    ProcessStatus status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第一批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_1.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第二批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_2.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第三批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_3.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第四批次结果";
    tester.printHypothesisDistribution(); // 打印当前假设节点分布
}

TEST_CASE("多目标测试", "[FunctionalityCheck][multi_track]")
{
    // 随机数种子
    unsigned int seed = Catch::getSeed();
    // unsigned int seed = 42;

    // 生成高斯分布的点迹
    auto points1 = generate_gaussian_points(10, 0,
                                            150, 10,
                                            150, 10,
                                            100.0, 50.0,
                                            seed++);
    auto points2 = generate_gaussian_points(10, 0,
                                            80, 10,
                                            80, 10,
                                            100.0, 50.0,
                                            seed++);
    auto points3 = generate_gaussian_points(10, 0,
                                            80, 10,
                                            150, 10,
                                            100.0, 50.0,
                                            seed++);
    auto points4 = generate_gaussian_points(10, 0,
                                            150, 10,
                                            80, 10,
                                            100.0, 50.0,
                                            seed++);
    std::vector<TrackPoint> points_all; // 所有点迹
    points_all.insert(points_all.end(), points1.begin(), points1.end());
    // points_all.insert(points_all.end(), points2.begin(), points2.end());
    // points_all.insert(points_all.end(), points3.begin(), points3.end());
    // points_all.insert(points_all.end(), points4.begin(), points4.end());

    // 创建算法实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);

    // 创建航迹管理器
    // track_project::ManagementService track_manager(0.3, 1.8, 0.3, 1.8); // 经纬度范围

    // // 绑定回调函数，显示航迹
    // alg.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
    //                        { track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    // 首次处理点迹
    // track_manager.clear_all_command(); // 清空状态

    // 处理点迹
    std::vector<std::array<TrackPoint, 4>> new_tracks;
    ProcessStatus status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第一批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_1.dat"); // 保存假设节点分布到DAT文件

    for (auto point : points_all)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第二批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_2.dat"); // 保存假设节点分布到DAT文件

    for (auto point : points_all)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第三批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_3.dat"); // 保存假设节点分布到DAT文件

    for (auto point : points_all)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第四批次结果";
    tester.printHypothesisDistribution(); // 打印当前假设节点分布
}

TEST_CASE("群目标测试", "[FunctionalityCheck][group_track]")
{
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);

    // 误差分布表格保存
    // REQUIRE(tester.saveErrorDistributionToDat("../error_distribution.dat") == true);

    // 生成目标航迹
    std::vector<std::array<double, 4>> params = {{10.0, 5.0, 100, 50},
                                                 {10.5, 5.5, 100, 50}}; // 初始化
    std::vector<TrackPoint> track_points = generate_target_points_xyv(0, params);

    // 处理点迹
    std::vector<std::array<TrackPoint, 4>> new_tracks;
    ProcessStatus status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第一批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_1.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第二批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_2.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第三批次结果";
    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_3.dat"); // 保存假设节点分布到DAT文件

    for (auto point : track_points)
    {
        point_update_cv(point, 1000); // 更新点迹位置，模拟1秒后的观测
    }
    status = initiator.process(track_points, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    LOG_INFO << "第四批次结果";
    tester.printHypothesisDistribution(); // 打印当前假设节点分布
}