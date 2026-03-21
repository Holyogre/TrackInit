#include <catch2/catch_all.hpp>
#include "../src/LogicBasedInitiator.hpp"
#include "../include/ManagementService.hpp"
#include "../include/defsystem.h"
#include "Func_pointgen.hpp"
#include "Func_noise_pointgen.hpp"
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
#include <unistd.h>

// 检测时间间隔，单位秒
constexpr double TIME_INTERVAL_S = 10.0;

std::atomic<bool> g_running{true};
std::condition_variable g_cv;
std::mutex g_mutex;

void signal_handler(int)
{
    g_running = false;
    g_cv.notify_one();
}

void wait_seconds(int seconds)
{
    std::unique_lock<std::mutex> lock(g_mutex);
    g_cv.wait_for(lock, std::chrono::seconds(seconds), []
                  { return !g_running; });
}

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
                        LOG_INFO << "   航迹索引" << "  depth=" << node->depth << ", confidence=" << node->confidence
                                 << ", 航向范围=" << " [" << node->heading_start / 3.14159 * 180 << ", " << node->heading_end / 3.14159 * 180
                                 << "] deg";

                        // 如果有多个节点，打印详细信息
                        if (nodes.size() == 1)
                        {
                            continue; // 只有一个节点，说明关联的不错，就不打印了
                        }
                        // 出现了重复航迹，向上追溯所有父节点，打印每个depth关联的点迹
                        const LogicBasedInitiator::HypothesisNode *current_node = node;
                        size_t node_depth = current_node->depth;
                        while (current_node != nullptr)
                        {
                            if (current_node->associated_point)
                            {
                                // 获取当前节点depth对应的点迹批次
                                size_t batch_index = node_depth - current_node->depth; // depth 对应点迹批次的索引
                                if (batch_index < initiator_.point_batches_.size() && !initiator_.point_batches_[batch_index].empty())
                                {
                                    const auto *begin = &initiator_.point_batches_[batch_index][0];
                                    const auto *current_point = current_node->associated_point;

                                    // 使用指针减法计算索引（避免std::distance的类型问题）
                                    ptrdiff_t point_index = current_point - begin;

                                    LOG_INFO << "   点迹索引:" << point_index
                                             << ", 经度=" << current_point->longitude
                                             << ", 纬度=" << current_point->latitude
                                             << ", x=" << current_point->x
                                             << ", y=" << current_point->y
                                             << ", sog=" << current_point->sog
                                             << ", cog=" << current_point->cog;
                                }
                                else
                                {
                                    LOG_INFO << "       无效的批次索引，depth=" << current_node->depth;
                                }
                            }
                            else
                            {
                                LOG_INFO << "       depth=" << current_node->depth << " 无关联点迹";
                            }

                            current_node = current_node->parent_node;
                        }
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

            for (size_t x = x_start; x <= x_end; ++x)
            {
                for (size_t y = y_start; y <= y_end; ++y)
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
    std::signal(SIGINT, signal_handler);

    double pos_noise_sigma_km = 2.0;  // 位置噪声标准差10米（0.01km）
    double doppler_noise_sigma = 0.5; // 多普勒噪声标准差0.1m/s

    unsigned int seed[3] = {4, 232, 3424}; // 三个随机数种子
    seed[0] = rand() % 10000;
    seed[1] = rand() % 10000;
    seed[2] = rand() % 10000;

    std::vector<std::array<double, 4>> params = {{70.0, 70.0, 1.0, 1.0}}; // 初始化
    std::vector<TrackPoint> track_points = generate_target_points_xyv(0, params);

    std::vector<TrackPoint> points_all = track_points; // 所有点迹

    // 创建算法实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);

    // 创建航迹管理器
    track_project::ManagementService track_manager(0.45, 0.75, 0.45, 0.75); // 经纬度范围

    // 绑定回调函数，显示航迹
    initiator.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                                 { track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    // 首次处理点迹
    track_manager.clear_all_command(); // 清空状态

    //*****************************************第一次处理数据***********************************************/
    LOG_INFO << "第一批次结果";
    track_manager.draw_point_command(points_all); // 绘制点迹

    std::vector<std::array<TrackPoint, 4>> new_tracks;
    ProcessStatus status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_1.dat"); // 保存假设节点分布到DAT文件

    //*****************************************第二次处理数据***********************************************/
    // 更新点迹位置并装载（带噪声）
    for (auto &p : points_all)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed[0], pos_noise_sigma_km, doppler_noise_sigma); // 使用带噪声的更新，位置噪声标准差10米（0.01km），多普勒噪声标准差0.1m/s
    }

    LOG_INFO << "第二批次结果";
    track_manager.draw_point_command(points_all); // 绘制点迹
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_2.dat"); // 保存假设节点分布到DAT文件

    //*****************************************第三次处理数据***********************************************/
    // 更新点迹位置并装载（带噪声）
    for (auto &p : points_all)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed[1], pos_noise_sigma_km, doppler_noise_sigma); // 使用带噪声的更新，位置噪声标准差10米（0.01km），多普勒噪声标准差0.1m/s
    }

    LOG_INFO << "第三批次结果";
    track_manager.draw_point_command(points_all); // 绘制点迹
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    tester.printHypothesisDistribution();                                       // 打印当前假设节点分布
    tester.saveHypothesisDistributionToDat("../hypothesis_distribution_3.dat"); // 保存假设节点分布到DAT文件

    //*****************************************第四次处理数据***********************************************/
    // 更新点迹位置并装载（带噪声）
    for (auto &p : points_all)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed[2], pos_noise_sigma_km, doppler_noise_sigma);
    }

    LOG_INFO << "第四批次结果";
    track_manager.draw_point_command(points_all); // 绘制点迹
    sleep(1);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    tester.printHypothesisDistribution(); // 打印当前假设节点分布

    // 保持窗口显示，直到用户按下CTRL+C
    while (g_running)
    {
        wait_seconds(5); // 可被 CTRL+C 中断的等待
    }
}

TEST_CASE("多目标测试", "[FunctionalityCheck][multi_track]")
{
    std::signal(SIGINT, signal_handler); // 添加信号处理

    // 随机数种子
    unsigned int seed = Catch::getSeed();
    // seed = 1343091126;
    unsigned int cluster_seed[4] = {seed + 1, seed + 2, seed + 3, seed + 4};

    // 目标数量
    std::vector<int> target_num = {10, 10, 10, 10};
    int target_num_sum = std::accumulate(target_num.begin(), target_num.end(), 0);
    int cluster_num = 200;

    // 目标创建
    auto points1 = generate_gaussian_points(target_num[0], 0, 250, 10, 150, 10, 100.0, 50.0, seed++);
    auto points2 = generate_gaussian_points(target_num[1], 0, 180, 10, 280, 10, 100.0, 50.0, seed++);
    auto points3 = generate_gaussian_points(target_num[2], 0, 180, 10, 150, 10, 100.0, 50.0, seed++);
    auto points4 = generate_gaussian_points(target_num[3], 0, 250, 10, 280, 10, 100.0, 50.0, seed++);
    auto points_cluster = std::vector<TrackPoint>(cluster_num);

    // 输入和输出
    std::vector<TrackPoint> points_all; // 所有点迹
    std::vector<std::array<TrackPoint, 4>> new_tracks(100);

    // 创建实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);                             // DEBUG容器
    TrackExtrapolator extrapolator(seed);                                   // 使用相同的随机数种子创建航迹外推器，保证噪声的一致性
    track_project::ManagementService track_manager(0.35, 1.65, 0.35, 1.65); // 恢复原来的范围

    // 绑定回调函数，显示航迹
    initiator.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                                 {
        LOG_INFO << "回调函数被调用，生成了 " << tracks.size() << " 条航迹";
        track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    //****************************************预处理数据***********************************************/
    track_manager.clear_all_command();
    points_all.insert(points_all.end(), points1.begin(), points1.end());
    points_all.insert(points_all.end(), points2.begin(), points2.end());
    points_all.insert(points_all.end(), points3.begin(), points3.end());
    points_all.insert(points_all.end(), points4.begin(), points4.end());
    for (size_t i = 0; i < points_all.size(); ++i)
    {
        auto [sigma_x, sigma_y] = extrapolator.getErrorDistribution(points_all[i].x, points_all[i].y); // 预热误差分布表格
        LOG_INFO << "目标点迹[" << i << "]: x=" << points_all[i].x << ", y=" << points_all[i].y
                 << ", longitude=" << points_all[i].longitude << ", latitude=" << points_all[i].latitude
                 << ", sigma_x=" << sigma_x << ", sigma_y=" << sigma_y;
    }

    points_cluster = generate_uniform_points(cluster_num, 0, 60, 170, 60, 170, 150, 100, cluster_seed[0]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    //*****************************************第一次处理数据***********************************************/
    LOG_INFO << "第一批次处理 - 时间片 0";
    track_manager.draw_point_command(points_all); // 绘制初始点迹
    ProcessStatus status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第二次处理数据***********************************************/
    LOG_INFO << "第二批次处理 - 时间片 1";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, 0, 60, 170, 60, 170, 150, 100, cluster_seed[1]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    track_manager.draw_point_command(points_all);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    //*****************************************第三次处理数据***********************************************/
    LOG_INFO << "第三批次处理 - 时间片 2";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, 0, 60, 170, 60, 170, 150, 100, cluster_seed[2]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    track_manager.draw_point_command(points_all);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    //*****************************************第四次处理数据***********************************************/
    LOG_INFO << "第四批次处理 - 时间片 3";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, 0, 60, 170, 60, 170, 150, 100, cluster_seed[3]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    // 绘制点迹
    track_manager.draw_point_command(points_all);

    // 画航迹，避免线程冲突
    sleep(1);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    tester.printHypothesisDistribution();

    //*****************************************统计数据***********************************************/
    points_all.resize(target_num_sum); // 只保留目标点迹，去掉杂波点迹
    //*****************************************统计数据***********************************************/
    // 初始化统计变量
    int true_track_count = 0;    // 正确航迹数
    int false_track_count = 0;   // 虚警航迹数
    int missed_target_count = 0; // 漏检目标数

    // 标记目标是否被匹配到
    std::vector<bool> target_matched(points_all.size(), false);

    // 遍历所有生成的航迹，统计虚警
    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &track = new_tracks[i];
        bool track_valid = false;

        for (size_t j = 0; j < points_all.size(); ++j)
        {
            size_t point_match = 0;
            for (size_t k = 0; k < 4; ++k)
            {
                // 作为AIS信息，我没改过，而且匀速直线，因此速度和航向不变，直接用来当匹配标识了
                if (track[k].sog < points_all[j].sog + 1e-4 && track[k].sog > points_all[j].sog - 1e-4 &&
                    track[k].cog < points_all[j].cog + 1e-4 && track[k].cog > points_all[j].cog - 1e-4)
                {
                    point_match++;
                }
            }

            if (point_match > 2) // 认为匹配成功，至少有3个点的速度和航向都匹配
            {
                track_valid = true;
                target_matched[j] = true; // 标记该目标已被匹配
                break;
            }
        }

        if (!track_valid)
        {
            false_track_count++;
            LOG_INFO << "第" << i << "条航迹为假航迹";
        }
        else
        {
            true_track_count++;
        }
    }

    // 统计漏检目标数
    for (size_t j = 0; j < points_all.size(); ++j)
    {
        if (!target_matched[j])
        {
            missed_target_count++;
            LOG_INFO << "第" << j << "个目标未被检测到";
        }
    }

    // 计算统计指标
    int total_targets = points_all.size();
    int total_tracks = new_tracks.size();

    double detection_rate = (total_targets - missed_target_count) * 100.0 / total_targets; // 检测率
    double false_alarm_rate = false_track_count * 100.0 / total_tracks;                    // 虚警率
    double missed_rate = missed_target_count * 100.0 / total_targets;                      // 漏检率

    // 输出统计结果
    LOG_INFO << "========== 统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << true_track_count;
    LOG_INFO << "虚警航迹数: " << false_track_count;
    LOG_INFO << "漏检目标数: " << missed_target_count;
    LOG_INFO << "检测率: " << detection_rate << "%";
    LOG_INFO << "虚警率: " << false_alarm_rate << "%";
    LOG_INFO << "漏检率: " << missed_rate << "%";

    // 保持窗口显示
    while (g_running)
    {
        wait_seconds(5);
    }
}

TEST_CASE("群目标测试", "[FunctionalityCheck][group_track]")
{
    std::signal(SIGINT, signal_handler);

    // 随机数种子
    unsigned int seed = Catch::getSeed();
    seed = 4053434269;

    // 群目标：4个同向同速、位置密集的目标
    std::vector<std::array<double, 4>> params = {
        {90.0, 10.0, 200, 10}, // 目标1
        {94.0, 10.0, 200, 10}, // 目标2，相距0.7km
        {98.0, 10.0, 200, 10}, // 目标3，相距0.5km
        {102.0, 10.0, 200, 10} // 目标4，相距0.7km
    };
    int target_num_sum = params.size();

    // 生成第一帧点迹
    auto points_all = generate_target_points_xyv(0, params);
    points_all[0].sog = 1.0;
    points_all[1].sog = 2.0;
    points_all[2].sog = 3.0;
    points_all[3].sog = 4.0;
    std::vector<std::array<TrackPoint, 4>> new_tracks(100);

    // 创建实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);
    TrackExtrapolator extrapolator(seed);
    track_project::ManagementService track_manager(0.60, 1.40, 0.08, 0.16);

    // 绑定回调函数，显示航迹
    initiator.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                                 {
        LOG_INFO << "回调函数被调用，生成了 " << tracks.size() << " 条航迹";
        track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    //****************************************预处理数据***********************************************/
    track_manager.clear_all_command();

    for (size_t i = 0; i < points_all.size(); ++i)
    {
        auto [sigma_x, sigma_y] = extrapolator.getErrorDistribution(points_all[i].x, points_all[i].y);
        // LOG_INFO << "目标点迹[" << i << "]: x=" << points_all[i].x << ", y=" << points_all[i].y
        //          << ", sog=" << points_all[i].sog << ", cog=" << points_all[i].cog
        //          << ", sigma_x=" << sigma_x << ", sigma_y=" << sigma_y;
    }

    //*****************************************第一次处理数据***********************************************/
    LOG_INFO << "第一批次处理 - 时间片 0";
    track_manager.draw_point_command(points_all);
    ProcessStatus status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    tester.printHypothesisDistribution();

    //*****************************************第二次处理数据***********************************************/
    LOG_INFO << "第二批次处理 - 时间片 1";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    track_manager.draw_point_command(points_all);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    tester.printHypothesisDistribution();

    //*****************************************第三次处理数据***********************************************/
    LOG_INFO << "第三批次处理 - 时间片 2";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    track_manager.draw_point_command(points_all);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    tester.printHypothesisDistribution();

    //*****************************************第四次处理数据***********************************************/
    LOG_INFO << "第四批次处理 - 时间片 3";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    track_manager.draw_point_command(points_all);
    sleep(1);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);
    tester.printHypothesisDistribution();

    //*****************************************统计数据***********************************************/
    int true_track_count = 0;
    int false_track_count = 0;
    int missed_target_count = 0;
    std::vector<bool> target_matched(points_all.size(), false);

    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &track = new_tracks[i];
        bool track_valid = false;

        for (size_t j = 0; j < points_all.size(); ++j)
        {
            size_t point_match = 0;
            for (size_t k = 0; k < 4; ++k)
            {
                if (std::abs(track[k].sog - points_all[j].sog) < 1e-4) // sog已经改成标签了，直接用来当匹配标识了
                {
                    point_match++;
                }
            }

            if (point_match > 2)
            {
                track_valid = true;
                target_matched[j] = true;
                break;
            }
        }

        if (!track_valid)
        {
            false_track_count++;
            LOG_INFO << "第" << i << "条航迹为假航迹";
        }
        else
        {
            true_track_count++;
        }
    }

    for (size_t j = 0; j < points_all.size(); ++j)
    {
        if (!target_matched[j])
        {
            missed_target_count++;
            LOG_INFO << "第" << j << "个目标未被检测到";
        }
    }

    // 遍历四条航迹的cog与点迹的cog对比，看误差是多少
    double standard_cog = points_all[0].cog; // COG是一致的，毕竟群目标
    double calculated_cog = 0.0;
    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &label_cog = new_tracks[i][3].cog;
        double cog_error = std::abs(label_cog - standard_cog);
        calculated_cog += cog_error;
    }

    int total_targets = points_all.size();
    int total_tracks = new_tracks.size();
    double detection_rate = (total_targets - missed_target_count) * 100.0 / total_targets;
    double false_alarm_rate = false_track_count * 100.0 / total_tracks;

    LOG_INFO << "========== 群目标统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << true_track_count;
    LOG_INFO << "虚警航迹数: " << false_track_count;
    LOG_INFO << "漏检目标数: " << missed_target_count;
    LOG_INFO << "检测率: " << detection_rate << "%";
    LOG_INFO << "虚警率: " << false_alarm_rate << "%";
    LOG_INFO << "航向一致性误差（平均）: " << calculated_cog / true_track_count << " deg";

    while (g_running)
    {
        wait_seconds(5);
    }
}

#include <chrono> // 添加时间相关头文件

TEST_CASE("多目标测试", "[Benchmark][multi_track]")
{
    std::signal(SIGINT, signal_handler); // 添加信号处理

    // 随机数种子
    unsigned int seed = Catch::getSeed();
    // seed = 1343091126;

    // 由seed,获取四个随机数用于生成杂波
    unsigned int cluster_seed[4] = {seed + 1, seed + 2, seed + 3, seed + 4};

    // 目标数量
    std::vector<int> target_num = {10, 10, 10, 10};
    int target_num_sum = std::accumulate(target_num.begin(), target_num.end(), 0);
    int cluster_num = 1960;

    // 目标创建
    auto points1 = generate_gaussian_points(target_num[0], 0, 280, 10, 280, 10, 100.0, 50.0, seed++);
    auto points2 = generate_gaussian_points(target_num[1], 0, 280, 10, 150, 10, 100.0, 50.0, seed++);
    auto points3 = generate_gaussian_points(target_num[2], 0, 150, 10, 280, 10, 100.0, 50.0, seed++);
    auto points4 = generate_gaussian_points(target_num[3], 0, 150, 10, 150, 10, 100.0, 50.0, seed++);
    auto points_cluster = std::vector<TrackPoint>(cluster_num);
    std::vector<double> cluster_param = {
        0,   // 统一时间戳
        00, // 最小距离
        400, // 最大距离
        00, // 最小距离
        400, // 最大距离
        00, // 速度中心
        50,  // 速度均值
    };

    // 输入和输出
    std::vector<TrackPoint> points_all; // 所有点迹
    std::vector<std::array<TrackPoint, 4>> new_tracks(100);

    // 创建实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator); // DEBUG容器
    TrackExtrapolator extrapolator(seed);       // 使用相同的随机数种子创建航迹外推器，保证噪声的一致性

    // 空回调加速
    initiator.set_track_callback([&](const std::vector<std::array<TrackPoint, 4>> &tracks) {});

    //****************************************预处理数据***********************************************/
    points_all.insert(points_all.end(), points1.begin(), points1.end());
    points_all.insert(points_all.end(), points2.begin(), points2.end());
    points_all.insert(points_all.end(), points3.begin(), points3.end());
    points_all.insert(points_all.end(), points4.begin(), points4.end());

    double error_mean = 0.0;
    for (size_t i = 0; i < points_all.size(); ++i)
    {
        auto [sigma_x, sigma_y] = extrapolator.getErrorDistribution(points_all[i].x, points_all[i].y); // 预热误差分布表格
        // LOG_INFO << "目标点迹[" << i << "]: x=" << points_all[i].x << ", y=" << points_all[i].y
        //          << ", longitude=" << points_all[i].longitude << ", latitude=" << points_all[i].latitude;
        error_mean += sqrt(sigma_x * sigma_x + sigma_y * sigma_y);
    }
    LOG_INFO << "ATTANTION!!!  平均位置误差: " << error_mean / points_all.size() << " km";

    points_cluster = generate_uniform_points(cluster_num, cluster_param[0], cluster_param[1], cluster_param[2],
                                             cluster_param[3], cluster_param[4], cluster_param[5], cluster_param[6], cluster_seed[0]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    // 存储每次处理的耗时（微秒）
    std::vector<long long> process_times_us;

    //*****************************************第一次处理数据***********************************************/
    LOG_INFO << "第一批次处理 - 时间片 0";
    auto start = std::chrono::high_resolution_clock::now();
    ProcessStatus status = initiator.process(points_all, new_tracks);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    process_times_us.push_back(duration_us);
    // LOG_INFO << "第一批次处理耗时: " << duration_us << " 微秒";
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第二次处理数据***********************************************/
    LOG_INFO << "第二批次处理 - 时间片 1";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, cluster_param[0], cluster_param[1], cluster_param[2],
                                             cluster_param[3], cluster_param[4], cluster_param[5], cluster_param[6], cluster_seed[1]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    start = std::chrono::high_resolution_clock::now();
    status = initiator.process(points_all, new_tracks);
    end = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    process_times_us.push_back(duration_us);
    // LOG_INFO << "第二批次处理耗时: " << duration_us << " 微秒";
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第三次处理数据***********************************************/
    LOG_INFO << "第三批次处理 - 时间片 2";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, cluster_param[0], cluster_param[1], cluster_param[2],
                                             cluster_param[3], cluster_param[4], cluster_param[5], cluster_param[6], cluster_seed[2]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    start = std::chrono::high_resolution_clock::now();
    status = initiator.process(points_all, new_tracks);
    end = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    process_times_us.push_back(duration_us);
    // LOG_INFO << "第三批次处理耗时: " << duration_us << " 微秒";
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第四次处理数据***********************************************/
    LOG_INFO << "第四批次处理 - 时间片 3";
    // 更新点迹位置（带噪声）
    points_all.resize(target_num_sum);
    extrapolator.update(points_all, TIME_INTERVAL_S);
    // 重置杂波点
    points_cluster = generate_uniform_points(cluster_num, cluster_param[0], cluster_param[1], cluster_param[2],
                                             cluster_param[3], cluster_param[4], cluster_param[5], cluster_param[6], cluster_seed[3]);
    points_all.insert(points_all.end(), points_cluster.begin(), points_cluster.end());

    // 绘制点迹

    // 画航迹，避免线程冲突
    start = std::chrono::high_resolution_clock::now();
    status = initiator.process(points_all, new_tracks);
    end = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    process_times_us.push_back(duration_us);
    // LOG_INFO << "第四批次处理耗时: " << duration_us << " 微秒";
    REQUIRE(status == ProcessStatus::SUCCESS);

    // 输出所有批次的统计信息
    LOG_INFO << "========== 处理时间统计 ==========";
    long long total_time_us = 0;
    for (size_t i = 0; i < process_times_us.size(); ++i)
    {
        // LOG_INFO << "批次 " << i + 1 << " 耗时: " << process_times_us[i] << " 微秒";
        total_time_us += process_times_us[i];
    }
    LOG_INFO << "总耗时: " << total_time_us << " 微秒 (" << total_time_us / 1000.0 << " 毫秒)";
    LOG_INFO << "平均耗时: " << total_time_us / process_times_us.size() << " 微秒";

    //*****************************************统计数据***********************************************/
    points_all.resize(target_num_sum); // 只保留目标点迹，去掉杂波点迹
    //*****************************************统计数据***********************************************/
    // 初始化统计变量
    int true_track_count = 0;    // 正确航迹数
    int false_track_count = 0;   // 虚警航迹数
    int missed_target_count = 0; // 漏检目标数

    // 标记目标是否被匹配到
    std::vector<bool> target_matched(points_all.size(), false);

    // 遍历所有生成的航迹，统计虚警
    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &track = new_tracks[i];
        bool track_valid = false;

        for (size_t j = 0; j < points_all.size(); ++j)
        {
            size_t point_match = 0;
            for (size_t k = 0; k < 4; ++k)
            {
                // 作为AIS信息，我没改过，而且匀速直线，因此速度和航向不变，直接用来当匹配标识了
                if (track[k].sog < points_all[j].sog + 1e-4 && track[k].sog > points_all[j].sog - 1e-4 &&
                    track[k].cog < points_all[j].cog + 1e-4 && track[k].cog > points_all[j].cog - 1e-4)
                {
                    point_match++;
                }
            }

            if (point_match > 2) // 认为匹配成功，至少有3个点的速度和航向都匹配
            {
                track_valid = true;
                target_matched[j] = true; // 标记该目标已被匹配
                break;
            }
        }

        if (!track_valid)
        {
            false_track_count++;
            // LOG_INFO << "第" << i << "条航迹为假航迹";
        }
        else
        {
            true_track_count++;
        }
    }

    // 统计漏检目标数
    for (size_t j = 0; j < points_all.size(); ++j)
    {
        if (!target_matched[j])
        {
            missed_target_count++;
            // LOG_INFO << "第" << j << "个目标未被检测到";
        }
    }

    // 计算统计指标
    int total_targets = points_all.size();
    int total_tracks = new_tracks.size();

    double detection_rate = (total_targets - missed_target_count) * 100.0 / total_targets; // 检测率
    double false_alarm_rate = false_track_count * 100.0 / total_tracks;                    // 虚警率
    double missed_rate = missed_target_count * 100.0 / total_targets;                      // 漏检率

    // 输出统计结果
    LOG_INFO << "========== 统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << true_track_count;
    LOG_INFO << "虚警航迹数: " << false_track_count;
    LOG_INFO << "漏检目标数: " << missed_target_count;
    LOG_INFO << "检测率: " << detection_rate << " %";
    LOG_INFO << "虚警率: " << false_alarm_rate << " %";
    LOG_INFO << "漏检率: " << missed_rate << " %";
    LOG_INFO << "SEED: " << seed;
}

TEST_CASE("群目标测试", "[Benchmark][group_track]")
{

    // 随机数种子
    unsigned int seed = Catch::getSeed();
    // seed = 42;

    // 群目标：4个同向同速、位置密集的目标
    std::vector<std::array<double, 4>> params = {
        {90.0, 10.0, 200, 10},  // 目标1
        {95.0, 10.0, 200, 10},  // 目标2，相距0.7km
        {100.0, 10.0, 200, 10}, // 目标3，相距0.5km
        {105.0, 10.0, 200, 10}  // 目标4，相距0.7km
    };
    int target_num_sum = params.size();

    // 生成第一帧点迹
    auto points_all = generate_target_points_xyv(0, params);
    points_all[0].sog = 1.0;
    points_all[1].sog = 2.0;
    points_all[2].sog = 3.0;
    points_all[3].sog = 4.0;
    std::vector<std::array<TrackPoint, 4>> new_tracks(100);

    // 创建实例
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);
    TrackExtrapolator extrapolator(seed);

    // 绑定回调函数，显示航迹
    initiator.set_track_callback([&](const std::vector<std::array<TrackPoint, 4>> &tracks)
                                 { LOG_INFO << "回调函数被调用，生成了 " << tracks.size() << " 条航迹"; });

    //****************************************预处理数据***********************************************/

    for (size_t i = 0; i < points_all.size(); ++i)
    {
        auto [sigma_x, sigma_y] = extrapolator.getErrorDistribution(points_all[i].x, points_all[i].y);
        // LOG_INFO << "目标点迹[" << i << "]: x=" << points_all[i].x << ", y=" << points_all[i].y
        //          << ", sog=" << points_all[i].sog << ", cog=" << points_all[i].cog
        //          << ", sigma_x=" << sigma_x << ", sigma_y=" << sigma_y;
    }

    //*****************************************第一次处理数据***********************************************/
    LOG_INFO << "第一批次处理 - 时间片 0";
    ProcessStatus status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第二次处理数据***********************************************/
    LOG_INFO << "第二批次处理 - 时间片 1";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第三次处理数据***********************************************/
    LOG_INFO << "第三批次处理 - 时间片 2";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    std::vector<TrackPoint> lose_group = {points_all[0], points_all[2], points_all[3]}; // 模拟群目标中的一个目标突然消失
    LOG_INFO << "模拟群目标中的一个目标突然消失，剩余点迹数: " << lose_group.size();
    status = initiator.process(lose_group, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************第四次处理数据***********************************************/
    LOG_INFO << "第四批次处理 - 时间片 3";
    extrapolator.update(points_all, TIME_INTERVAL_S);
    status = initiator.process(points_all, new_tracks);
    REQUIRE(status == ProcessStatus::SUCCESS);

    //*****************************************统计数据***********************************************/
    int true_track_count = 0;
    int false_track_count = 0;
    int missed_target_count = 0;
    std::vector<bool> target_matched(points_all.size(), false);

    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &track = new_tracks[i];
        bool track_valid = false;

        for (size_t j = 0; j < points_all.size(); ++j)
        {
            size_t point_match = 0;
            for (size_t k = 0; k < 4; ++k)
            {
                if (std::abs(track[k].sog - points_all[j].sog) < 1e-4) // sog已经改成标签了，直接用来当匹配标识了
                {
                    point_match++;
                }
            }

            if (point_match > 2)
            {
                track_valid = true;
                target_matched[j] = true;
                break;
            }
        }

        if (!track_valid)
        {
            false_track_count++;
            LOG_INFO << "第" << i << "条航迹为假航迹";
        }
        else
        {
            true_track_count++;
        }
    }

    for (size_t j = 0; j < points_all.size(); ++j)
    {
        if (!target_matched[j])
        {
            missed_target_count++;
            LOG_INFO << "第" << j << "个目标未被检测到";
        }
    }

    // 遍历四条航迹的cog与点迹的cog对比，看误差是多少
    double standard_cog = points_all[0].cog; // COG是一致的，毕竟群目标
    double calculated_cog = 0.0;
    for (size_t i = 0; i < new_tracks.size(); ++i)
    {
        const auto &label_cog = new_tracks[i][3].cog;
        double cog_error = std::abs(label_cog - standard_cog);
        calculated_cog += cog_error;
    }

    int total_targets = points_all.size();
    int total_tracks = new_tracks.size();
    double detection_rate = (total_targets - missed_target_count) * 100.0 / total_targets;
    double false_alarm_rate = false_track_count * 100.0 / total_tracks;

    LOG_INFO << "========== 群目标统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << true_track_count;
    LOG_INFO << "虚警航迹数: " << false_track_count;
    LOG_INFO << "漏检目标数: " << missed_target_count;
    LOG_INFO << "检测率: " << detection_rate << " %";
    LOG_INFO << "虚警率: " << false_alarm_rate << " %";
    LOG_INFO << "航向一致性误差（平均）: " << calculated_cog / true_track_count << " deg";
}