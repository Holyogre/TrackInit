#include <catch2/catch_all.hpp>
#include "../src/HoughSlice.hpp"
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

#include "Func_pointgen.hpp" //测试点迹生成以及更新函数
#include "Func_cout.hpp"     //一些<<重载，方便输出
#include "../utils/Logger.hpp"

// 检测时间间隔，单位秒
constexpr double TIME_INTERVAL_S = 10.0;

using track_project::TrackPoint;
using track_project::trackinit::HoughSlice;
using track_project::trackinit::HOUGHSLICE_CLUSTER_RADIUS_KM;

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

// 友元访问器，实现对私有聚类函数的调用
namespace track_project::trackinit
{
    struct test_HoughSlice
    {
        static void run_cluster_process(HoughSlice &alg, const std::vector<TrackPoint> &pts)
        {
            alg.process_cluster_generation(pts);
        }
    };
}

using track_project::trackinit::test_HoughSlice;

TEST_CASE("功能测试", "[FunctionalityCheck]")
{
    std::signal(SIGINT, signal_handler);

    // 随机数种子
    // unsigned int seed = Catch::getSeed();
    unsigned int seed = 42;

    // 生成高斯分布的点迹
    auto points1 = generate_gaussian_points(2, 0,
                                            150, 10,
                                            150, 10,
                                            200.0, 30.0,
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

    std::vector<std::array<TrackPoint, 4>> new_tracks;

    // 创建HoughSlice对象
    HoughSlice alg;

    // 创建航迹管理器
    track_project::ManagementService track_manager(0.3, 1.8, 0.3, 1.8); // 经纬度范围

    // 绑定回调函数，显示航迹
    alg.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                           { track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    // 首次处理点迹
    track_manager.clear_all_command(); // 清空状态

    //*****************************************第一次处理数据***********************************************/
    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************第二次处理数据***********************************************/
    // seed++;
    // 更新点迹位置并装在
    for (auto &p : points_all)
    {
        // point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
        LOG_DEBUG << "更新后点迹: " << p;
    }

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************第三次处理数据***********************************************/
    // 更新点迹位置并装在
    // seed++;
    for (auto &p : points_all)
    {
        // point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
        LOG_DEBUG << "更新后点迹: " << p;
    }

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************最终输出处理数据***********************************************/
    // 更新点迹位置并装在
    // seed++;
    for (auto &p : points_all)
    {
        // point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
        LOG_DEBUG << "更新后点迹: " << p;
    }

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    // std::this_thread::sleep_for(std::chrono::seconds(5));
    while (g_running)
    {
        wait_seconds(5); // 可被 CTRL+C 中断的等待
    }
}

TEST_CASE("抗杂波能力", "[ClutterResistance]")
{
    size_t cluster_point_num = 40;
    std::signal(SIGINT, signal_handler);

    // 随机数种子
    unsigned int seed = Catch::getSeed();
    // unsigned int seed = 42;

    // 生成高斯分布的点迹
    auto points1 = generate_gaussian_points(10, 0,
                                            150, 10,
                                            150, 10,
                                            100.0, 50.0,
                                            seed++);

    auto point_cluster = generate_gaussian_points(4 * cluster_point_num, 0,
                                                  150, 10,
                                                  150, 10,
                                                  100.0, 50.0,
                                                  seed++);

    std::vector<TrackPoint> points_all; // 所有点迹
    points_all.insert(points_all.end(), points1.begin(), points1.end());
    points_all.insert(points_all.end(), point_cluster.begin(), point_cluster.begin() + cluster_point_num);

    std::vector<std::array<TrackPoint, 4>> new_tracks;

    // 创建HoughSlice对象
    HoughSlice alg;

    // 创建航迹管理器
    track_project::ManagementService track_manager(0.8, 1.8, 0.8, 1.8); // 经纬度范围

    // 绑定回调函数，显示航迹
    alg.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                           { track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    // 首次处理点迹
    track_manager.clear_all_command(); // 清空状态

    //*****************************************第一次处理数据***********************************************/
    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************第二次处理数据***********************************************/
    seed++;
    // 更新点迹位置并装在
    for (auto &p : points_all)
    {
        point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        // point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
    }
    // 刷新杂波点，保持杂波点数量不变
    points_all.resize(points1.size());
    points_all.insert(points_all.end(), point_cluster.begin() + cluster_point_num, point_cluster.begin() + 2 * cluster_point_num);

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************第三次处理数据***********************************************/
    // 更新点迹位置并装在
    seed++;
    for (auto &p : points_all)
    {
        point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        // point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
    }
    points_all.resize(points1.size());
    points_all.insert(points_all.end(), point_cluster.begin() + 2 * cluster_point_num, point_cluster.begin() + 3 * cluster_point_num);

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    //*****************************************最终输出处理数据***********************************************/
    // 更新点迹位置并装在
    seed++;
    for (auto &p : points_all)
    {
        point_update_cv(p, TIME_INTERVAL_S); // 无噪声
        // point_update_cv_with_noise(p, TIME_INTERVAL_S, seed);
    }
    points_all.resize(points1.size());
    points_all.insert(points_all.end(), point_cluster.begin() + 3 * cluster_point_num, point_cluster.begin() + 4 * cluster_point_num);

    track_manager.draw_point_command(points_all); // 绘制点迹
    alg.process(points_all, new_tracks);

    // std::this_thread::sleep_for(std::chrono::seconds(5));
    while (g_running)
    {
        wait_seconds(5); // 可被 CTRL+C 中断的等待
    }
}

TEST_CASE("抗杂波能力 Benchmark", "[Benchmark]")
{
    std::signal(SIGINT, signal_handler);
    size_t cluster_point_num = 40; // 每批杂波点数
    size_t target_num = 10;        // 目标点数
    size_t num_rounds = 4;         // 测试轮数（4轮）

    std::vector<double> processing_times; // 存储每轮处理时间
    std::vector<size_t> track_counts;     // 存储每轮产生的航迹数

    // 新增统计量
    std::vector<double> detection_rates;   // 检测率
    std::vector<double> false_alarm_rates; // 虚警率
    std::vector<double> missed_rates;      // 漏检率

    // 每轮用不同的种子
    unsigned int base_seed = Catch::getSeed();
    base_seed = 3201435970; // 复现问题

    // 生成目标点迹（10个目标）- 保存真实值用于对比
    std::vector<TrackPoint> ground_truth_targets; // 真实目标状态（每个目标一个点，用于匹配标识）
    auto target_points = generate_gaussian_points(
        target_num, 0,
        150, 10, 150, 10,
        100.0, 50.0,
        base_seed);

    // 保存目标的真实SOG和COG（用于匹配）
    for (const auto &p : target_points)
    {
        TrackPoint truth = p;
        ground_truth_targets.push_back(truth);
    }

    // 生成杂波点（足够多，供4轮使用）
    auto clutter_points = generate_gaussian_points(
        4 * cluster_point_num, 0,
        150, 10, 150, 10,
        100.0, 50.0,
        base_seed + 1000);

    std::vector<TrackPoint> all_points;
    std::vector<std::array<TrackPoint, 4>> new_tracks;

    // 创建算法实例
    HoughSlice alg;
    track_project::ManagementService track_manager(1.1, 1.5, 1.1, 1.5);

    // 存储本轮检测到的航迹
    std::vector<std::array<TrackPoint, 4>> detected_tracks;

    // 计数回调 - 同时收集检测到的航迹
    std::atomic<size_t> track_count{0};
    alg.set_track_callback([&](const std::vector<std::array<TrackPoint, 4>> &tracks)
                           {
            track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks));
            track_count += tracks.size();
            
            // 保存检测到的航迹用于后续分析
            for (const auto& track : tracks) {
                detected_tracks.push_back(track);
            } });

    // 开始计时
    auto start_time = std::chrono::high_resolution_clock::now();

    // 第1轮处理
    track_manager.clear_all_command();
    all_points.clear();
    detected_tracks.clear(); // 清空上一轮的航迹
    track_count = 0;

    all_points.insert(all_points.end(), target_points.begin(), target_points.end());
    all_points.insert(all_points.end(),
                      clutter_points.begin(),
                      clutter_points.begin() + cluster_point_num);

    track_manager.draw_point_command(all_points);
    alg.process(all_points, new_tracks);

    // 第2轮处理
    for (auto &p : target_points)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, base_seed + 500);
    }
    all_points.clear();
    all_points.insert(all_points.end(), target_points.begin(), target_points.end());
    all_points.insert(all_points.end(),
                      clutter_points.begin() + cluster_point_num,
                      clutter_points.begin() + 2 * cluster_point_num);

    // track_manager.draw_point_command(all_points);
    alg.process(all_points, new_tracks);

    // 第3轮处理
    for (auto &p : target_points)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, base_seed + 1000);
    }
    all_points.clear();
    all_points.insert(all_points.end(), target_points.begin(), target_points.end());
    all_points.insert(all_points.end(),
                      clutter_points.begin() + 2 * cluster_point_num,
                      clutter_points.begin() + 3 * cluster_point_num);

    // track_manager.draw_point_command(all_points);
    alg.process(all_points, new_tracks);

    // 第4轮处理
    for (auto &p : target_points)
    {
        point_update_cv_with_noise(p, TIME_INTERVAL_S, base_seed + 2000); // 添加噪声参数
    }
    all_points.clear();
    all_points.insert(all_points.end(), target_points.begin(), target_points.end());
    all_points.insert(all_points.end(),
                      clutter_points.begin() + 3 * cluster_point_num,
                      clutter_points.begin() + 4 * cluster_point_num);

    track_manager.draw_point_command(all_points);
    alg.process(all_points, new_tracks);

    // 结束计时
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    processing_times.push_back(duration.count() / 1000.0);
    track_counts.push_back(track_count);

    // ========== 统计检测率和虚警率（基于 SOG/COG）==========
    const double COG_TOLERANCE_DEG = 5.0; // 航向容差（度）
    const double SOG_TOLERANCE_KN = 2.0;  // 速度容差（节）
    const int MIN_MATCH_POINTS = 3;       // 最少匹配点数

    // 统计变量
    std::set<int> matched_targets; // 已匹配的真实目标索引
    int false_alarms = 0;          // 虚警数
    int missed_targets = 0;        // 漏检数
    int total_valid_matches = 0;   // 有效匹配数

    // 标记每个真实目标是否被匹配（用于漏检统计）
    std::vector<bool> target_matched(ground_truth_targets.size(), false);

    // 遍历所有检测到的航迹，统计虚警
    for (size_t i = 0; i < detected_tracks.size(); ++i)
    {
        const auto &track = detected_tracks[i];
        bool track_matched = false;

        // 找到匹配的真实目标
        int matched_target_idx = -1;
        int max_match_count = 0;

        for (size_t j = 0; j < ground_truth_targets.size(); ++j)
        {
            const auto &target = ground_truth_targets[j];
            int match_count = 0;

            // 用航迹的所有点进行匹配
            for (size_t k = 0; k < 4 && k < track.size(); ++k)
            {
                // 检查速度和航向是否匹配
                if (std::abs(track[k].sog - target.sog) < SOG_TOLERANCE_KN &&
                    std::abs(track[k].cog - target.cog) < COG_TOLERANCE_DEG)
                {
                    match_count++;
                }
            }

            if (match_count > max_match_count)
            {
                max_match_count = match_count;
                matched_target_idx = j;
            }
        }

        // 判断是否匹配成功
        if (matched_target_idx >= 0 && max_match_count >= MIN_MATCH_POINTS)
        {
            track_matched = true;
            matched_targets.insert(matched_target_idx);
            target_matched[matched_target_idx] = true;
            total_valid_matches++;

            LOG_DEBUG << "航迹 " << i << " 匹配到目标 " << matched_target_idx
                      << "，匹配点数: " << max_match_count;
        }
        else
        {
            false_alarms++;
            LOG_INFO << "航迹 " << i << " 为虚警，最大匹配点数: " << max_match_count;
        }
    }

    // 统计漏检目标
    for (size_t j = 0; j < ground_truth_targets.size(); ++j)
    {
        if (!target_matched[j])
        {
            missed_targets++;
            LOG_INFO << "目标 " << j << " 未被检测到";
        }
    }

    // 计算统计指标
    int total_targets = ground_truth_targets.size();
    int total_tracks = detected_tracks.size();

    double detection_rate = (total_targets - missed_targets) * 100.0 / total_targets;
    double false_alarm_rate = (total_tracks > 0) ? false_alarms * 100.0 / total_tracks : 0.0;
    double missed_rate = missed_targets * 100.0 / total_targets;

    // 存储统计结果
    detection_rates.push_back(detection_rate);
    false_alarm_rates.push_back(false_alarm_rate);
    missed_rates.push_back(missed_rate);

    // 输出统计结果
    LOG_INFO << "========== 统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << (total_tracks - false_alarms);
    LOG_INFO << "虚警航迹数: " << false_alarms;
    LOG_INFO << "漏检目标数: " << missed_targets;
    LOG_INFO << "检测率: " << detection_rate << "%";
    LOG_INFO << "虚警率: " << false_alarm_rate << "%";
    LOG_INFO << "漏检率: " << missed_rate << "%";

    // 可选：输出平均值
    LOG_INFO << "========== 性能指标 ==========";
    LOG_INFO << "平均处理时间: " << processing_times.back() << " ms";
    LOG_INFO << "平均检测率: " << detection_rate << "%";
    LOG_INFO << "平均虚警率: " << false_alarm_rate << "%";

    // 断言（根据期望值调整）
    // REQUIRE(detection_rate > 80.0);   // 检测率应大于80%
    // REQUIRE(false_alarm_rate < 20.0); // 虚警率应小于20%
    while (g_running)
    {
        wait_seconds(5); // 可被 CTRL+C 中断的等待
    }
}
