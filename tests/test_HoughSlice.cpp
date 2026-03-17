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
                                            100.0, 30.0,
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
    size_t cluster_point_num = 40; // 每批杂波点数
    size_t target_num = 10;        // 目标点数
    size_t num_rounds = 10;        // 测试轮数

    std::vector<double> processing_times; // 存储每轮处理时间
    std::vector<size_t> track_counts;     // 存储每轮产生的航迹数

    // 新增统计量
    std::vector<double> detection_rates;     // 检测率
    std::vector<double> false_alarm_rates;   // 虚警率
    std::vector<double> doppler_error_rates; // Doppler误差

    for (size_t round = 0; round < num_rounds; ++round)
    {
        // 每轮用不同的种子
        unsigned int base_seed = Catch::getSeed() + round * 100;

        // 生成目标点迹（10个目标）- 保存真实值用于对比
        std::vector<TrackPoint> ground_truth_targets; // 真实目标状态
        auto target_points = generate_gaussian_points(
            target_num, 0,
            150, 10, 150, 10,
            100.0, 50.0,
            base_seed);

        // 保存目标的真实Doppler
        for (const auto &p : target_points)
        {
            TrackPoint truth = p;
            // 计算真实Doppler（径向速度）
            double range = std::sqrt(p.x * p.x + p.y * p.y);
            if (range > 1e-6)
            {
                double los_x = -p.x / range;
                double los_y = -p.y / range;
                truth.doppler = p.vx * los_x + p.vy * los_y;
            }
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
        track_project::ManagementService track_manager(0.8, 1.8, 0.8, 1.8);

        // 存储本轮检测到的目标
        std::set<int> detected_target_indices;
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
        all_points.insert(all_points.end(), target_points.begin(), target_points.end());
        all_points.insert(all_points.end(),
                          clutter_points.begin(),
                          clutter_points.begin() + cluster_point_num);

        track_manager.draw_point_command(all_points);
        alg.process(all_points, new_tracks);

        // 第2轮处理
        for (auto &p : target_points)
        {
            point_update_cv(p, TIME_INTERVAL_S);
        }
        all_points.clear();
        all_points.insert(all_points.end(), target_points.begin(), target_points.end());
        all_points.insert(all_points.end(),
                          clutter_points.begin() + cluster_point_num,
                          clutter_points.begin() + 2 * cluster_point_num);

        track_manager.draw_point_command(all_points);
        alg.process(all_points, new_tracks);

        // 第3轮处理
        for (auto &p : target_points)
        {
            point_update_cv(p, TIME_INTERVAL_S);
        }
        all_points.clear();
        all_points.insert(all_points.end(), target_points.begin(), target_points.end());
        all_points.insert(all_points.end(),
                          clutter_points.begin() + 2 * cluster_point_num,
                          clutter_points.begin() + 3 * cluster_point_num);

        track_manager.draw_point_command(all_points);
        alg.process(all_points, new_tracks);

        // 第4轮处理
        for (auto &p : target_points)
        {
            point_update_cv(p, TIME_INTERVAL_S);
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

        // ========== 计算检测率和虚警率（基于Doppler）==========
        const double POSITION_TOLERANCE_KM = 10.0; // 1公里匹配阈值
        const double DOPPLER_TOLERANCE = 2.0;      // 2 m/s 多普勒容差

        std::set<int> matched_targets;
        int false_alarms = 0;

        // 统计Doppler误差
        double total_doppler_error = 0.0;
        int valid_matches = 0;

        for (const auto &track : detected_tracks)
        {
            bool matched = false;
            const TrackPoint &last_point = track[3]; // 用最新的点

            // 找最近的真实目标
            double min_dist = POSITION_TOLERANCE_KM;
            int best_target_idx = -1;

            for (size_t i = 0; i < ground_truth_targets.size(); ++i)
            {
                double dx = last_point.x - ground_truth_targets[i].x;
                double dy = last_point.y - ground_truth_targets[i].y;
                double dist = std::sqrt(dx * dx + dy * dy);

                if (dist < min_dist)
                {
                    min_dist = dist;
                    best_target_idx = i;
                }
            }

            if (best_target_idx >= 0)
            {
                // 找到了匹配的真实目标
                matched = true;
                matched_targets.insert(best_target_idx);

                // 计算Doppler误差
                double doppler_error = std::abs(last_point.doppler - ground_truth_targets[best_target_idx].doppler);
                total_doppler_error += doppler_error;
                valid_matches++;

                // 如果Doppler误差太大，也算虚警
                if (doppler_error > DOPPLER_TOLERANCE)
                {
                    false_alarms++;
                }
            }
            else
            {
                false_alarms++;
            }
        }

        // 计算本轮的检测率、虚警率和Doppler误差
        double detection_rate = static_cast<double>(matched_targets.size()) / target_num;
        double false_alarm_rate = static_cast<double>(false_alarms) / track_count;
        double doppler_error = (valid_matches > 0) ? (total_doppler_error / valid_matches) : 0.0;

        detection_rates.push_back(detection_rate);
        false_alarm_rates.push_back(false_alarm_rate);
        doppler_error_rates.push_back(doppler_error);

        LOG_INFO << "第 " << round + 1 << " 轮: 处理时间 = "
                 << std::fixed << std::setprecision(2) << processing_times.back()
                 << " ms, 航迹数 = " << track_counts.back()
                 << ", 检测率 = " << std::setprecision(1) << detection_rate * 100 << "%"
                 << ", 虚警率 = " << std::setprecision(1) << false_alarm_rate * 100 << "%"
                 << ", Doppler误差 = " << std::setprecision(2) << doppler_error << " m/s";
    }

    // 计算统计信息
    auto calc_mean = [](const std::vector<double> &v)
    {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };

    auto calc_variance = [](const std::vector<double> &v, double mean)
    {
        double sum = 0.0;
        for (double x : v)
        {
            sum += (x - mean) * (x - mean);
        }
        return sum / v.size();
    };

    auto calc_stddev = [](double variance)
    {
        return std::sqrt(variance);
    };

    // 输出统计结果
    double time_mean = calc_mean(processing_times);
    double time_var = calc_variance(processing_times, time_mean);
    double time_std = calc_stddev(time_var);

    double track_mean = calc_mean(std::vector<double>(track_counts.begin(), track_counts.end()));
    double track_var = calc_variance(std::vector<double>(track_counts.begin(), track_counts.end()), track_mean);
    double track_std = calc_stddev(track_var);

    LOG_INFO << "\n========== 抗杂波能力 Benchmark 结果 ==========";
    LOG_INFO << "测试配置: " << num_rounds << " 轮, 每轮 " << target_num
             << " 目标, " << cluster_point_num << " 杂波/帧";

    LOG_INFO << "\n处理时间统计 (ms):";
    LOG_INFO << "  平均值: " << std::fixed << std::setprecision(2) << time_mean << " ms";
    LOG_INFO << "  标准差: " << std::fixed << std::setprecision(2) << time_std << " ms";

    LOG_INFO << "\n检测性能统计:";
    double det_mean = calc_mean(detection_rates);
    double fa_mean = calc_mean(false_alarm_rates);
    double dop_mean = calc_mean(doppler_error_rates);

    LOG_INFO << "  平均检测率: " << std::fixed << std::setprecision(1) << det_mean * 100 << "%";
    LOG_INFO << "  平均虚警率: " << std::fixed << std::setprecision(2) << fa_mean * 100 << "%";
    LOG_INFO << "  平均Doppler误差: " << std::fixed << std::setprecision(2) << dop_mean << " m/s";

    LOG_INFO << "\n航迹数统计:";
    LOG_INFO << "  平均值: " << std::fixed << std::setprecision(2) << track_mean;
    LOG_INFO << "  标准差: " << std::fixed << std::setprecision(2) << track_std;
    LOG_INFO << "  最小值: " << *std::min_element(track_counts.begin(), track_counts.end());
    LOG_INFO << "  最大值: " << *std::max_element(track_counts.begin(), track_counts.end());
    LOG_INFO << "================================================";

    // // 性能断言
    // CHECK(time_mean < 100.0);
    // CHECK(det_mean > 0.8); // 检测率 > 80%
    // CHECK(fa_mean < 0.3);  // 虚警率 < 30%
    // CHECK(dop_mean < 2.0); // Doppler误差 < 2 m/s
}

// TEST_CASE("聚类压测（废弃，另一个函数已经被删除了，log有）", "[clust_gen]")
// {
//     // 随机数种子
//     unsigned int seed = Catch::getSeed();
//     std::mt19937 gen(seed);

//     // 参数定义
//     const size_t N_POINTS = 1000;
//     const size_t N_CLUSTERS_5 = 5;
//     const size_t POINTS_PER_CLUSTER = N_POINTS / N_CLUSTERS_5;
//     const double x_min = -400.0, x_max = 400.0;
//     const double y_min = -400.0, y_max = 400.0;
//     const double v_mean = 10.0, v_stddev = 5.0;
//     const int x_cell = 1 + static_cast<int>((x_max - x_min) / HOUGHSLICE_CLUSTER_RADIUS_KM);
//     const int y_cell = 1 + static_cast<int>((y_max - y_min) / HOUGHSLICE_CLUSTER_RADIUS_KM);

//     // 生成数据
//     std::vector<TrackPoint> uniform_point, gauss_point, rayleigh_point;
//     std::vector<std::pair<double, double>> centers;
//     std::set<std::pair<int, int>> used_cells;

//     // 生成噪声点迹
//     auto noise_point = generate_uniform_points(N_POINTS, 0, x_min, x_max, y_min, y_max, v_mean, v_stddev, seed);

//     // 生成聚类中心 - 使用普通随机数生成器
//     std::uniform_int_distribution<> dist_x(0, x_cell - 1); // 注意边界
//     std::uniform_int_distribution<> dist_y(0, y_cell - 1);

//     // 抽调一个聚类中心用于构造新航迹
//     used_cells.insert({0, 0});
//     auto single_clust_point = generate_gaussian_points(
//         POINTS_PER_CLUSTER, TIME_INTERVAL_S,
//         0, HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//         0, HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//         v_mean, v_stddev, seed);

//     while (centers.size() < N_CLUSTERS_5)
//     {
//         int cell_x = dist_x(gen);
//         int cell_y = dist_y(gen);

//         if (used_cells.find({cell_x, cell_y}) == used_cells.end())
//         {
//             used_cells.insert({cell_x, cell_y});
//             double center_x = x_min + (cell_x + 0.5) * HOUGHSLICE_CLUSTER_RADIUS_KM;
//             double center_y = y_min + (cell_y + 0.5) * HOUGHSLICE_CLUSTER_RADIUS_KM;
//             centers.emplace_back(center_x, center_y);
//         }
//     }

//     LOG_INFO << "生成聚类中心点完成，聚类中心位置" << centers;

//     // 为每个聚类中心生成三种分布的点
//     for (const auto &center : centers)
//     {
//         // 均匀分布点
//         auto temp_uniform = generate_uniform_points(
//             POINTS_PER_CLUSTER, 0,
//             center.first - HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//             center.first + HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//             center.second - HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//             center.second + HOUGHSLICE_CLUSTER_RADIUS_KM / 2,
//             v_mean, v_stddev, seed);
//         uniform_point.insert(uniform_point.end(), temp_uniform.begin(), temp_uniform.end());

//         // 高斯分布点
//         auto temp_gauss = generate_gaussian_points(
//             POINTS_PER_CLUSTER, 0,
//             center.first, HOUGHSLICE_CLUSTER_RADIUS_KM / 4,
//             center.second, HOUGHSLICE_CLUSTER_RADIUS_KM / 4,
//             v_mean, v_stddev, seed);
//         gauss_point.insert(gauss_point.end(), temp_gauss.begin(), temp_gauss.end());

//         // 瑞利分布点
//         auto temp_rayleigh = generate_rayleigh_points(
//             POINTS_PER_CLUSTER, 0,
//             center.first, center.second,
//             HOUGHSLICE_CLUSTER_RADIUS_KM / 4, 1.0,
//             v_mean, v_stddev, seed);
//         rayleigh_point.insert(rayleigh_point.end(), temp_rayleigh.begin(), temp_rayleigh.end());
//     }

//     // 存入Accessor构造第一批聚类
//     SECTION("完全随机情况下的点迹更新测试")
//     {
//         HoughSlice alg;
//         test_HoughSlice::run_cluster_process(alg, noise_point); // 预先放入一批数据

//         // 更新所有点迹
//         for (auto &p : noise_point)
//         {
//             point_update(p, TIME_INTERVAL_S);
//         }

//         BENCHMARK("cluster_process")
//         {
//             test_HoughSlice::run_cluster_process(alg, noise_point);
//         };
//     }

//     SECTION("五个均匀分布聚类测试")
//     {
//         HoughSlice alg;
//         test_HoughSlice::run_cluster_process(alg, uniform_point); // 预先放入一批数据

//         // 更新所有点迹
//         for (auto &p : uniform_point)
//         {
//             point_update(p, TIME_INTERVAL_S);
//         }

//         // 添加噪声点和聚类点
//         uniform_point.insert(uniform_point.end(), noise_point.begin(), noise_point.end());
//         uniform_point.insert(uniform_point.end(), single_clust_point.begin(), single_clust_point.end());

//         BENCHMARK("cluster_process")
//         {
//             test_HoughSlice::run_cluster_process(alg, uniform_point);
//         };
//     }

//     SECTION("五个高斯分布聚类测试")
//     {
//         HoughSlice alg;
//         test_HoughSlice::run_cluster_process(alg, gauss_point); // 预先放入一批数据

//         // 更新所有点迹
//         for (auto &p : gauss_point)
//         {
//             point_update(p, TIME_INTERVAL_S);
//         }

//         // 添加噪声点和聚类点
//         gauss_point.insert(gauss_point.end(), noise_point.begin(), noise_point.end());
//         gauss_point.insert(gauss_point.end(), single_clust_point.begin(), single_clust_point.end());

//         BENCHMARK("cluster_process")
//         {
//             test_HoughSlice::run_cluster_process(alg, gauss_point);
//         };
//     }

//     SECTION("五个瑞利分布聚类测试")
//     {
//         HoughSlice alg;
//         test_HoughSlice::run_cluster_process(alg, rayleigh_point); // 预先放入一批数据

//         // 更新所有点迹
//         for (auto &p : rayleigh_point)
//         {
//             point_update(p, TIME_INTERVAL_S);
//         }

//         // 添加噪声点和聚类点
//         rayleigh_point.insert(rayleigh_point.end(), noise_point.begin(), noise_point.end());
//         rayleigh_point.insert(rayleigh_point.end(), single_clust_point.begin(), single_clust_point.end());

//         BENCHMARK("cluster_process")
//         {
//             test_HoughSlice::run_cluster_process(alg, rayleigh_point);
//         };
//     }
// }