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
    auto points = generate_gaussian_points(1, 0,
                                           0, 10,
                                           0, 10,
                                           100.0, 50.0,
                                           seed);
    std::vector<std::array<TrackPoint, 4>> new_tracks;

    // //debug测试
    // auto point=points[1];
    // points.clear();
    // points.push_back(point);

    // 创建HoughSlice对象
    HoughSlice alg;

    // 创建航迹管理器
    track_project::ManagementService track_manager(-0.5, 0.5, -0.5, 0.5); // 经纬度范围

    // 绑定回调函数，显示航迹
    alg.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
                           { track_manager.create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks)); });

    // 首次处理点迹
    track_manager.clear_all_command(); // 清空状态

    for (auto &p : points)
    {
        LOG_INFO << "初始点迹位置：" << p;
    }
    track_manager.draw_point_command(points); // 绘制点迹
    alg.process(points, new_tracks);

    // 更新点迹位置
    for (auto &p : points)
    {
        point_update(p, TIME_INTERVAL_S);
        LOG_INFO << "第一次更新后点迹位置：" << p;
    }
    track_manager.draw_point_command(points); // 绘制点迹
    alg.process(points, new_tracks);

    // 第三次更新和第二次本质上是一样的，所以其实没必要继续更新了
    for (auto &p : points)
    {
        point_update(p, TIME_INTERVAL_S);
        LOG_INFO << "第二次更新后点迹位置：" << p;
    }
    track_manager.draw_point_command(points); // 绘制点迹
    alg.process(points, new_tracks);

    // 第四次更新的时候可能会发送数据，作benchmark
    for (auto &p : points)
    {
        point_update(p, TIME_INTERVAL_S);
        LOG_INFO << "第三次更新后点迹位置：" << p;
    }
    track_manager.draw_point_command(points); // 绘制点迹

    alg.process(points, new_tracks);

    // std::this_thread::sleep_for(std::chrono::seconds(5));
    while (g_running) {
        wait_seconds(5);  // 可被 CTRL+C 中断的等待
    }
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