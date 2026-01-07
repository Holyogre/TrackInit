#include <catch2/catch_all.hpp>
#include "../src/HoughSlice.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <memory>
#include <algorithm>
#include <cstring>

#include "Func_pointgen.hpp" //测试点迹生成以及更新函数
#include "Func_cout.hpp"     //一些<<重载，方便输出
#include "../utils/Logger.hpp"

// 检测时间间隔，单位秒
constexpr double TIME_INTERVAL_S = 10.0;

using track_project::TrackPoint;
using track_project::trackinit::SliceHough;
using track_project::trackinit::SLICEHOUGH_CLUSTER_RADIUS_KM;

// 友元访问器，实现对私有聚类函数的调用
namespace track_project::trackinit
{
    struct SliceHoughBenchAccessor
    {
        static void run_gen_1(SliceHough &alg, const std::vector<TrackPoint> &pts) { alg.clust_gen(pts); } //原为clust_gen_1，经过测试，保留
        static void run_gen_2(SliceHough &alg, const std::vector<TrackPoint> &pts) { alg.clust_gen(pts); } //原为clust_gen_2，经过测试，删除gen_2
    };
}

using track_project::trackinit::SliceHoughBenchAccessor;

TEST_CASE("性能比较", "[clust_gen][bench_mark]")
{
    // 随机数种子
    unsigned int seed = Catch::getSeed();
    std::mt19937 gen(seed);

    // 参数定义
    const size_t N_POINTS = 1000;
    const size_t N_CLUSTERS_5 = 5;
    const size_t POINTS_PER_CLUSTER = N_POINTS / N_CLUSTERS_5;
    const double x_min = -400.0, x_max = 400.0;
    const double y_min = -400.0, y_max = 400.0;
    const double v_mean = 10.0, v_stddev = 5.0;
    const int x_cell = 1 + static_cast<int>((x_max - x_min) / SLICEHOUGH_CLUSTER_RADIUS_KM);
    const int y_cell = 1 + static_cast<int>((y_max - y_min) / SLICEHOUGH_CLUSTER_RADIUS_KM);

    // 生成数据
    std::vector<TrackPoint> uniform_point, gauss_point, rayleigh_point;
    std::vector<std::pair<double, double>> centers;
    std::set<std::pair<int, int>> used_cells;

    // 生成噪声点迹
    auto noise_point = generate_uniform_points(N_POINTS, 0, x_min, x_max, y_min, y_max, v_mean, v_stddev, seed);

    // 生成聚类中心 - 使用普通随机数生成器
    std::uniform_int_distribution<> dist_x(0, x_cell - 1); // 注意边界
    std::uniform_int_distribution<> dist_y(0, y_cell - 1);

    // 抽调一个聚类中心用于构造新航迹
    used_cells.insert({0, 0});
    auto single_clust_point = generate_gaussian_points(
        POINTS_PER_CLUSTER, TIME_INTERVAL_S,
        0, SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
        0, SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
        v_mean, v_stddev, seed);

    while (centers.size() < N_CLUSTERS_5)
    {
        int cell_x = dist_x(gen);
        int cell_y = dist_y(gen);

        if (used_cells.find({cell_x, cell_y}) == used_cells.end())
        {
            used_cells.insert({cell_x, cell_y});
            double center_x = x_min + (cell_x + 0.5) * SLICEHOUGH_CLUSTER_RADIUS_KM;
            double center_y = y_min + (cell_y + 0.5) * SLICEHOUGH_CLUSTER_RADIUS_KM;
            centers.emplace_back(center_x, center_y);
        }
    }

    LOG_INFO << "生成聚类中心点完成，聚类中心位置" << centers;

    // 为每个聚类中心生成三种分布的点
    for (const auto &center : centers)
    {
        // 均匀分布点
        auto temp_uniform = generate_uniform_points(
            POINTS_PER_CLUSTER, 0,
            center.first - SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
            center.first + SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
            center.second - SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
            center.second + SLICEHOUGH_CLUSTER_RADIUS_KM / 2,
            v_mean, v_stddev, seed);
        uniform_point.insert(uniform_point.end(), temp_uniform.begin(), temp_uniform.end());

        // 高斯分布点
        auto temp_gauss = generate_gaussian_points(
            POINTS_PER_CLUSTER, 0,
            center.first, SLICEHOUGH_CLUSTER_RADIUS_KM / 4,
            center.second, SLICEHOUGH_CLUSTER_RADIUS_KM / 4,
            v_mean, v_stddev, seed);
        gauss_point.insert(gauss_point.end(), temp_gauss.begin(), temp_gauss.end());

        // 瑞利分布点
        auto temp_rayleigh = generate_rayleigh_points(
            POINTS_PER_CLUSTER, 0,
            center.first, center.second,
            SLICEHOUGH_CLUSTER_RADIUS_KM / 4, 1.0,
            v_mean, v_stddev, seed);
        rayleigh_point.insert(rayleigh_point.end(), temp_rayleigh.begin(), temp_rayleigh.end());
    }

    // 存入Accessor构造第一批聚类
    SECTION("完全随机情况下的点迹更新测试")
    {
        SliceHough alg;
        SliceHoughBenchAccessor::run_gen_1(alg, noise_point); // 预先放入一批数据

        // 更新所有点迹
        for (auto &p : noise_point)
        {
            point_update(p, TIME_INTERVAL_S);
        }

        BENCHMARK("gen1函数")
        {
            SliceHoughBenchAccessor::run_gen_1(alg, noise_point);
        };

        BENCHMARK("gen2函数")
        {
            SliceHoughBenchAccessor::run_gen_2(alg, noise_point);
        };
    }

    SECTION("五个均匀分布聚类测试")
    {
        SliceHough alg;
        SliceHoughBenchAccessor::run_gen_1(alg, uniform_point); // 预先放入一批数据

        // 更新所有点迹
        for (auto &p : uniform_point)
        {
            point_update(p, TIME_INTERVAL_S);
        }

        // 添加噪声点和聚类点
        uniform_point.insert(uniform_point.end(), noise_point.begin(), noise_point.end());
        uniform_point.insert(uniform_point.end(), single_clust_point.begin(), single_clust_point.end());

        BENCHMARK("gen1函数")
        {
            SliceHoughBenchAccessor::run_gen_1(alg, uniform_point);
        };

        BENCHMARK("gen2函数")
        {
            SliceHoughBenchAccessor::run_gen_2(alg, uniform_point);
        };
    }

    SECTION("五个高斯分布聚类测试")
    {
        SliceHough alg;
        SliceHoughBenchAccessor::run_gen_1(alg, gauss_point); // 预先放入一批数据

        // 更新所有点迹
        for (auto &p : gauss_point)
        {
            point_update(p, TIME_INTERVAL_S);
        }

        // 添加噪声点和聚类点
        gauss_point.insert(gauss_point.end(), noise_point.begin(), noise_point.end());
        gauss_point.insert(gauss_point.end(), single_clust_point.begin(), single_clust_point.end());

        BENCHMARK("gen1函数")
        {
            SliceHoughBenchAccessor::run_gen_1(alg, gauss_point);
        };

        BENCHMARK("gen2函数")
        {
            SliceHoughBenchAccessor::run_gen_2(alg, gauss_point);
        };
    }

    SECTION("五个瑞利分布聚类测试")
    {
        SliceHough alg;
        SliceHoughBenchAccessor::run_gen_1(alg, rayleigh_point); // 预先放入一批数据

        // 更新所有点迹
        for (auto &p : rayleigh_point)
        {
            point_update(p, TIME_INTERVAL_S);
        }

        // 添加噪声点和聚类点
        rayleigh_point.insert(rayleigh_point.end(), noise_point.begin(), noise_point.end());
        rayleigh_point.insert(rayleigh_point.end(), single_clust_point.begin(), single_clust_point.end());

        BENCHMARK("gen1函数")
        {
            SliceHoughBenchAccessor::run_gen_1(alg, rayleigh_point);
        };

        BENCHMARK("gen2函数")
        {
            SliceHoughBenchAccessor::run_gen_2(alg, rayleigh_point);
        };
    }
}