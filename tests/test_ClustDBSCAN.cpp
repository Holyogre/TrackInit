#include <catch2/catch_all.hpp>
#include "../src/Func_dbscan.hpp"
#include "../include/defstruct.h"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace track_project;
using namespace track_project::trackinit;

// 辅助函数：创建简单的 TrackPoint
track_project::TrackPoint create_point(double x, double y)
{
    track_project::TrackPoint point;
    point.x = x;
    point.y = y;
    point.longitude = 0.0;
    point.latitude = 0.0;
    point.sog = 0.0;
    point.cog = 0.0;
    point.vx = 0.0;
    point.vy = 0.0;
    point.is_associated = false;
    return point;
}

// 辅助函数：检查聚类结果是否与预期匹配
bool compare_clusters(const std::vector<std::vector<size_t>>& result,
                     const std::vector<std::vector<size_t>>& expected)
{
    if (result.size() != expected.size())
        return false;
    
    // 对每个聚类进行排序（因为 dbscan 内部已经排序了）
    auto sorted_result = result;
    auto sorted_expected = expected;
    
    // 对聚类内的索引排序
    for (auto& cluster : sorted_result)
        std::sort(cluster.begin(), cluster.end());
    for (auto& cluster : sorted_expected)
        std::sort(cluster.begin(), cluster.end());
    
    // 对聚类本身排序（按第一个元素）
    std::sort(sorted_result.begin(), sorted_result.end());
    std::sort(sorted_expected.begin(), sorted_expected.end());
    
    return sorted_result == sorted_expected;
}

TEST_CASE("DBSCAN - 基础功能测试", "[dbscan][basic]")
{
    SECTION("空数据集")
    {
        std::vector<TrackPoint> points;
        auto clusters = dbscan(points, 1.0, 3);
        
        REQUIRE(clusters.empty());
    }
    
    SECTION("单个点 - 无法形成聚类")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0)
        };
        
        auto clusters = dbscan(points, 1.0, 2);
        
        REQUIRE(clusters.empty());
    }
    
    SECTION("两个距离较远的点 - 无法形成聚类")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(10.0, 10.0)
        };
        
        auto clusters = dbscan(points, 1.0, 2);
        
        REQUIRE(clusters.empty());
    }
    
    SECTION("三个紧密的点 - 形成一个聚类")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(0.5, 0.0),
            create_point(0.0, 0.5)
        };
        
        auto clusters = dbscan(points, 1.0, 2);
        
        REQUIRE(clusters.size() == 1);
        REQUIRE(clusters[0].size() == 3);
        
        // 检查是否包含所有点
        std::vector<size_t> cluster = clusters[0];
        std::sort(cluster.begin(), cluster.end());
        REQUIRE(cluster == std::vector<size_t>{0, 1, 2});
    }
    
    SECTION("两个分离的聚类")
    {
        // 第一个聚类：点 (0,0), (0.5,0), (0,0.5)
        // 第二个聚类：点 (10,10), (10.5,10), (10,10.5)
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),    // 0
            create_point(0.5, 0.0),    // 1
            create_point(0.0, 0.5),    // 2
            create_point(10.0, 10.0),  // 3
            create_point(10.5, 10.0),  // 4
            create_point(10.0, 10.5)   // 5
        };
        
        auto clusters = dbscan(points, 1.0, 2);
        
        REQUIRE(clusters.size() == 2);
        
        // 对聚类排序以便比较
        std::sort(clusters.begin(), clusters.end(),
                 [](const auto& a, const auto& b) { return a[0] < b[0]; });
        
        REQUIRE(clusters[0].size() == 3);
        REQUIRE(clusters[1].size() == 3);
        
        // 检查聚类内容
        std::vector<size_t> cluster1 = clusters[0];
        std::vector<size_t> cluster2 = clusters[1];
        std::sort(cluster1.begin(), cluster1.end());
        std::sort(cluster2.begin(), cluster2.end());
        
        REQUIRE((cluster1 == std::vector<size_t>{0, 1, 2} || cluster1 == std::vector<size_t>{3, 4, 5}));
        REQUIRE((cluster2 == std::vector<size_t>{0, 1, 2} || cluster2 == std::vector<size_t>{3, 4, 5}));
        REQUIRE(cluster1 != cluster2);
    }
    
    SECTION("包含噪声点的聚类")
    {
        // 核心点：点 (0,0), (0.5,0), (0,0.5) 形成聚类
        // 噪声点：点 (10,10) 远离其他点
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),    // 0 - 聚类
            create_point(0.5, 0.0),    // 1 - 聚类
            create_point(0.0, 0.5),    // 2 - 聚类
            create_point(10.0, 10.0)   // 3 - 噪声
        };
        
        auto clusters = dbscan(points, 1.0, 2);
        
        // 应该只有一个聚类（三个点）
        REQUIRE(clusters.size() == 1);
        REQUIRE(clusters[0].size() == 3);
        
        std::vector<size_t> cluster = clusters[0];
        std::sort(cluster.begin(), cluster.end());
        REQUIRE(cluster == std::vector<size_t>{0, 1, 2});
    }
}

TEST_CASE("DBSCAN - 参数测试", "[dbscan][parameters]")
{
    SECTION("不同的 eps 值")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(0.8, 0.0),  // 距离 0.8
            create_point(0.0, 0.8)   // 距离 0.8
        };
        
        // eps=0.5，点之间距离 0.8 > 0.5，无法形成聚类
        auto clusters1 = dbscan(points, 0.5, 2);
        REQUIRE(clusters1.empty());
        
        // eps=1.0，点之间距离 0.8 < 1.0，可以形成聚类
        auto clusters2 = dbscan(points, 1.0, 2);
        REQUIRE(clusters2.size() == 1);
        REQUIRE(clusters2[0].size() == 3);
    }
    
    SECTION("不同的 min_pts 值")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(0.5, 0.0),
            create_point(0.0, 0.5)
        };
        
        // min_pts=3，刚好满足条件
        auto clusters1 = dbscan(points, 1.0, 3);
        REQUIRE(clusters1.size() == 1);
        REQUIRE(clusters1[0].size() == 3);
        
        // min_pts=4，无法满足条件
        auto clusters2 = dbscan(points, 1.0, 4);
        REQUIRE(clusters2.empty());
    }
    
    SECTION("min_pts=1 的特殊情况")
    {
        // 当 min_pts=1 时，每个点都是核心点
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(10.0, 10.0)
        };
        
        auto clusters = dbscan(points, 1.0, 1);
        
        // 每个点都应该形成自己的聚类
        REQUIRE(clusters.size() == 2);
        REQUIRE(clusters[0].size() == 1);
        REQUIRE(clusters[1].size() == 1);
    }
}

TEST_CASE("DBSCAN - 边界情况测试", "[dbscan][edge]")
{
    SECTION("所有点都相同")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(0.0, 0.0),
            create_point(0.0, 0.0),
            create_point(0.0, 0.0)
        };
        
        auto clusters = dbscan(points, 0.1, 2);
        
        // 所有相同的点应该形成一个聚类
        REQUIRE(clusters.size() == 1);
        REQUIRE(clusters[0].size() == 4);
    }
    
    SECTION("链式结构")
    {
        // 点形成一条链：每个点只与相邻的点距离足够近
        // 点0-点1距离0.6，点1-点2距离0.6，但点0-点2距离1.2
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),   // 0
            create_point(0.6, 0.0),   // 1
            create_point(1.2, 0.0)    // 2
        };
        
        // eps=0.7, min_pts=2
        // 点0和点1相互可达，点1和点2相互可达
        // 因此所有点应该属于同一个聚类
        auto clusters = dbscan(points, 0.7, 2);
        
        REQUIRE(clusters.size() == 1);
        REQUIRE(clusters[0].size() == 3);
    }
    
    SECTION("eps 非常小")
    {
        std::vector<TrackPoint> points = {
            create_point(0.0, 0.0),
            create_point(0.001, 0.001)
        };
        
        // eps 太小，无法检测到邻近点
        auto clusters = dbscan(points, 0.0001, 2);
        REQUIRE(clusters.empty());
        
        // eps 足够大，可以检测到
        auto clusters2 = dbscan(points, 0.01, 2);
        REQUIRE(clusters2.size() == 1);
        REQUIRE(clusters2[0].size() == 2);
    }
}

TEST_CASE("DBSCAN - 性能测试", "[dbscan][benchmark]")
{
    SECTION("大量点聚类")
    {
        const int NUM_POINTS = 100;
        std::vector<TrackPoint> points;
        points.reserve(NUM_POINTS);
        
        // 创建多个小聚类
        for (int i = 0; i < NUM_POINTS / 4; ++i)
        {
            // 第一个聚类
            points.push_back(create_point(i * 0.1, 0.0));
            // 第二个聚类
            points.push_back(create_point(i * 0.1, 10.0));
            // 第三个聚类
            points.push_back(create_point(i * 0.1, 20.0));
            // 第四个聚类
            points.push_back(create_point(i * 0.1, 30.0));
        }
        
        BENCHMARK("聚类 " + std::to_string(NUM_POINTS) + " 个点")
        {
            auto clusters = dbscan(points, 0.5, 2);
            return clusters.size();
        };
    }
}
