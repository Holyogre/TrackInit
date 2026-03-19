#include "Func_dbscan.hpp"

#include <cstddef>
//傻逼代码报错多
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#include "../utils/nanoflann.hpp"
#pragma GCC diagnostic pop
#include "../utils/Logger.hpp"

#include <type_traits>
#include <vector>
#include <algorithm>

constexpr bool enable_single_point_clust = true; // 是否启用单点聚类

/*****************************************************************************
 * @brief nanoflann的适配器，需要提供（名字得一样）：
 * 1. points，原始数据的引用
 * 2. kdtree_get_point_count()，获取点的数量
 * 3. kdtree_get_pt(const std::size_t idx, const std::size_t dim),核心访问器，用于访问K纬的数据
 * 可选：
 * 1. kdtree_get_bbox()，边界框计算，详情参考nanoflann的实现
 * 自定义拓展
 * 1. elem_ptr()，获取指定的指针
 *****************************************************************************/
namespace track_project::trackinit
{
    // 给聚类用的适配器
    struct Cluster_adaptor
    {
        const std::vector<TrackPoint> &points; // 存储点的集合（原始数据的引用）

        // 通过引用绑定实现快速初始化
        Cluster_adaptor(const std::vector<TrackPoint> &points) : points(points) {}

        // 必须返回数据点的数量
        inline std::size_t kdtree_get_point_count() const { return points.size(); }

        // 返回类中第 idx 个点的 dim 维的值：
        inline double kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
        {
            if (dim == 0)
                return points[idx].x;
            return points[idx].y;
        }

        // 可选的边界框计算：如果返回 false 则会使用默认的 bbox 计算循环。
        // 如果返回 true，表示边界框已由类计算，并将其存储在 "bb" 中，这样就可以避免重复计算。
        // 可以查看 bb.size() 来查看期望的维度（例如，点云的 2 维或 3 维）
        template <class BBOX>
        bool kdtree_get_bbox(BBOX & /*bb*/) const { return false; }

        // 获取指定点的指针
        auto const *elem_ptr(const std::size_t idx) const
        {
            return &points[idx].x;
        }
    };

    // 对聚类结果进行排序
    auto sort_clusters(std::vector<std::vector<size_t>> &clusters)
    {
        for (auto &cluster : clusters)
        {
            std::sort(cluster.begin(), cluster.end());
        }
    }

    // DBSCAN 聚类算法
    template <int n_cols, typename Adaptor>
    auto dbscan(const Adaptor &adapt, double eps, int min_pts)
    {
        eps *= eps; // 将半径 eps 转为平方距离
        using namespace nanoflann;
        using my_kd_tree_t = KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, Adaptor>, // 距离度量适配器
            Adaptor,                            // 数据源适配器
            n_cols>;                            // 维度

        // 使用 k-d 树创建索引
        my_kd_tree_t index(n_cols, adapt, KDTreeSingleIndexAdaptorParams(10));
        index.buildIndex();

        const auto n_points = adapt.kdtree_get_point_count();        // 获取点的数量
        auto visited = std::vector<bool>(n_points);                  // 访问标志数组
        auto clusters = std::vector<std::vector<size_t>>();          // 聚类结果
        auto matches = std::vector<std::pair<size_t, double>>();     // 存储匹配的点对
        auto sub_matches = std::vector<std::pair<size_t, double>>(); // 存储子匹配的点对

        // 遍历所有点
        for (size_t i = 0; i < n_points; i++)
        {
            if (visited[i]) // 如果该点已经访问过，跳过
                continue;

            // 进行半径搜索，寻找与当前点距离小于 eps 的所有点
            index.radiusSearch(adapt.elem_ptr(i), eps, matches, nanoflann::SearchParams(32, 0.f, false));

            if (matches.size() < static_cast<size_t>(min_pts)) // 如果邻域内点数小于 min_pts，则跳过
                continue;

            visited[i] = true; // 标记当前点为已访问

            auto cluster = std::vector<size_t>({i}); // 创建一个新的聚类，包含当前点
            LOG_DEBUG << "形成新聚类，初始点索引：" << i << "，邻域点数量：" << matches.size();

            // 扩展聚类：对匹配的点进行遍历，直到没有新的点可以加入聚类
            while (!matches.empty())
            {
                auto nb_idx = matches.back().first; // 获取邻域点的索引
                matches.pop_back();
                if (visited[nb_idx]) // 如果邻域点已访问，跳过
                    continue;
                visited[nb_idx] = true; // 标记邻域点为已访问

                // 对邻域点进行半径搜索，寻找其邻域点
                index.radiusSearch(adapt.elem_ptr(nb_idx), eps, sub_matches, nanoflann::SearchParams(32, 0.f, false));

                // 如果邻域内点数大于 min_pts，则将这些点加入当前的 matches 中
                if (sub_matches.size() >= static_cast<size_t>(min_pts))
                {
                    std::copy(sub_matches.begin(), sub_matches.end(), std::back_inserter(matches));
                }
                cluster.push_back(nb_idx); // 将邻域点加入当前聚类
            }
            clusters.emplace_back(std::move(cluster)); // 将当前聚类添加到聚类列表中
        }
        sort_clusters(clusters); // 对聚类结果进行排序
        return clusters;
    }

    // 特殊的dbscan，用于单点聚类
    template <int n_cols, typename Adaptor>
    auto dbscan_single_point(const Adaptor &adapt, double eps, int min_pts)
        -> std::vector<std::vector<size_t>>
    {
        eps *= eps; // 将半径 eps 转为平方距离
        using namespace nanoflann;
        using my_kd_tree_t = KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, Adaptor>, // 距离度量适配器
            Adaptor,                            // 数据源适配器
            n_cols>;                            // 维度

        // 使用 k-d 树创建索引
        my_kd_tree_t index(n_cols, adapt, KDTreeSingleIndexAdaptorParams(10));
        index.buildIndex();

        const auto n_points = adapt.kdtree_get_point_count();        // 获取点的数量
        auto visited = std::vector<bool>(n_points);                  // 访问标志数组
        auto clusters = std::vector<std::vector<size_t>>();          // 聚类结果
        auto matches = std::vector<std::pair<size_t, double>>();     // 存储匹配的点对
        auto sub_matches = std::vector<std::pair<size_t, double>>(); // 存储子匹配的点对

        // 遍历所有点
        for (size_t i = 0; i < n_points; i++)
        {
            if (visited[i]) // 如果该点已经访问过，跳过
                continue;

            // 进行半径搜索，寻找与当前点距离小于 eps 的所有点
            index.radiusSearch(adapt.elem_ptr(i), eps, matches, nanoflann::SearchParams(32, 0.f, false));

            if (matches.size() >= static_cast<size_t>(min_pts))
            {
                // 如果邻域内点数大于等于 min_pts，则进行聚类扩展
                visited[i] = true; // 标记当前点为已访问

                auto cluster = std::vector<size_t>({i}); // 创建一个新的聚类，包含当前点
                LOG_DEBUG << "形成新聚类，初始点索引：" << i << "，邻域点数量：" << matches.size();

                // 扩展聚类：对匹配的点进行遍历，直到没有新的点可以加入聚类
                while (!matches.empty())
                {
                    auto nb_idx = matches.back().first; // 获取邻域点的索引
                    matches.pop_back();
                    if (visited[nb_idx]) // 如果邻域点已访问，跳过
                        continue;
                    visited[nb_idx] = true; // 标记邻域点为已访问

                    // 对邻域点进行半径搜索，寻找其邻域点
                    index.radiusSearch(adapt.elem_ptr(nb_idx), eps, sub_matches, nanoflann::SearchParams(32, 0.f, false));

                    // 如果邻域内点数大于 min_pts，则将这些点加入当前的 matches 中
                    if (sub_matches.size() >= static_cast<size_t>(min_pts))
                    {
                        std::copy(sub_matches.begin(), sub_matches.end(), std::back_inserter(matches));
                    }
                    cluster.push_back(nb_idx); // 将邻域点加入当前聚类
                }
                clusters.emplace_back(std::move(cluster)); // 将当前聚类添加到聚类列表中
            }
            else
            {
                // 如果邻域内点数小于 min_pts，则自己形成一个单点聚类
                visited[i] = true;                         // 标记当前点为已访问
                auto cluster = std::vector<size_t>({i});   // 创建单点聚类
                clusters.emplace_back(std::move(cluster)); // 将单点聚类添加到聚类列表中
                LOG_DEBUG << "形成单点聚类，点索引：" << i << "，邻域点数量：" << matches.size();
            }
        }
        sort_clusters(clusters); // 对聚类结果进行排序
        return clusters;
    }

    // 二维DBSCAN算法实现
    auto dbscan(const std::vector<TrackPoint> &data, double eps, int min_pts) -> std::vector<std::vector<size_t>>
    {
        const auto adapt = Cluster_adaptor(data); // 创建适配器

        if (!enable_single_point_clust)
        {
            return dbscan<2, Cluster_adaptor>(adapt, eps, min_pts); // 调用 DBSCAN 聚类算法，n_cols 为 2，表示二维点
        }
        else
        {
            return dbscan_single_point<2, Cluster_adaptor>(adapt, eps, min_pts); // 调用支持单点聚类的 DBSCAN 聚类算法
        }
    }
}