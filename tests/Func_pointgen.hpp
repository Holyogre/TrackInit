#include "../src/HoughSlice.hpp"
#include "../include/defsystem.h"
#include <vector>
#include <random>
#include <cmath>

using track_project::TrackPoint;

constexpr double KM_PER_DEGREE = 111.0;

// ——————————————————————————————————坐标转换函数————————————————————————————————//
void xy2ll(double x, double y, double &lon, double &lat)
{
    // 反向转换：假设1度经纬对应固定km比例
    lon = x / KM_PER_DEGREE + track_project::base_longitude;

    double cos_base = std::cos(track_project::base_latitude * M_PI / 180.0);
    if (std::abs(cos_base) < 1e-9)
    {
        lat = track_project::base_latitude;
    }
    else
    {
        lat = y / (KM_PER_DEGREE * cos_base) + track_project::base_latitude;
    }
}

void ll2xy(double lon, double lat, double &x, double &y)
{
    // 简单的线性转换，假设1度经度约等于111km，1度纬度约等于111km
    x = (lon - track_project::base_longitude) * KM_PER_DEGREE;
    y = (lat - track_project::base_latitude) * KM_PER_DEGREE * std::cos(track_project::base_latitude * M_PI / 180.0);
}

void sync_lon_lat(TrackPoint &p)
{
    xy2ll(p.x, p.y, p.longitude, p.latitude);
}

// ——————————————————————————————————外推函数————————————————————————————————//
void point_update(TrackPoint &p, double time_interval_s)
{
    // 更新位置
    p.x += p.vx * time_interval_s / 1000.0; // vx单位m/s，转换为km
    p.y += p.vy * time_interval_s / 1000.0;

    // 更新经纬度
    sync_lon_lat(p);

    // 更新时间戳
    p.time.milliseconds += static_cast<int64_t>(time_interval_s * 1000);
}

// ——————————————————————————————————点迹生成函数————————————————————————————————//
/*****************************************************************************
 * @brief 用于生成符合点迹生成参数的均匀分布的点迹群
 *
 * @param n_points 要求生成的点迹数量
 * @param time 时间戳，单位毫秒
 * @param x_min,x_max 点迹x坐标范围
 * @param y_min,y_max 点迹y坐标范围
 * @param v_mean,v_stddev 点迹速度均值方差，据参考文献来看，这两个均符合高斯分布，没必要强行均匀分布
 * @param seed 随机数种子
 * @return std::vector<TrackPoint> 生成的点迹群
 *****************************************************************************/
std::vector<TrackPoint> generate_uniform_points(size_t n_points, int64_t time,
                                                double x_min, double x_max,
                                                double y_min, double y_max,
                                                double v_mean, double v_stddev,
                                                unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist_x(x_min, x_max);
    std::uniform_real_distribution<double> dist_y(y_min, y_max);
    std::normal_distribution<double> dist_sog(v_mean, v_stddev);
    std::uniform_real_distribution<double> dist_cog(0.0, 360.0); // 航向均匀分布

    std::vector<TrackPoint> pts;
    pts.reserve(n_points);

    for (size_t i = 0; i < n_points; ++i)
    {
        TrackPoint p{};

        // 位置 (均匀分布)
        p.x = dist_x(rng);
        p.y = dist_y(rng);

        // 航速SOG (高斯分布)
        p.sog = std::abs(dist_sog(rng));

        // 航向COG (均匀分布)
        p.cog = dist_cog(rng);
        double cog_rad = p.cog * M_PI / 180.0;

        // 速度分量
        p.vx = p.sog * std::cos(cog_rad);
        p.vy = p.sog * std::sin(cog_rad);

        // 多普勒速度 (假设雷达在原点)
        double point_angle = std::atan2(p.y, p.x);
        p.doppler = p.sog * std::cos(point_angle - cog_rad);

        // 时间戳
        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);

        pts.push_back(p);
    }
    return pts;
}

/*****************************************************************************
 * @brief 用于生成符合高斯分布的均匀分布的点迹群
 *
 * @param n_points 要求生成的点迹数量
 * @param time 时间戳，单位毫秒
 * @param x_mean,x_stddev 点迹x坐标均值和标准差
 * @param y_mean,y_stddev 点迹y坐标均值和标准差
 * @param v_mean,v_stddev 点迹速度均值和标准差
 * @param seed 随机数种子
 * @return std::vector<TrackPoint> 生成的点迹群
 *****************************************************************************/
std::vector<TrackPoint> generate_gaussian_points(size_t n_points, int64_t time,
                                                 double x_mean, double x_stddev,
                                                 double y_mean, double y_stddev,
                                                 double v_mean, double v_stddev,
                                                 unsigned seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist_x(x_mean, x_stddev);
    std::normal_distribution<double> dist_y(y_mean, y_stddev);
    std::normal_distribution<double> dist_sog(v_mean, v_stddev);
    std::uniform_real_distribution<double> dist_cog(0.0, 360.0);

    std::vector<TrackPoint> pts;
    pts.reserve(n_points);

    for (size_t i = 0; i < n_points; ++i)
    {
        TrackPoint p{};

        // 位置 (高斯分布)
        p.x = dist_x(rng);
        p.y = dist_y(rng);

        // 航速航向
        p.sog = std::abs(dist_sog(rng));
        p.cog = dist_cog(rng);
        double cog_rad = p.cog * M_PI / 180.0;

        // 速度分量
        p.vx = p.sog * std::cos(cog_rad);
        p.vy = p.sog * std::sin(cog_rad);

        // 多普勒速度
        double point_angle = std::atan2(p.y, p.x);
        p.doppler = p.sog * std::cos(point_angle - cog_rad);

        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);
        pts.push_back(p);
    }
    return pts;
}

/*****************************************************************************
 * @brief 用于生成符合瑞利分布的均匀分布的点迹群
 *
 * @param n_points 要求生成的点迹数量
 * @param time 时间戳，单位毫秒
 * @param center_x,center_y 瑞利分布中心点坐标
 * @param radius_mean,radius_std 半径的均值和标准差
 * @param v_mean,v_stddev 点迹速度均值和标准差
 * @param seed 随机数种子
 * @return std::vector<TrackPoint> 生成的点迹群
 *****************************************************************************/
std::vector<TrackPoint> generate_rayleigh_points(size_t n_points, int64_t time,
                                                 double center_x, double center_y,
                                                 double radius_mean, double radius_std,
                                                 double v_mean, double v_stddev,
                                                 unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::normal_distribution<double> normal_dist(0.0, radius_std);
    std::normal_distribution<double> dist_sog(v_mean, v_stddev);
    std::uniform_real_distribution<double> dist_cog(0.0, 360.0);

    std::vector<TrackPoint> pts;
    pts.reserve(n_points);

    for (size_t i = 0; i < n_points; ++i)
    {
        TrackPoint p{};

        // 位置 (瑞利分布 - 环状)
        double angle = angle_dist(rng);
        double u1 = normal_dist(rng);
        double u2 = normal_dist(rng);
        double radius = radius_mean + std::sqrt(u1 * u1 + u2 * u2);
        p.x = center_x + radius * std::cos(angle);
        p.y = center_y + radius * std::sin(angle);

        // 航速航向
        p.sog = std::abs(dist_sog(rng));
        p.cog = dist_cog(rng);
        double cog_rad = p.cog * M_PI / 180.0;

        // 速度分量
        p.vx = p.sog * std::cos(cog_rad);
        p.vy = p.sog * std::sin(cog_rad);

        // 多普勒速度
        double point_angle = std::atan2(p.y, p.x);
        p.doppler = p.sog * std::cos(point_angle - cog_rad);

        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);
        pts.push_back(p);
    }
    return pts;
}
