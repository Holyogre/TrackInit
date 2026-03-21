#pragma once
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
    lon = x / KM_PER_DEGREE + track_project::BASE_LONGITUDE;

    double cos_base = std::cos(track_project::BASE_LATITUDE * M_PI / 180.0);
    if (std::abs(cos_base) < 1e-9)
    {
        lat = track_project::BASE_LATITUDE;
    }
    else
    {
        lat = y / (KM_PER_DEGREE * cos_base) + track_project::BASE_LATITUDE;
    }
}

void ll2xy(double lon, double lat, double &x, double &y)
{
    // 简单的线性转换，假设1度经度约等于111km，1度纬度约等于111km
    x = (lon - track_project::BASE_LONGITUDE) * KM_PER_DEGREE;
    y = (lat - track_project::BASE_LATITUDE) * KM_PER_DEGREE * std::cos(track_project::BASE_LATITUDE * M_PI / 180.0);
}

void sync_lon_lat(TrackPoint &p)
{
    xy2ll(p.x, p.y, p.longitude, p.latitude);
}

// ——————————————————————————————————外推函数————————————————————————————————//
void point_update_cv(TrackPoint &p, double time_interval_s)
{
    // 更新位置
    p.x += p.vx * time_interval_s / 1000.0; // vx单位m/s，转换为km
    p.y += p.vy * time_interval_s / 1000.0;

    // 更新DOPPLER
    double range = std::sqrt(p.x * p.x + p.y * p.y);
    if (range > 1e-6)
    {
        double los_x = -p.x / range;
        double los_y = -p.y / range;
        p.doppler = p.vx * los_x + p.vy * los_y;
    }

    // 更新经纬度
    sync_lon_lat(p);

    // 更新时间戳
    p.time.milliseconds += static_cast<int64_t>(time_interval_s * 1000);
}

// 带噪声的外推函数
void point_update_cv_with_noise(
    TrackPoint &p,
    double time_interval_s,
    unsigned seed,
    double pos_noise_sigma_km = 0.1,  // 位置噪声标准差 (km)，默认10米
    double doppler_noise_sigma = 0.1) // 多普勒噪声标准差 (m/s)
{
    // 用seed创建随机数生成器
    std::mt19937 gen(seed);
    std::normal_distribution<> gauss{0.0, 1.0};

    // 保存原始速度
    double vx_orig = p.vx;
    double vy_orig = p.vy;

    // 更新位置（先不加噪声，保持运动模型纯净）
    p.x += p.vx * time_interval_s / 1000.0;
    p.y += p.vy * time_interval_s / 1000.0;

    // 加位置测量噪声
    p.x += gauss(gen) * pos_noise_sigma_km;
    p.y += gauss(gen) * pos_noise_sigma_km;

    // 更新DOPPLER（基于真实速度）
    double range = std::sqrt(p.x * p.x + p.y * p.y);
    if (range > 1e-6)
    {
        double los_x = -p.x / range;
        double los_y = -p.y / range;
        p.doppler = p.vx * los_x + p.vy * los_y;

        // 加多普勒测量噪声
        p.doppler += gauss(gen) * doppler_noise_sigma;
    }

    // 恢复速度（如果不希望速度被改变）
    p.vx = vx_orig;
    p.vy = vy_orig;

    // 更新经纬度
    sync_lon_lat(p);

    // 更新时间戳
    p.time.milliseconds += static_cast<int64_t>(time_interval_s * 1000);
}

// ——————————————————————————————————点迹生成函数————————————————————————————————//
/*****************************************************************************
 * @brief 使用直角坐标参数生成点迹：x, y, vx, vy
 *
 * @param time 时间戳，单位毫秒
 * @param params 参数列表，每个元素为[x, y, vx, vy]
 * @return std::vector<TrackPoint> 生成的点迹群
 *****************************************************************************/
std::vector<TrackPoint> generate_target_points_xyv(int64_t time, const std::vector<std::array<double, 4>> &params)
{
    std::vector<TrackPoint> pts;
    pts.reserve(params.size());

    for (const auto &param : params)
    {
        TrackPoint p{};

        // 直接使用传入的参数
        p.x = param[0];
        p.y = param[1];
        p.vx = param[2];
        p.vy = param[3];

        // 计算SOG和COG
        p.sog = std::sqrt(p.vx * p.vx + p.vy * p.vy);
        p.cog = std::atan2(p.vx, p.vy) * 180.0 / M_PI; // 北偏东，atan2(y,x)给出与x轴夹角
        if (p.cog < 0)
            p.cog += 360.0;

        // 多普勒速度 (假设雷达在原点)
        p.doppler = -(p.vx * p.x + p.vy * p.y) / hypot(p.x, p.y); // 投影到航向上

        // 时间戳
        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);

        pts.push_back(p);
    }
    return pts;
}

/*****************************************************************************
 * @brief 使用极坐标参数生成点迹：theta(deg), rho(km), sog(knot), cog(deg)
 *        theta和cog均为北偏东
 *
 * @param time 时间戳，单位毫秒
 * @param params 参数列表，每个元素为[theta, rho, sog, cog]
 * @return std::vector<TrackPoint> 生成的点迹群
 *****************************************************************************/
std::vector<TrackPoint> generate_target_points_polar(int64_t time,
                                                     const std::vector<std::array<double, 4>> &params)
{
    std::vector<TrackPoint> pts;
    pts.reserve(params.size());

    for (const auto &param : params)
    {
        TrackPoint p{};

        double theta_deg = param[0]; // 方位角，度，北偏东
        double rho = param[1];       // 距离，千米
        p.sog = param[2];            // 航速，节
        p.cog = param[3];            // 航向，度，北偏东

        // 极坐标转直角坐标 (北偏东)
        double theta_rad = theta_deg * M_PI / 180.0;
        p.x = rho * std::sin(theta_rad); // x = rho * sin(theta)
        p.y = rho * std::cos(theta_rad); // y = rho * cos(theta)

        // 航向转速度分量 (北偏东)
        double cog_rad = p.cog * M_PI / 180.0;
        p.vx = p.sog * std::sin(cog_rad); // vx = sog * sin(cog)
        p.vy = p.sog * std::cos(cog_rad); // vy = sog * cos(cog)

        // 多普勒速度 (假设雷达在原点)
        p.doppler = -(p.vx * p.x + p.vy * p.y) / hypot(p.x, p.y); // 投影到航向上

        // 时间戳
        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);

        pts.push_back(p);
    }

    return pts;
}

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
        p.doppler = -(p.vx * p.x + p.vy * p.y) / hypot(p.x, p.y);

        // 时间戳
        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);

        pts.push_back(p);
    }
    return pts;
}

/*****************************************************************************
 * @brief 用于生成符合高斯分布的均匀分布的点迹群，默认法线方向指向正北方向,
 * 雷达在x=0,y=0处
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

        // 速度分量,cog为北偏东
        p.vx = p.sog * std::sin(cog_rad); // vx = sog * sin(cog)
        p.vy = p.sog * std::cos(cog_rad); // vy = sog * cos(cog)

        // 2. 视线方向的单位向量（从雷达到目标）
        double range = std::sqrt(p.x * p.x + p.y * p.y);
        if (range < 1e-6)
        {
            LOG_INFO << "生成点迹时，点迹位置过于接近雷达站，已调整位置以避免除零错误";
            p.x = 1e-5;
            p.y = 0;
            range = 1e-5;
        } // 避免除零
        double los_x = -p.x / range; // 指向航向的向量的x分量（东）
        double los_y = -p.y / range; // 指向航向的向量的y分量（北）
        p.doppler = p.vx * los_x + p.vy * los_y;

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
        p.doppler = -(p.vx * p.x + p.vy * p.y) / hypot(p.x, p.y);

        p.time.milliseconds = time;

        // 经纬度
        sync_lon_lat(p);
        pts.push_back(p);
    }
    return pts;
}
