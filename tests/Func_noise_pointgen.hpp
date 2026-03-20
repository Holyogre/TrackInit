#pragma once
#include "../include/defsystem.h"
#include "../include/defstruct.h"
#include <vector>
#include <random>
#include <cmath>
#include <memory>
#include <algorithm>

using track_project::TrackPoint;

/**
 * @brief 航迹外推类
 *
 * 对外只暴露update接口，内部维护误差分布表格和随机数生成器。
 * 支持设置随机数种子以保证测试结果的可重复性。
 */
class TrackExtrapolator
{
public:
    /**
     * @brief 构造函数
     * @param seed 随机数种子，用于控制噪声生成的随机性，保证测试结果稳定
     */
    explicit TrackExtrapolator(unsigned int seed = std::random_device{}())
        : rng_(seed)
    {
        // 初始化误差分布表格
        initErrorDistributionTable();
    }

    void update(std::vector<TrackPoint> &points, double time_interval_s)
    {
        for (auto &point : points)
        {
            update_single_point(point, time_interval_s);
        }
    }

    /*****************************************************************************
     * @brief 获取指定位置的误差分布（x和y方向的标准差）
     *
     * @param x x坐标（千米）
     * @param y y坐标（千米）
     * @return std::pair<double, double> (sigma_x, sigma_y) 单位：千米
     *****************************************************************************/
    std::pair<double, double> getErrorDistribution(double x, double y) const
    {
        size_t bin = location_to_bin_index(x, y);

        if (bin < error_distribution_table_.size())
        {
            return error_distribution_table_[bin];
        }
        return {0.0, 0.0}; // 默认值
    }

private:
    /*****************************************************************************
     * @brief 单点更新（带噪声）
     *
     * @param p 待更新的航迹点
     * @param time_interval_s 时间间隔（秒）
     *****************************************************************************/
    void update_single_point(TrackPoint &p,
                             double time_interval_s)
    {
        // 保存原始速度
        double vx_orig = p.vx;
        double vy_orig = p.vy;

        // 运动学更新（无噪声）
        p.x += p.vx * time_interval_s / 1000.0;
        p.y += p.vy * time_interval_s / 1000.0;

        // 提取目标位置的误差分布
        auto [sigma_x, sigma_y] = getErrorDistribution(p.x, p.y);
        if (sigma_x > 0 || sigma_y > 0)
        {
            std::normal_distribution<> gauss_x{0.0, sigma_x};
            std::normal_distribution<> gauss_y{0.0, sigma_y};
            p.x += gauss_x(rng_);
            p.y += gauss_y(rng_);
        }

        // 更新多普勒（基于真实速度）
        double range = std::sqrt(p.x * p.x + p.y * p.y);
        if (range > 1e-6)
        {
            double los_x = -p.x / range;
            double los_y = -p.y / range;
            p.doppler = p.vx * los_x + p.vy * los_y;

            // 添加多普勒噪声
            std::normal_distribution<> doppler_noise{0.0, track_project::BASE_DOPPLER_TOLERANCE_M_S};
            p.doppler += doppler_noise(rng_);
        }

        // 恢复速度（保持不变）
        p.vx = vx_orig;
        p.vy = vy_orig;

        // 更新经纬度
        sync_lon_lat(p);

        // 更新时间戳
        p.time.milliseconds += static_cast<int64_t>(time_interval_s * 1000);
    }

    // 坐标转换函数（保持不变）
    void xy2ll(double x, double y, double &lon, double &lat)
    {
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
        x = (lon - track_project::BASE_LONGITUDE) * KM_PER_DEGREE;
        y = (lat - track_project::BASE_LATITUDE) * KM_PER_DEGREE * std::cos(track_project::BASE_LATITUDE * M_PI / 180.0);
    }

    void sync_lon_lat(TrackPoint &p)
    {
        xy2ll(p.x, p.y, p.longitude, p.latitude);
    }

    // 初始化误差分布表格
    void initErrorDistributionTable()
    {
        error_distribution_table_.resize(MAX_BINS);

        for (size_t i = 0; i < MAX_BINS; ++i)
        {
            // 计算bin对应的x,y坐标
            auto [x_index, y_index] = bin_index_to_xy_index(i);
            double sbcpp_x = static_cast<double>(x_index);
            double sbcpp_y = static_cast<double>(y_index);

            double x = (sbcpp_x + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_X / LOGIC_BASED_NUM_X_BINS) - LOGIC_BASED_MAX_ABS_X;
            double y = (sbcpp_y + 0.5) * static_cast<double>(2 * LOGIC_BASED_MAX_ABS_Y / LOGIC_BASED_NUM_Y_BINS) - LOGIC_BASED_MAX_ABS_Y;

            // 转换为极坐标
            double rho = std::sqrt(x * x + y * y);
            double theta = std::atan2(y, x); // 弧度，范围 -π 到 π

            // 常数
            double sigma_r = track_project::BASE_R_SIGMA_KM;
            double sigma_theta_rad = track_project::BASE_THETA_SIGMA_DEG * M_PI / 180.0; // 度转弧度

            // 计算x方向上的方差：sigma_x^2 = sigma_r^2 * cos^2(theta) + rho^2 * sigma_theta^2 * sin^2(theta)
            double sigma_x2 = std::pow(sigma_r * std::cos(theta), 2) +
                              std::pow(rho * sigma_theta_rad * std::sin(theta), 2);

            // 计算y方向上的方差：sigma_y^2 = sigma_r^2 * sin^2(theta) + rho^2 * sigma_theta^2 * cos^2(theta)
            double sigma_y2 = std::pow(sigma_r * std::sin(theta), 2) +
                              std::pow(rho * sigma_theta_rad * std::cos(theta), 2);

            double sigma_x = std::sqrt(sigma_x2);
            double sigma_y = std::sqrt(sigma_y2);

            error_distribution_table_[i] = std::make_pair(sigma_x, sigma_y);
        }
    }

    // 计算x,y坐标对应的离散x,y索引
    std::pair<size_t, size_t> location_to_xy_index(double x, double y) const
    {
        //  计算分辨率
        double x_bin_size = (2 * LOGIC_BASED_MAX_ABS_X) / LOGIC_BASED_NUM_X_BINS; // X轴每个bin的大小
        double y_bin_size = (2 * LOGIC_BASED_MAX_ABS_Y) / LOGIC_BASED_NUM_Y_BINS; // Y轴每个bin的大小

        // 离散化输入坐标，四舍五入，计算对应的x,y索引
        size_t x_index = static_cast<size_t>(std::round((x + LOGIC_BASED_MAX_ABS_X) / x_bin_size));
        size_t y_index = static_cast<size_t>(std::round((y + LOGIC_BASED_MAX_ABS_Y) / y_bin_size));
        x_index = std::clamp(x_index, static_cast<size_t>(0), LOGIC_BASED_NUM_X_BINS - 1);
        y_index = std::clamp(y_index, static_cast<size_t>(0), LOGIC_BASED_NUM_Y_BINS - 1);

        return std::make_pair(x_index, y_index);
    }

    // 计算BIN值对应的索引
    size_t location_to_bin_index(double x, double y) const
    {
        // 调用location_to_xy_index获取x_index和y_index
        auto [x_index, y_index] = location_to_xy_index(x, y);
        size_t bin_index = xy_index_to_bin_index(x_index, y_index);

        return bin_index;
    }

    size_t xy_index_to_bin_index(size_t x_index, size_t y_index) const
    {
        return x_index + y_index * LOGIC_BASED_NUM_X_BINS;
    }

    std::pair<size_t, size_t> bin_index_to_xy_index(size_t bin_index) const
    {
        size_t x_index = bin_index % LOGIC_BASED_NUM_X_BINS;
        size_t y_index = bin_index / LOGIC_BASED_NUM_X_BINS;
        return std::make_pair(x_index, y_index);
    }

private:
    // 随机数生成器
    std::mt19937 rng_;

    // 误差分布表格：每个bin存储(sigma_x, sigma_y)
    std::vector<std::pair<double, double>> error_distribution_table_;

    // 表格维度常量（需要根据实际定义调整）
    static constexpr size_t LOGIC_BASED_NUM_X_BINS = 800;  // X轴bin数量
    static constexpr size_t LOGIC_BASED_NUM_Y_BINS = 800;  // Y轴bin数量
    static constexpr double LOGIC_BASED_MAX_ABS_X = 400.0; // X轴最大绝对值（千米）
    static constexpr double LOGIC_BASED_MAX_ABS_Y = 400.0; // Y轴最大绝对值（千米）
    static constexpr size_t MAX_BINS = LOGIC_BASED_NUM_X_BINS * LOGIC_BASED_NUM_Y_BINS;

    static constexpr double KM_PER_DEGREE = 111.0;
};