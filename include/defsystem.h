#pragma once

namespace track_project
{
    // 雷达站基础信息
    constexpr double BASE_LONGITUDE = 0.0;
    constexpr double BASE_LATITUDE = 0.0;
    constexpr double BASE_NORMAL_DIRECTION = 0.0;      // 法线方向：北偏东，单位°，范围仅允许填写[0,360)，不然结果错误
    constexpr double BASE_THETA_SIGMA_DEG = 0.4;       // 先验航向误差标准差，单位°，
    constexpr double BASE_R_SIGMA_KM = 1.0;            // 先验距离误差标准差，单位KM，
    constexpr double BASE_DOPPLER_TOLERANCE_M_S = 0.5; // DOPPLER容差，单位m/s，超过这个值就认为不满足条件了

    // 目标基础信息
    constexpr double VELOCITY_MAX = 250.0; // 理论目标最大速度，单位m/s

    // 跟踪基本信息
    constexpr size_t MAX_INPUT_POINTS = 500; // 可以接受的最大点迹输入数量，超过这个数量就直接放弃，避免算力浪费
    constexpr size_t MAX_TRACK_NUM = 1000; // 最大航迹数量，超过这个数量就直接放弃，避免算力浪费
}