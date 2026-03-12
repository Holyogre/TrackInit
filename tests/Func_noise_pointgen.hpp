#pragma once
#include <../tests/Func_pointgen.hpp>

// 超参数设置，避免混用
namespace pointgen
{
    // 允许出现的最远的点迹范围
    constexpr double MAX_ABS_X = 400.0; // 最大的X坐标的绝对值上限（KM），以雷达站为原点，正东向右
    constexpr double MIN_ABS_X = 10.0;  // 最小的X坐标的绝对值下限（KM），以雷达站为原点，正东向右
    constexpr double MAX_ABS_Y = 400.0; // 最大的Y坐标的绝对值上限（KM），以雷达站为原点，正北向上
    constexpr double MIN_ABS_Y = 10.0;  // 最小的Y坐标的绝对值下限（KM），以雷达站为原点，正北向上

    // 噪声参数设置
    constexpr double POS_NOISE_SIGMA_KM = 0.01;  // 位置噪声标准差 (km)，默认10米
    constexpr double DOPPLER_NOISE_SIGMA = 0.8; // 多普勒噪声标准差 (m/s)


};