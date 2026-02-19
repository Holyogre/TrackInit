#include <cstdint>
#ifndef _TRACK_PROJECT_DEF_INIT_H
#define _TRACK_PROJECT_DEF_INIT_H
namespace track_project::trackinit
{
    // ============ SLICEHOUGH 聚类参数 ============
    constexpr double SLICEHOUGH_CLUSTER_RADIUS_KM = 20.0;    // 聚类切片直径（KM）
    constexpr double SLICEHOUGH_CORE_POINT_RADIUS_KM = 5.0; // 核心点点迹邻域半径（KM）
    constexpr int SLICEHOUGH_MAX_POINTS_PER_CLUSTER = 4;     // 每个聚类中最多可以拥有的点迹数量
    constexpr int SLICEHOUGH_MIN_POINTS_PER_CORE = 3;        // 每个核心点最少需要拥有的点迹数量
    // ============ SLICEHOUGH 霍夫变换参数 ============
    constexpr double SLICEHOUGH_THETA_RESOLUTION_DEG = 1.0; // 霍夫角度分辨率（度）,得是180的公因数
    constexpr double SLICEHOUGH_RHO_RESOLUTION_KM = 0.01;    // 霍夫距离分辨率（公里）

    /*****************************************************************************
     * @brief 状态码
     * 1. 0：处理成功，有航迹输出
     * 2. 1000-1100：霍夫变换处理成功，但推荐更换处理方式
     * 3.
     *****************************************************************************/
    enum class ProcessStatus : std::int32_t
    {
        SUCCESS = 0, // 处理成功，有航迹输出

        // 霍夫变换类状态码
        POINTS_TOO_DISPERSED = 1001, // 点迹过于分散,部分点迹被忽略
        TOO_MANY_CLUSTER = 1002,     // 聚类数量过多，可以更换成其他算法
        TOO_LARGE_CLUSTER = 1003,    // 过大的聚类区域，群目标存在丢失情况，考虑更换参数

        // 其他
        //  VELOCITY_OUT_OF_RANGE, // 速度约束不满足
        //  ACCELERATION_TOO_HIGH, // 机动性太强
        //  LOW_CONFIDENCE,        // 置信度过低
        //  AMBIGUOUS_TRACKS,      // 存在多个相似航迹假设

    };

}

#endif // _TRACK_PROJECT_DEF_INIT_H