#include <cstdint>
#ifndef _TRACK_PROJECT_DEF_INIT_H
#define _TRACK_PROJECT_DEF_INIT_H
namespace track_project::trackinit
{
    // ============ HOUGHSLICE 聚类参数 ============
    constexpr double HOUGHSLICE_CLUSTER_RADIUS_KM = 30.0;    // 聚类切片直径（KM）
    constexpr double HOUGHSLICE_CORE_POINT_RADIUS_KM = 25.0; // 核心点点迹邻域半径（KM）
    constexpr int HOUGHSLICE_MAX_POINTS_PER_CLUSTER = 4;     // 每个聚类中最多可以拥有的点迹数量
    constexpr int HOUGHSLICE_MIN_POINTS_PER_CORE = 3;        // 每个核心点最少需要拥有的点迹数量
    // ============ HOUGHSLICE 霍夫变换参数 ============
    constexpr double HOUGHSLICE_THETA_RESOLUTION_DEG = 0.1;  // 霍夫角度分辨率（度）,得是180的公因数
    constexpr double HOUGHSLICE_RHO_RESOLUTION_KM = 0.1;     // 霍夫距离分辨率（公里）
    constexpr size_t HOUGHSLICE_DOPPLER_BIT_NUM = 64;        // 霍夫变换中多普勒速度位数，必须是32的整数倍
    constexpr size_t HOUGHSLICE_BATCH_NUM = 4;               // 批次数量，必须大于等于4,不然不能回溯四批次数据
    constexpr size_t HOUGHSLICE_THETA_CLUSTER_TOL_DEG = 1.0; // 霍夫空间中检测到的直线参数凝聚时的角度阈值（度），必须可以被角度分辨率整除
    constexpr double HOUGHSLICE_RHO_CLUSTER_TOL_KM = 0.1;    // 霍夫空间中检测到的直线参数凝聚时的距离阈值（公里），必须可以被距离分辨率整除

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
        NO_POINT = 2001, // 没有点迹，无法处理
        //  VELOCITY_OUT_OF_RANGE, // 速度约束不满足
        //  ACCELERATION_TOO_HIGH, // 机动性太强
        //  LOW_CONFIDENCE,        // 置信度过低
        //  AMBIGUOUS_TRACKS,      // 存在多个相似航迹假设

    };
}

#endif // _TRACK_PROJECT_DEF_INIT_H