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
    constexpr double HOUGHSLICE_RHO_CLUSTER_TOL_KM = 0.3;    // 霍夫空间中检测到的直线参数凝聚时的距离阈值（公里），必须可以被距离分辨率整除

    // ============ LOGIC BASED 空间参数 ============
    constexpr double LOGIC_BASED_MAX_ABS_X = 400.0; // 最大的X坐标的绝对值上限（KM），以雷达站为原点，正东向右
    constexpr double LOGIC_BASED_MAX_ABS_Y = 400.0; // 最大的Y坐标的绝对值上限（KM），以雷达站为原点，正北向上
    constexpr size_t LOGIC_BASED_NUM_X_BINS = 800;  // X轴方向离散单元数量，O(N)严重影响性能和内存占用
    constexpr size_t LOGIC_BASED_NUM_Y_BINS = 800;  // Y轴方向离散单元数量，O(N)严重影响性能和内存占用
    // ============ LOGIC BASED 外推分辨率参数 ============
    constexpr double LOGIC_BASED_HEADING_RESOLUTION_DEG = 2.0; // 航向离散分辨率（度），迅速提高虚假航迹分辨率，但是计算复杂度也会迅速提升，O(N)
    constexpr size_t LOGIC_BASED_MAX_NODE_PER_BINS = 100;      // 每个波门中允许拥有的最大假设数量,建议设置为n^2，减枝一次O(N)
    constexpr double LOGIC_BASED_PROTECTIVE_RADIUS_KM = 0.1;   // 保护半径（KM），误差分布函数不一定准确，以防万一设置的扩张值，谨慎修改，非DEBUG最好是0
    constexpr double LOGIC_BASED_PROTECTIVE_DOPPLER_M_S = 0.01; // doppler保护半径,这个保护半径用于抑制波动值，可以适当给大点
    // ============ LOGIC BASED 减枝参数 ============
    constexpr double LOGIC_BASED_CONFLICT_THETA_RANGE = 5;   // 冲突假设航向范围，务必远大于雷达站角度分辨率，不然无减枝作用
    constexpr double LOGIC_BASED_CONFLICT_RHO_RANGE = 5;     // 冲突假设角度范围，务必远大于雷达站距离分辨率，不然无减枝作用

    /*****************************************************************************
     * @brief 状态码
     * 1. 0：处理成功，有航迹输出
     * 2. 1000-1100：霍夫变换处理成功，但推荐更换处理方式
     * 3. 2000-2100：逻辑法处理成功，但推荐更换处理方式
     *****************************************************************************/
    enum class ProcessStatus : std::int32_t
    {
        SUCCESS = 0, // 处理成功，有航迹输出

        // 霍夫变换类状态码
        POINTS_TOO_DISPERSED = 1001, // 点迹过于分散,部分点迹被忽略
        TOO_MANY_CLUSTER = 1002,     // 聚类数量过多，可以更换成其他算法
        TOO_LARGE_CLUSTER = 1003,    // 过大的聚类区域，群目标存在丢失情况，考虑更换参数

        // 逻辑法状态码
        NO_HYPOTHESIS = 2001,      // 没有生成任何假设
        TOO_MANY_HYPOTHSIS = 2002, // 单波门中假设超过上限，有假设因此被舍弃
        WRONG_PREV_NODE = 2003,    // 前置假设节点存在异常

        // 其他
        NO_POINT = 9001, // 没有点迹，无法处理
        //  VELOCITY_OUT_OF_RANGE, // 速度约束不满足
        //  ACCELERATION_TOO_HIGH, // 机动性太强
        //  LOW_CONFIDENCE,        // 置信度过低
        //  AMBIGUOUS_TRACKS,      // 存在多个相似航迹假设

    };
}

#endif // _TRACK_PROJECT_DEF_INIT_H