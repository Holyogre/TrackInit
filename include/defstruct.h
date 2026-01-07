/*****************************************************************************
 * @file defstruct.h
 * @author xjl (xjl20011009@126.com)
 *
 * @brief 包含三个部分
 * 1. 通用结构体
 * 2，流水线架构专用结构体
 * 3.通信模块结构体
 *
 * @version 0.1
 * @date 2025-11-03
 *
 * @copyright Copyright (c) 2025
 *****************************************************************************/
#ifndef _DEF_TRACK_STRUCT_H_
#define _DEF_TRACK_STRUCT_H_

// 基础结构
#include <vector>
#include <string>
#include <cstdint> //精确位数整数
// 文件控制
#include <fstream>
// 输出格式控制
#include <iostream>
#include <iomanip>
// 时间戳
#include <ctime>
#include <chrono>
// 结构体控制
#include <type_traits>
// 断言
#include <assert.h>

namespace track_project
{
    // 时间戳结构（毫秒精度）
    struct Timestamp
    {
        int64_t milliseconds; // 从1970-01-01开始的毫秒数

        Timestamp()
        {
            auto now = std::chrono::system_clock::now();
            milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
                               now.time_since_epoch())
                               .count();
        }

        // 静态函数：获取当前时间戳
        static Timestamp now()
        {
            return Timestamp(); // 直接调用构造函数创建当前时间戳
        }

        // 重载 << 操作符用于调试输出
        friend std::ostream &operator<<(std::ostream &os, const Timestamp &ts)
        {
            time_t seconds = ts.milliseconds / 1000;
            int64_t ms = ts.milliseconds % 1000;

            std::tm tm_info = *std::localtime(&seconds);
            os << std::put_time(&tm_info, "%Y-%m-%d %H:%M:%S")
               << "." << std::setfill('0') << std::setw(3) << ms;
            return os;
        }
    };

    // 通用联合体，用于不同类型转换
    union IntFloatUnion
    {
        int ri;
        float rf;
    };

    /*****************************************************************************
     * @brief 流水线架构：TrackingBuffer，计划设计为环形缓冲区结构体
     * 特性：
     * 1. 表示航迹，点迹大小情况统一使用std::uint32_t确保数据长度一致性和非负性
     * 2. 运动类参数统一使用double
     * 3. is类型参数统一使用bool
     *****************************************************************************/
    struct TrackerHeader
    {
        std::uint32_t track_id = 0;
        std::uint32_t extrapolation_count = 0;
        std::uint32_t point_num = 0;
        int state = 3; // 0表示航迹正常，1表示航迹外推，2表示航迹终结

        void start(std::uint32_t id)
        {
            track_id = id;
            extrapolation_count = 0;
            point_num = 0;
            state = 0;
        };

        void clear()
        {
            track_id = 0;
            extrapolation_count = 0;
            point_num = 0;
            state = -1;
        };
    };

    // 基础点迹结构体，构造的时候统一初始化为0
    struct TrackPoint
    {
        // // 对照信息
        // double angle;    // 雷达观测角度，雷达法线顺时针方向，度
        // double distance; // 雷达观测距离，距离雷达站距离，km

        // 雷达点迹基础信息，传入前确定
        double longitude = 0.0; // 经度，单位度
        double latitude = 0.0;  // 纬度，单位度
        double doppler = 0.0;   // 径向速度，m/s

        // 有AIS数据引入后的点迹信息
        double sog = 0.0; // 对地速度,m/s
        double cog = 0.0; // 对地航向,北偏东角度，度

        // 内部信息，传入后计算得到
        double x = 0.0;  // 单位km
        double y = 0.0;  // 单位km
        double vx = 0.0; // 单位m/s
        double vy = 0.0; // 单位m/s

        // 航迹关联
        bool is_associated = false; // 是否关联

        // 时间戳
        Timestamp time;
    };

} // track_project

#endif
