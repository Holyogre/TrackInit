/*****************************************************************************
 * @file TrackAlgorithmBase.hpp
 * @author xjl (xjl20011009@126.com)
 *
 * @brief 航迹起始算法抽象基类
 * 用于约束readme中三种算法的实现：
 * 1. 基础算法 (M/N逻辑确认、最小二乘状态估计、物理约束过滤)
 * 2. 霍夫变换家族 (标准霍夫变换、概率霍夫变换、约束霍夫变换)
 * 3. RANSAC系列 (基本RANSAC、多模型RANSAC、自适应RANSAC)
 *
 * @version 0.1
 * @date 2025-12-10
 *
 * @copyright Copyright (c) 2025
 *****************************************************************************/
#ifndef _TRACK_INIT_BASE_HPP_
#define _TRACK_INIT_BASE_HPP_

#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <functional>
#include "../include/defstruct.h"
#include "../include/def_init.h"

namespace track_project::trackinit
{
    // 定义回调类型
    using TrackCallback = std::function<void(std::vector<std::array<TrackPoint, 4>> &)>;

    /**
     * @brief 航迹起始算法抽象基类
     *
     * 输入：std::vector<TrackPoint> - 原始点迹数据
     * 输出：std::vector<std::array<TrackPoint, 4>> &new_track - 新生成的航迹
     *
     * 每个航迹由4个TrackPoint组成，表示一个完整的航迹段。
     */
    class TrackInitBase
    {
    public:
        virtual ~TrackInitBase() = default;

        /**
         * @brief 处理点迹数据，生成新航迹
         *
         * @param points 输入点迹数据，包含所有待处理的TrackPoint
         * @param new_track 输出参数，存储生成的新航迹
         * @return ErrorCode 错误码，SUCCESS表示成功，其他值表示具体错误
         */
        virtual ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track) = 0;

        /**
         * @brief 重置算法状态
         *
         * 当config里面请求的算法不变的时候，调用clear_all，不然，调用析构函数
         */
        virtual void clear_all() = 0;

        /**
         * @brief 获取算法名称（用于日志和调试）
         *
         * @return std::string 算法名称
         */
        virtual std::string get_name() const = 0;

    protected:
        TrackCallback trackCallback_; // 用于发送航迹数据
    };
} // namespace track_project::track_init

#endif // _TRACK_INIT_BASE_HPP_
