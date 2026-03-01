/*****************************************************************************
 * @file LogicBased.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 特点，逻辑法，我计划引入一点MHT的思路，这就要求初始航迹的结构可能需要调整
 * @version 0.1
 * @date 2026-03-01
 *
 * @copyright Copyright (c) 2026
 *
 *****************************************************************************/

#ifndef _LOGIC_BASED_HPP_
#define _LOGIC_BASED_HPP_

#include "TrackInitBase.hpp"

namespace track_project::trackinit
{
    class LogicBasedTracker : public TrackInitBase
    {
    private:
        // 宏定义，结构体定义
    public:
        LogicBasedTracker();
        virtual ~LogicBasedTracker() noexcept = default;
        
        /*****************************************************************************
         * @brief 检测直线，在cpp文件中
         * 当点迹过于分散，开辟了过多的空间，就反馈点迹过多
         *****************************************************************************/
        ProcessStatus process(const std::vector<TrackPoint> &points, std::vector<std::array<TrackPoint, 4>> &new_track) override;

        /*****************************************************************************
         * @brief 重置算法状态，清空所有数据
         *****************************************************************************/
        void clear_all() override;

        /*****************************************************************************
         * @brief 获取 name 对象
         *****************************************************************************/
        std::string get_name() const override { return "LogicBasedTracker"; }

    private:
        // 写函数
    };
}
#endif //_LOGIC_BASED_HPP_