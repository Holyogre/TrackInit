/*****************************************************************************
 * @file ManagementService.hpp
 * @brief 航迹管理服务 - 对外统一接口，内部使用多线程处理指令
 *
 * 主要功能：
 * 1. 航迹创建、更新、融合、删除全生命周期管理
 * 2. 实时航迹与点迹可视化
 * 3. 多线程优先级指令处理
 * 4. 线程安全的数据缓冲区管理
 *
 * @version 1.0
 * @date 2025-12-06
 *
 * @copyright Copyright (c) 2025
 *****************************************************************************/

#ifndef _MANAGEMENT_SERVICE_HPP_
#define _MANAGEMENT_SERVICE_HPP_

#include <memory>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <cstdint>

#include "defstruct.h"
#include "../TrackManager/TrackerManager.hpp"
#include "../TrackManager/TrackerVisualizer.hpp"

namespace track_project
{

    /**
     * @class ManagementService
     * @brief 航迹管理服务主类 - 对外统一接口
     *
     * 使用工作线程循环处理指令，指令处理优先级：DRAW -> MERGE -> CREATE -> ADD -> CLEAR_ALL
     */
    class ManagementService
    {
    public:
        /*****************************************************************************
         * @brief 构造函数
         *
         * @param track_size 航迹容量上限
         * @param point_size 点迹容量上限
         *****************************************************************************/
        ManagementService(std::uint32_t track_size = 2000, std::uint32_t point_size = 2000);

        /*****************************************************************************
         * @brief 析构函数，停止工作线程并清理资源
         *****************************************************************************/
        virtual ~ManagementService();

        // 禁止拷贝和移动
        ManagementService(const ManagementService &) = delete;
        ManagementService &operator=(const ManagementService &) = delete;
        ManagementService(ManagementService &&) = delete;
        ManagementService &operator=(ManagementService &&) = delete;

        /*****************************************************************************
         * @brief 航迹生成请求
         *
         * @param new_track 新航迹结构
         *****************************************************************************/
        void create_track_command(std::vector<std::array<TrackPoint, 4>> &new_track);

        /*****************************************************************************
         * @brief 航迹添加请求
         *
         * @param updated_track 卡尔曼滤波结果
         *****************************************************************************/
        void add_track_command(std::vector<std::pair<TrackerHeader, TrackPoint>> &updated_track);

        /*****************************************************************************
         * @brief 航迹融合请求
         *
         * @param source_track_id
         * @param target_track_id
         *****************************************************************************/
        void merge_command(std::uint32_t source_track_id, std::uint32_t target_track_id);

        /*****************************************************************************
         * @brief 点迹绘制请求
         *
         * @param point 请求绘制的点迹
         *****************************************************************************/
        void draw_point_command(std::vector<TrackPoint> &point);

        /*****************************************************************************
         * @brief 清空数据区
         *****************************************************************************/
        void clear_all_command();

        /*****************************************************************************
         * @brief 获取TrackerManager引用（只读）
         *****************************************************************************/
        const trackmanager::TrackerManager &get_tracker_manager() const { return tracker_manager_; }

    private:
        // 指令类型枚举
        enum class CommandType
        {
            DRAW,     // 点迹绘制指令，优先级最高
            MERGE,    // 航迹融合指令
            CREATE,   // 航迹创建指令
            ADD,      // 航迹添加指令
            CLEAR_ALL // 清空所有指令
        };

        // 指令结构体
        struct Command
        {
            CommandType type;

            // 根据不同指令类型存储不同数据
            union
            {
                struct
                {
                    std::vector<TrackPoint> *point_data; // DRAW指令数据
                } draw_data;

                struct
                {
                    std::uint32_t source_track_id;
                    std::uint32_t target_track_id;
                } merge_data;

                struct
                {
                    std::vector<std::array<TrackPoint, 4>> *new_track;
                } create_data;

                struct
                {
                    std::vector<std::pair<TrackerHeader, TrackPoint>> *updated_track;
                } add_data;
            };

            // 构造函数
            Command(CommandType t) : type(t) {}
        };

        /*****************************************************************************
         * @brief 工作线程函数，循环处理指令
         *****************************************************************************/
        void worker_thread();

        /*****************************************************************************
         * @brief 处理单个指令
         *
         * @param cmd 要处理的指令
         *****************************************************************************/
        void process_command(const Command &cmd);

        /*****************************************************************************
         * @brief 处理merge指令
         *
         * @param source_track_id 源航迹ID
         * @param target_track_id 目标航迹ID
         *****************************************************************************/
        void process_merge(std::uint32_t source_track_id, std::uint32_t target_track_id);

        /*****************************************************************************
         * @brief 处理create指令
         *
         * @param new_track 新航迹数据
         *****************************************************************************/
        void process_create(std::vector<std::array<TrackPoint, 4>> &new_track);

        /*****************************************************************************
         * @brief 处理add指令
         *
         * @param updated_track 更新的航迹数据
         *****************************************************************************/
        void process_add(std::vector<std::pair<TrackerHeader, TrackPoint>> &updated_track);

        /*****************************************************************************
         * @brief 处理clear_all指令
         *****************************************************************************/
        void process_clear_all();

        /*****************************************************************************
         * @brief 处理draw指令
         *
         * @param point_data 点迹数据
         *****************************************************************************/
        void process_draw(std::vector<TrackPoint> &point_data);

        /*****************************************************************************
         * @brief 处理指定类型的所有指令
         *
         * @param type 要处理的指令类型
         * @return bool 是否处理了任何指令
         *****************************************************************************/
        bool process_commands_by_type(CommandType type);

    private:
        // Tracker管理器
        trackmanager::TrackerManager tracker_manager_;
        trackmanager::TrackerVisualizer track_visualizer_;

        // 线程控制
        std::thread worker_thread_;
        std::atomic<bool> stop_flag_;

        // 指令队列
        std::queue<Command> command_queue_;
        std::mutex queue_mutex_;
        std::condition_variable queue_cv_;

        // 数据存储（用于避免数据在队列中失效）
        std::vector<std::array<TrackPoint, 4>> create_buffer_;
        std::vector<std::pair<TrackerHeader, TrackPoint>> add_buffer_;
        std::vector<TrackPoint> draw_buffer_; // DRAW指令数据缓冲区
        std::mutex buffer_mutex_;

        // 用于执行对应操作的回调函数
        std::function<void()> merge_result_callback; // 航迹融合操作结果反馈函数
    };

} // namespace track_project

#endif // _MANAGEMENT_SERVICE_HPP_
