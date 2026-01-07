/*****************************************************************************
 * @file ManagementService.cpp
 * @brief 管理服务实现类 - 实现文件
 *
 * 实现要求：
 * 1. 对外暴露接口不变（继承自TrackManagementAPI）
 * 2. 内部是一个循环，每次先处理merge，然后处理create，然后处理add指令
 * 3. 初始化的时候生成线程，一直循环
 * 4. 析构的时候释放线程
 *
 * @version 0.1
 * @date 2025-12-05
 *
 * @copyright Copyright (c) 2025
 *****************************************************************************/

#include "../include/ManagementService.hpp"
#include <iostream>
#include <algorithm>

#include "../utils/Logger.hpp"

namespace track_project
{

    /*****************************************************************************
     * @brief 构造函数
     *
     * @param track_size 航迹容量上限
     * @param point_size 点迹容量上限
     *****************************************************************************/
    ManagementService::ManagementService(std::uint32_t track_size, std::uint32_t point_size)
        : tracker_manager_(track_size, point_size),
          track_visualizer_(119.9, 120.1, 29.9, 30.1, track_size, point_size),
          stop_flag_(false)
    {
        // 启动工作线程
        worker_thread_ = std::thread(&ManagementService::worker_thread, this);
        std::cout << "ManagementService: 工作线程已启动" << std::endl;
    }

    /*****************************************************************************
     * @brief 析构函数，停止工作线程并清理资源
     *****************************************************************************/
    ManagementService::~ManagementService()
    {
        // 设置停止标志
        stop_flag_ = true;

        // 通知工作线程
        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            queue_cv_.notify_all();
        }

        // 等待工作线程结束
        if (worker_thread_.joinable())
        {
            worker_thread_.join();
            std::cout << "ManagementService: 工作线程已停止" << std::endl;
        }

        // 清理缓冲区
        {
            std::lock_guard<std::mutex> lock(buffer_mutex_);
            create_buffer_.clear();
            add_buffer_.clear();
            draw_buffer_.clear();
        }
    }

    /*****************************************************************************
     * @brief 航迹生成请求
     *
     * @param new_track 新航迹结构
     *****************************************************************************/
    void ManagementService::create_track_command(std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        std::lock_guard<std::mutex> lock(buffer_mutex_);

        // 复制数据到缓冲区
        create_buffer_ = new_track;

        // 创建指令并加入队列
        Command cmd(CommandType::CREATE);
        cmd.create_data.new_track = &create_buffer_;

        {
            std::lock_guard<std::mutex> queue_lock(queue_mutex_);
            command_queue_.push(cmd);
            queue_cv_.notify_one();
        }

        std::cout << "ManagementService: 创建航迹指令已加入队列，数量: " << new_track.size() << std::endl;
    }

    /*****************************************************************************
     * @brief 航迹添加请求
     *
     * @param updated_track 卡尔曼滤波结果
     *****************************************************************************/
    void ManagementService::add_track_command(std::vector<std::pair<TrackerHeader, TrackPoint>> &updated_track)
    {
        std::lock_guard<std::mutex> lock(buffer_mutex_);

        // 复制数据到缓冲区
        add_buffer_ = updated_track;

        // 创建指令并加入队列
        Command cmd(CommandType::ADD);
        cmd.add_data.updated_track = &add_buffer_;

        {
            std::lock_guard<std::mutex> queue_lock(queue_mutex_);
            command_queue_.push(cmd);
            queue_cv_.notify_one();
        }

        std::cout << "ManagementService: 添加航迹指令已加入队列，数量: " << updated_track.size() << std::endl;
    }

    /*****************************************************************************
     * @brief 航迹融合请求
     *
     * @param source_track_id
     * @param target_track_id
     *****************************************************************************/
    void ManagementService::merge_command(std::uint32_t source_track_id, std::uint32_t target_track_id)
    {
        // 创建指令并加入队列
        Command cmd(CommandType::MERGE);
        cmd.merge_data.source_track_id = source_track_id;
        cmd.merge_data.target_track_id = target_track_id;

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            command_queue_.push(cmd);
            queue_cv_.notify_one();
        }

        std::cout << "ManagementService: 融合航迹指令已加入队列，源ID: "
                  << source_track_id << ", 目标ID: " << target_track_id << std::endl;
    }

    /*****************************************************************************
     * @brief 清空数据区
     *****************************************************************************/
    void ManagementService::clear_all_command()
    {
        // 创建指令并加入队列
        Command cmd(CommandType::CLEAR_ALL);

        {
            std::lock_guard<std::mutex> lock(queue_mutex_);
            command_queue_.push(cmd);
            queue_cv_.notify_one();
        }

        std::cout << "ManagementService: 清空所有指令已加入队列" << std::endl;
    }

    /*****************************************************************************
     * @brief 工作线程函数，循环处理指令
     *
     * 按照优先级顺序处理指令：DRAW -> MERGE -> CREATE -> ADD -> CLEAR_ALL
     * DRAW指令优先级最高，必须优先执行
     * 每个循环迭代中，先处理所有DRAW指令，然后是MERGE，最后是CREATE和ADD
     *****************************************************************************/
    void ManagementService::worker_thread()
    {
        std::cout << "ManagementService: 工作线程开始运行" << std::endl;

        while (!stop_flag_)
        {
            // 按照优先级顺序处理指令
            bool processed = false;

            // 处理所有DRAW指令（优先级最高）
            processed |= process_commands_by_type(CommandType::DRAW);

            // 处理所有MERGE指令
            processed |= process_commands_by_type(CommandType::MERGE);

            // 处理所有CREATE指令
            processed |= process_commands_by_type(CommandType::CREATE);

            // 处理所有ADD指令
            processed |= process_commands_by_type(CommandType::ADD);

            // 处理CLEAR_ALL指令（如果有）
            processed |= process_commands_by_type(CommandType::CLEAR_ALL);

            // 绘制航迹（显示当前状态）
            track_visualizer_.draw_track(tracker_manager_);

            // 如果没有指令处理，等待新指令
            if (!processed && !stop_flag_)
            {
                std::unique_lock<std::mutex> lock(queue_mutex_);
                queue_cv_.wait(lock, [this]()
                               { return !command_queue_.empty() || stop_flag_; });
            }

            // 短暂休眠避免CPU空转
            if (!processed && !stop_flag_)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }

        std::cout << "ManagementService: 工作线程结束运行" << std::endl;
    }

    /*****************************************************************************
     * @brief 处理指定类型的所有指令
     *
     * @param type 要处理的指令类型
     * @return bool 是否处理了任何指令
     *****************************************************************************/
    bool ManagementService::process_commands_by_type(CommandType type)
    {
        bool processed = false;

        while (!stop_flag_)
        {
            Command cmd(type); // 默认值，实际类型会被检查

            // 检查队列中是否有指定类型的指令
            {
                std::lock_guard<std::mutex> lock(queue_mutex_);

                // 查找第一个指定类型的指令
                std::queue<Command> temp_queue;
                bool found = false;

                while (!command_queue_.empty())
                {
                    Command current_cmd = command_queue_.front();
                    command_queue_.pop();

                    if (!found && current_cmd.type == type)
                    {
                        cmd = current_cmd;
                        found = true;
                        processed = true;
                    }
                    else
                    {
                        temp_queue.push(current_cmd);
                    }
                }

                // 将其他指令放回队列
                while (!temp_queue.empty())
                {
                    command_queue_.push(temp_queue.front());
                    temp_queue.pop();
                }

                if (!found)
                {
                    break; // 没有找到指定类型的指令
                }
            }

            // 处理找到的指令
            try
            {
                process_command(cmd);
            }
            catch (const std::exception &e)
            {
                std::cerr << "ManagementService: 处理指令时发生异常: " << e.what() << std::endl;
            }
        }

        return processed;
    }

    /*****************************************************************************
     * @brief 处理单个指令
     *
     * @param cmd 要处理的指令
     *****************************************************************************/
    void ManagementService::process_command(const Command &cmd)
    {
        switch (cmd.type)
        {
        case CommandType::DRAW:
            if (cmd.draw_data.point_data != nullptr)
            {
                process_draw(*cmd.draw_data.point_data);
            }
            break;

        case CommandType::MERGE:
            process_merge(cmd.merge_data.source_track_id, cmd.merge_data.target_track_id);
            break;

        case CommandType::CREATE:
            if (cmd.create_data.new_track != nullptr)
            {
                process_create(*cmd.create_data.new_track);
            }
            break;

        case CommandType::ADD:
            if (cmd.add_data.updated_track != nullptr)
            {
                process_add(*cmd.add_data.updated_track);
            }
            break;

        case CommandType::CLEAR_ALL:
            process_clear_all();
            break;

        default:
            std::cerr << "ManagementService: 未知指令类型" << std::endl;
            break;
        }
    }

    /*****************************************************************************
     * @brief 处理merge指令
     *
     * @param source_track_id 源航迹ID
     * @param target_track_id 目标航迹ID
     *****************************************************************************/
    void ManagementService::process_merge(std::uint32_t source_track_id, std::uint32_t target_track_id)
    {
        LOG_DEBUG << "ManagementService: 处理融合指令，源ID: " << source_track_id << ", 目标ID: " << target_track_id << std::endl;

        bool success = tracker_manager_.merge_tracks(source_track_id, target_track_id);

        if (!success)
        {
            LOG_ERROR << "TODO！这一部分还没做完"
                      << "ManagementService: 处理融合指令将" << source_track_id << "融合到" << target_track_id << "失败" << std::endl;
            // TODO准备发送一些请求，告诉操作者请求失败
        }
    }

    /*****************************************************************************
     * @brief 处理create指令
     *
     * @param new_track 新航迹数据
     *****************************************************************************/
    void ManagementService::process_create(std::vector<std::array<TrackPoint, 4>> &new_track)
    {
        LOG_DEBUG << "ManagementService: 处理创建指令，航迹数量: " << new_track.size() << std::endl;

        for (const auto &track_array : new_track)
        {
            // 为每个航迹数组创建航迹
            std::uint32_t track_id = tracker_manager_.create_track();

            if (track_id == 0)
            {
                LOG_ERROR << "ManagementService: 创建航迹 " << track_id << " 失败:航迹已满" << std::endl;
            }

            // 将点迹添加到航迹中
            for (const auto &point : track_array)
            {
                bool push_success = tracker_manager_.push_track_point(track_id, point);

                if (!push_success) // 添加失败撤回整个流程
                {
                    LOG_ERROR << "ManagementService: 创建航迹 " << track_id << " 失败:unknown" << std::endl;
                    // 如果添加失败，删除该航迹
                    tracker_manager_.delete_track(track_id);
                    return;
                }
            }

            // TODO，向其他线程发送最新一个航迹点
        }
    }

    /*****************************************************************************
     * @brief 处理add指令
     *
     * @param updated_track 更新的航迹数据
     *****************************************************************************/
    void ManagementService::process_add(std::vector<std::pair<TrackerHeader, TrackPoint>> &updated_track)
    {
        LOG_DEBUG << "ManagementService: 处理添加指令，更新数量: " << updated_track.size() << std::endl;

        for (auto &pair : updated_track)
        {
            TrackerHeader &header = pair.first;
            TrackPoint &point = pair.second;

            // 添加点迹到航迹
            bool success = tracker_manager_.push_track_point(header.track_id, point);

            if (!success)
            {
                LOG_ERROR << "ManagementService: 添加点迹到航迹 " << header.track_id << " 失败:unknown" << std::endl;
            }
        }
    }

    /*****************************************************************************
     * @brief 处理clear_all指令
     *****************************************************************************/
    void ManagementService::process_clear_all()
    {
        LOG_INFO << "ManagementService: 全部清空";
        tracker_manager_.clear_all();
        track_visualizer_.clear_all();
    }

    /*****************************************************************************
     * @brief 点迹绘制请求
     *
     * @param point 请求绘制的点迹
     *****************************************************************************/
    void ManagementService::draw_point_command(std::vector<TrackPoint> &point)
    {
        std::lock_guard<std::mutex> lock(buffer_mutex_);

        // 复制数据到缓冲区
        draw_buffer_ = point;

        // 创建指令并加入队列
        Command cmd(CommandType::DRAW);
        cmd.draw_data.point_data = &draw_buffer_;

        {
            std::lock_guard<std::mutex> queue_lock(queue_mutex_);
            command_queue_.push(cmd);
            queue_cv_.notify_one();
        }

        std::cout << "ManagementService: 点迹绘制指令已加入队列，数量: " << point.size() << std::endl;
    }

    /*****************************************************************************
     * @brief 处理draw指令
     *
     * @param point_data 点迹数据
     *****************************************************************************/
    void ManagementService::process_draw(std::vector<TrackPoint> &point_data)
    {
        LOG_DEBUG << "ManagementService: 处理点迹绘制指令，数量: " << point_data.size() << std::endl;
        
        // 调用TrackerVisualizer的draw_point_cloud函数
        track_visualizer_.draw_point_cloud(point_data);
        
    }

} // namespace track_project
