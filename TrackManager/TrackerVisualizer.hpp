/*****************************************************************************
 * @file TrackerVisualizer.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 航迹管理器可视化模块
 * 1. 需求：trackmanager的航迹结构
 * 2. 提供航迹绘制、点迹绘制、清零操作
 * 3. 要求航迹后到达，不然只能显示点迹
 *
 * @version 0.1
 * @date 2025-11-29
 *
 * @copyright Copyright (c) 2025 xjl
 *
 *****************************************************************************/
#ifndef TRACKER_VISUALIZER_HPP
#define TRACKER_VISUALIZER_HPP

#include <opencv2/opencv.hpp>
#include "TrackerManager.hpp"
#include "Logger.hpp"

namespace track_project::trackmanager
{

    class TrackerVisualizer
    {
    public:
        /*****************************************************************************
         * @brief 构造可视化对象
         * @param lon_min 经度最小值
         * @param lon_max 经度最大值
         * @param lat_min 纬度最小值
         * @param lat_max 纬度最大值
         * @param track_size 预分配航迹ID数量
         * @param track_length 预分配航迹点数量
         *****************************************************************************/
        TrackerVisualizer(double lon_min, double lon_max, double lat_min, double lat_max,
                          std::uint32_t track_size = 2000, std::uint32_t track_length = 2000);

        ~TrackerVisualizer() = default;

        /*****************************************************************************
         * @brief 读取航迹管理器并绘制航迹
         * 航迹用渐变黑色线条绘制，新航迹点透明度较高，历史航迹点透明度较低，标注航迹号
         * 当航迹号错误，如内存超过范围、航迹点为0、航迹点超出边界，跳过该航迹，并报相关错误
         *
         * @param manager 航迹管理器对象
         *****************************************************************************/
        void draw_track(const TrackerManager &manager);

        /*****************************************************************************
         * @brief 绘制点云
         *****************************************************************************/
        void draw_point_cloud(std::vector<TrackPoint> x);

        /*****************************************************************************
         * @brief 清楚画布上的所有航迹
         *****************************************************************************/
        void clear_all();

        /*****************************************************************************
         * @brief 打印航迹管理器完整状态到日志，包括统计信息（INFO）和内存池详情（DEBUG）
         *
         * @param manager 航迹管理器对象
         *****************************************************************************/
        void print_full_state(const TrackerManager &manager);

    private:
        // 打印统计信息
        void print_statistics(const TrackerManager &manager, std::stringstream &ss);

        // 打印内存池详情
        void print_memory_pool(const TrackerManager &manager, std::stringstream &ss);

        // 状态转换辅助函数
        const char *state_to_string(int state);

        // 坐标转换：经纬度到图像像素
        cv::Point convert_to_image_coords(double longitude, double latitude) const;

        // 绘制单个航迹
        void draw_single_track(std::uint32_t track_id, const TrackerManager &manager);

        // 绘制航迹线条
        void draw_track_lines(const std::vector<cv::Point> &points);

        // 绘制航迹标签
        void draw_track_label(std::uint32_t track_id, const cv::Point &position);

    private:
        // 画布参数
        cv::Mat img;                               // 画布
        cv::Mat bg_img;                            // 背景图像，存储点迹
        double lon_min, lon_max, lat_min, lat_max; // 经纬度范围
        std::uint32_t height, width;               // 画布高度和宽度

        // 航迹点存放空间,为提高速度采用预分配方式，仅被draw_track使用
        std::vector<std::uint32_t> active_track_ids;
        std::vector<cv::Point> track_points;
    };

} // namespace track_project::trackmanager
#endif // TRACKER_VISUALIZER_HPP
