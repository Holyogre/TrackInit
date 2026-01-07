#include "TrackerVisualizer.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace track_project::trackmanager
{

    TrackerVisualizer::TrackerVisualizer(double lon_min, double lon_max,
                                         double lat_min, double lat_max,
                                         std::uint32_t track_size, std::uint32_t track_length)
        : img(1440, 2560, CV_8UC3, cv::Scalar(255, 255, 255)),    // 彩色RGB画布，背景是白色
          bg_img(1440, 2560, CV_8UC3, cv::Scalar(255, 255, 255)), // 彩色RGB画布，背景是白色，存储背景
          lon_min(lon_min), lon_max(lon_max),
          lat_min(lat_min), lat_max(lat_max)
    {
        height = img.rows;
        width = img.cols;
        active_track_ids.reserve(track_size);
        track_points.reserve(track_length);

        LOG_DEBUG << "TrackerVisualizer初始化完成: 画布"
                  << width << "x" << height
                  << ", 范围[" << lon_min << "," << lon_max << "]x["
                  << lat_min << "," << lat_max << "]" << std::endl;
    }

    void TrackerVisualizer::draw_track(const TrackerManager &manager)
    {
        bg_img.copyTo(img); // 显示点迹结果

        // 使用公开只读接口获取活跃航迹ID
        active_track_ids = manager.get_active_track_ids();

        for (auto track_id : active_track_ids)
        {
            draw_single_track(track_id, manager);
        }

        cv::imshow("Track Visualizer", img);
        cv::waitKey(10);
    }

    void TrackerVisualizer::draw_point_cloud(std::vector<TrackPoint> x)
    {
        // 重置背景为白色
        bg_img.setTo(cv::Scalar(255, 255, 255));

        if (x.empty())
        {
            LOG_ERROR << "TrackerVisualizer: 点迹向量为空，跳过绘制" << std::endl;
            return;
        }

        LOG_DEBUG << "TrackerVisualizer: 开始绘制点迹，数量: " << x.size() << std::endl;

        // 绘制每个点迹
        for (const auto &point : x)
        {
            // 将经纬度坐标转换为图像坐标
            cv::Point img_point = convert_to_image_coords(point.longitude, point.latitude);

            // 检查点是否在图像范围内
            if (img_point.x < 0 || img_point.x >= width ||
                img_point.y < 0 || img_point.y >= height)
            {
                LOG_DEBUG << "TrackerVisualizer: 点迹坐标超出图像范围，跳过绘制 ("
                          << point.longitude << ", " << point.latitude << ")" << std::endl;
                continue;
            }

            // 根据点迹属性选择颜色
            cv::Scalar point_color;
            if (point.is_associated)
            {
                // 已关联的点迹使用蓝色
                point_color = cv::Scalar(255, 0, 0); // BGR格式：蓝色
            }
            else
            {
                // 未关联的点迹使用红色
                point_color = cv::Scalar(0, 0, 255); // BGR格式：红色
            }

            // 在背景图像上绘制点迹（使用小圆点）
            cv::circle(bg_img, img_point, 3, point_color, -1); // 半径为3的实心圆
        }
    }

    cv::Point TrackerVisualizer::convert_to_image_coords(double longitude, double latitude) const
    {
        int x = static_cast<int>(((longitude - lon_min) / (lon_max - lon_min)) * width);
        int y = static_cast<int>(((lat_max - latitude) / (lat_max - lat_min)) * height);
        return cv::Point(x, y);
    }

    void TrackerVisualizer::clear_all()
    {
        // 将画布颜色重置为初始状态（白色背景）
        img.setTo(cv::Scalar(255, 255, 255));

        // 清空显示窗口（如果存在）
        if (cv::getWindowProperty("Track Visualizer", cv::WND_PROP_VISIBLE) >= 0)
        {
            cv::imshow("Track Visualizer", img);
            cv::waitKey(1);
        }

        LOG_DEBUG << "TrackerVisualizer: 画布已清空，重置为初始状态" << std::endl;
    }

    void TrackerVisualizer::draw_single_track(std::uint32_t track_id, const TrackerManager &manager)
    {
        using TrackPoint = track_project::TrackPoint;
        using TrackerHeader = track_project::TrackerHeader;

        // 使用零拷贝只读引用获取头部和数据
        const TrackerHeader *header_ptr = manager.get_header_ref(track_id);
        const LatestKBuffer<TrackPoint> *data_ptr = manager.get_data_ref(track_id);

        if (!header_ptr || !data_ptr)
        {
            LOG_ERROR << "TrackerVisualizer: 无法获取航迹" << track_id << "的只读引用，跳过";
            return;
        }

        if (data_ptr->size() == 0)
        {
            LOG_ERROR << "TrackerVisualizer: 航迹ID" << track_id << "的航迹点为空，跳过该航迹绘制";
            return;
        }

        // 坐标转换，超界点跳过
        track_points.clear();
        for (size_t i = 0; i < data_ptr->size(); ++i)
        {
            const auto &point = (*data_ptr)[i];
            cv::Point img_point = convert_to_image_coords(point.longitude, point.latitude);

            if (img_point.x < 0 || img_point.x >= width ||
                img_point.y < 0 || img_point.y >= height)
            {
                LOG_ERROR << "航迹ID" << track_id << "点" << i << "坐标超出图像范围，跳过该点";
                continue;
            }

            track_points.push_back(img_point);
        }

        if (track_points.size() < 2)
        {
            LOG_ERROR << "航迹ID" << track_id << "有效点少于2个，无法绘制线条";
            return;
        }

        // 绘制航迹线条和标签
        draw_track_lines(track_points);
        draw_track_label(track_id, track_points.back());
    }

    void TrackerVisualizer::draw_track_lines(const std::vector<cv::Point> &points)
    {
        // 绘制航迹线条，使用颜色渐变表示时间
        for (int i = 1; i < points.size(); ++i)
        {
            double age_ratio = static_cast<double>(i) / points.size();
            int gray_value = static_cast<int>(255 * (1.0 - age_ratio));
            cv::Scalar color(gray_value, gray_value, gray_value);
            cv::line(img, points[i - 1], points[i], color, 2);
        }
    }

    void TrackerVisualizer::draw_track_label(std::uint32_t track_id, const cv::Point &position)
    {
        std::string track_text = std::to_string(track_id);

        // 设置文字参数
        const int font_face = cv::FONT_HERSHEY_SIMPLEX;
        const double font_scale = 0.6;
        const int thickness = 2;
        const cv::Scalar text_color(0, 0, 0);

        // 计算文字位置（确保在图像范围内）
        cv::Point text_org(position.x + 5, position.y - 5);

        // 边界检查：确保文字不会超出图像边界
        int text_width = track_text.length() * 12; // 估算文字宽度
        int text_height = 20;                      // 估算文字高度

        if (text_org.x + text_width >= width) // 横坐标判断
        {
            text_org.x = width - text_width - 5;
        }
        else if (text_org.x < 0)
        {
            text_org.x = 5;
        }
        if (text_org.y - text_height < 0) // 纵坐标判断
        {
            text_org.y = text_height + 5;
        }
        else if (text_org.y >= height)
        {
            text_org.y = height - 5;
        }

        // 直接绘制文字
        cv::putText(img, track_text, text_org, font_face, font_scale,
                    text_color, thickness, cv::LINE_AA);
    }

    void TrackerVisualizer::print_full_state(const TrackerManager &manager)
    {
        std::stringstream ss;
        ss << "\n"
           << std::string(60, '=') << std::endl;
        ss << "              TRACKER MANAGER 完整状态" << std::endl;
        ss << std::string(60, '=') << std::endl;

        // 打印航迹管理类统计信息
        print_statistics(manager, ss);
        ss << std::endl;
        LOG_INFO << ss.str();

        // 打印内存池详情
        ss.clear();
        ss.str("");
        print_memory_pool(manager, ss);
        ss << std::endl;

        ss << std::string(60, '=') << "\n"
           << std::endl;

        // 输出到日志
        LOG_DEBUG << ss.str();
    }

    void TrackerVisualizer::print_statistics(const TrackerManager &manager, std::stringstream &ss)
    {
        ss << "系统统计:" << std::endl;
        ss << std::string(50, '-') << std::endl;
        ss << "  总容量: " << manager.get_total_capacity() << " 个航迹" << std::endl;
        ss << "  使用中: " << manager.get_used_count() << " 个航迹" << std::endl;
        ss << "  下个ID: " << manager.get_next_track_id() << std::endl;
    }

    void TrackerVisualizer::print_memory_pool(const TrackerManager &manager, std::stringstream &ss)
    {
        ss << "内存池详情 (" << manager.get_total_capacity() << "个槽位):" << std::endl;
        ss << std::string(50, '-') << std::endl;

        size_t active_count = 0;
        std::vector<std::uint32_t> active_ids = manager.get_active_track_ids();

        for (auto track_id : active_ids)
        {
            const track_project::TrackerHeader *header_ptr = manager.get_header_ref(track_id);
            const LatestKBuffer<track_project::TrackPoint> *data_ptr = manager.get_data_ref(track_id);
            if (!header_ptr || !data_ptr)
                continue;

            active_count++;
            ss << "  航迹" << std::setw(4) << header_ptr->track_id
               << " [状态:" << std::setw(4) << state_to_string(header_ptr->state)
               << ", 外推:" << std::setw(1) << header_ptr->extrapolation_count
               << ", 点数:" << std::setw(3) << data_ptr->size() << "]";

            // 显示最近点的时间（如果有）
            if (data_ptr->size() > 0)
            {
                const auto &latest_point = (*data_ptr)[data_ptr->size() - 1];
                ss << " 最新时间:" << latest_point.time;
            }
            ss << std::endl;
            ss << std::endl;
        }

        if (active_count == 0)
        {
            ss << "  [无活跃航迹]" << std::endl;
        }
    }

    const char *TrackerVisualizer::state_to_string(int state)
    {
        switch (state)
        {
        case 0:
            return "正常";
        case 1:
            return "外推";
        case 2:
            return "终结";
        default:
            return "未知";
        }
    }

} // namespace track_project::trackmanager
