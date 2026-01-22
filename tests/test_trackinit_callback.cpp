// 回调函数绑定与触发测试
#include <catch2/catch_all.hpp>
#include <thread>
#include <atomic>
#include <csignal>
#include <iostream>

#include "../include/ManagementService.hpp"
#include "../src/HoughSlice.hpp"

using track_project::ManagementService;
using track_project::TrackPoint;
using namespace track_project::trackinit;

namespace track_project::trackinit
{
    class testTrackInitBase
    {
    public:
        static void test_callback(TrackInitBase &obj, std::vector<std::array<TrackPoint, 4>> &new_track)
        {
            if (obj.trackCallback_)
            {
                obj.trackCallback_(new_track);
            }
        }
    };
}
TEST_CASE("回调函数测试", "[basic]")
{
    SECTION("基本回调测试")
    {
        // 创建HoughSlice实例
        SliceHough hough_slice;

        // 创建航迹管理器实例
        ManagementService mgmt_service(10, 50);

        // 绑定回调函数到create_track_command
        hough_slice.set_track_callback(
            [&mgmt_service](std::vector<std::array<TrackPoint, 4>> &new_tracks)
            {
                // 直接调用ManagementService的create_track_command
                mgmt_service.create_track_command(new_tracks);
            });

        // 生成测试点迹数据
        std::vector<std::array<TrackPoint, 4>> mockTracks;
        std::array<TrackPoint, 4> track1; // 测试航迹
        for (int i = 0; i < 4; i++)
        {
            track1[i].longitude = 0.1 + static_cast<float>(i * 0.01); // x坐标
            track1[i].latitude = 0.1 + static_cast<float>(i * 0.01);  // y坐标
            track1[i].doppler = 10.0;                                 // 多普勒速度
        }
        mockTracks.push_back(track1);

        // 强行触发回调函数
        testTrackInitBase::test_callback(hough_slice, mockTracks);

        // 延时，观察结果
        std::this_thread::sleep_for(std::chrono::milliseconds(5000));
    }
}