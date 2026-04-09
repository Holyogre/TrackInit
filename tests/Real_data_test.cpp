#include <catch2/catch_all.hpp>

#include "../include/ManagementService.hpp"
#include "../include/def_init.h"
#include "../include/defstruct.h"
#include "../src/LogicBasedInitiator.hpp"
#include "../utils/Logger.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace track_project;
using track_project::trackinit::LogicBasedInitiator;

namespace
{
	constexpr double RADAR_LATITUDE_DEG = 37.0;
	constexpr double RADAR_LONGITUDE_DEG = 120.0;
	constexpr double KM_PER_DEGREE = 111.0;
	constexpr double UNIFIED_SIGMA_X_KM = 1.0;
	constexpr double UNIFIED_SIGMA_Y_KM = 1.0;
	constexpr double MIN_RANGE_KM = 50.0;
	constexpr double MAX_RANGE_KM = 300.0;
	constexpr double MAX_ABS_ANGLE_DEG = 30.0;
	constexpr auto PLAYBACK_INTERVAL = std::chrono::seconds(1);
	constexpr auto FINAL_FRAME_HOLD = std::chrono::seconds(300);
	constexpr std::array<int64_t, 29> RADAR_TIMES_S = {
		21, 23, 25, 27, 30, 32, 34, 36, 38, 40,
		42, 44, 46, 48, 50, 52, 54, 57, 59, 61,
		63, 65, 67, 69, 71, 73, 75, 77, 79};

	struct RealDataFrame
	{
		std::filesystem::path source_file;
		std::vector<TrackPoint> points;
		int64_t timestamp_ms = 0;
	};

	struct GeoBounds
	{
		double lon_min = std::numeric_limits<double>::max();
		double lon_max = std::numeric_limits<double>::lowest();
		double lat_min = std::numeric_limits<double>::max();
		double lat_max = std::numeric_limits<double>::lowest();

		void include(const TrackPoint &point)
		{
			lon_min = std::min(lon_min, point.longitude);
			lon_max = std::max(lon_max, point.longitude);
			lat_min = std::min(lat_min, point.latitude);
			lat_max = std::max(lat_max, point.latitude);
		}

		void expand(double margin_deg)
		{
			lon_min -= margin_deg;
			lon_max += margin_deg;
			lat_min -= margin_deg;
			lat_max += margin_deg;
		}
	};

	double deg_to_rad(double degrees)
	{
		return degrees * M_PI / 180.0;
	}

	int extract_sheet_index(const std::filesystem::path &file_path)
	{
		static const std::regex pattern(R"(sheet(\d+)_traceDatShip_\d+\.csv)");
		std::smatch match;
		const std::string filename = file_path.filename().string();
		if (!std::regex_search(filename, match, pattern))
		{
			throw std::runtime_error("无法从文件名中解析 sheet 序号: " + filename);
		}

		return std::stoi(match[1].str());
	}

	std::vector<std::filesystem::path> collect_real_data_files(const std::filesystem::path &data_dir)
	{
		std::vector<std::filesystem::path> files;
		for (const auto &entry : std::filesystem::directory_iterator(data_dir))
		{
			if (!entry.is_regular_file())
			{
				continue;
			}

			const auto filename = entry.path().filename().string();
			if (filename.find("radar_sheet") == std::string::npos)
			{
				continue;
			}

			files.push_back(entry.path());
		}

		std::sort(files.begin(), files.end(), [](const auto &lhs, const auto &rhs)
				  { return extract_sheet_index(lhs) < extract_sheet_index(rhs); });
		return files;
	}

	std::vector<int64_t> build_frame_timestamps_ms(size_t frame_count)
	{
		std::vector<int64_t> timestamps_ms;
		timestamps_ms.reserve(frame_count);

		if (frame_count == RADAR_TIMES_S.size())
		{
			for (int64_t time_s : RADAR_TIMES_S)
			{
				timestamps_ms.push_back(time_s * 1000);
			}
			return timestamps_ms;
		}

		if (frame_count == RADAR_TIMES_S.size() + 1)
		{
			const int64_t first_gap_s = RADAR_TIMES_S[1] - RADAR_TIMES_S[0];
			timestamps_ms.push_back((RADAR_TIMES_S[0] - first_gap_s) * 1000);
			for (int64_t time_s : RADAR_TIMES_S)
			{
				timestamps_ms.push_back(time_s * 1000);
			}
			return timestamps_ms;
		}

		throw std::runtime_error("实测数据帧数与 radar_time 长度不匹配");
	}

	TrackPoint make_track_point(double range_km, double angle_deg, double doppler_m_s, int64_t timestamp_ms)
	{
		TrackPoint point;
		const double angle_rad = deg_to_rad(angle_deg);
		const double radar_lat_rad = deg_to_rad(RADAR_LATITUDE_DEG);

		point.x = range_km * std::sin(angle_rad);
		point.y = range_km * std::cos(angle_rad);
		point.doppler = doppler_m_s;
		point.time.milliseconds = timestamp_ms;

		point.longitude = RADAR_LONGITUDE_DEG + point.x / (KM_PER_DEGREE * std::cos(radar_lat_rad));
		point.latitude = RADAR_LATITUDE_DEG + point.y / KM_PER_DEGREE;

		return point;
	}

	bool should_keep_detection(double range_km, double angle_deg)
	{
		if (std::abs(angle_deg) > MAX_ABS_ANGLE_DEG)
		{
			return false;
		}

		if (range_km < MIN_RANGE_KM || range_km > MAX_RANGE_KM)
		{
			return false;
		}

		return true;
	}

	std::vector<TrackPoint> load_points_from_csv(const std::filesystem::path &csv_path, int64_t timestamp_ms)
	{
		std::ifstream input(csv_path);
		if (!input.is_open())
		{
			throw std::runtime_error("无法打开实测数据文件: " + csv_path.string());
		}

		std::vector<TrackPoint> points;
		points.reserve(512);

		std::string line;
		std::getline(input, line);

		while (std::getline(input, line))
		{
			if (line.empty())
			{
				continue;
			}

			std::stringstream line_stream(line);
			std::string range_token;
			std::string angle_token;
			std::string doppler_token;

			if (!std::getline(line_stream, range_token, ',') ||
				!std::getline(line_stream, angle_token, ',') ||
				!std::getline(line_stream, doppler_token, ','))
			{
				continue;
			}

			const double range_km = std::stod(range_token);
			const double angle_deg = std::stod(angle_token);
			const double doppler_m_s = std::stod(doppler_token);

			if (!std::isfinite(range_km) || !std::isfinite(angle_deg) || !std::isfinite(doppler_m_s))
			{
				continue;
			}

			if (range_km <= 0.0)
			{
				continue;
			}

			if (!should_keep_detection(range_km, angle_deg))
			{
				continue;
			}

			TrackPoint point = make_track_point(range_km, angle_deg, doppler_m_s, timestamp_ms);
			if (std::abs(point.x) > track_project::trackinit::LOGIC_BASED_MAX_ABS_X ||
				std::abs(point.y) > track_project::trackinit::LOGIC_BASED_MAX_ABS_Y)
			{
				continue;
			}

			points.push_back(point);
		}

		return points;
	}

	std::vector<RealDataFrame> load_real_data_frames(const std::filesystem::path &data_dir)
	{
		const auto files = collect_real_data_files(data_dir);
		const auto timestamps_ms = build_frame_timestamps_ms(files.size());

		std::vector<RealDataFrame> frames;
		frames.reserve(files.size());

		for (size_t index = 0; index < files.size(); ++index)
		{
			RealDataFrame frame;
			frame.source_file = files[index];
			frame.timestamp_ms = timestamps_ms[index];
			frame.points = load_points_from_csv(files[index], timestamps_ms[index]);
			frames.push_back(std::move(frame));
		}

		return frames;
	}

	GeoBounds compute_bounds(const std::vector<RealDataFrame> &frames)
	{
		GeoBounds bounds;
		for (const auto &frame : frames)
		{
			for (const auto &point : frame.points)
			{
				bounds.include(point);
			}
		}

		bounds.expand(0.05);
		return bounds;
	}

	bool enable_visualization()
	{
		if (std::getenv("DISPLAY") == nullptr && std::getenv("WAYLAND_DISPLAY") == nullptr)
		{
			return false;
		}

		const char *visualization_flag = std::getenv("TRACKINIT_REALDATA_VIS");
		return visualization_flag == nullptr || std::string(visualization_flag) != "0";
	}
}

namespace track_project::trackinit
{
	class test_LogicBasedInitiator
	{
	public:
		explicit test_LogicBasedInitiator(LogicBasedInitiator &initiator) : initiator_(initiator) {}

		void printHypothesisDistribution() const
		{
			const auto &index = initiator_.current_hypothesis_index_;

			LOG_INFO << "========== 当前假设节点分布 ==========";
			for (size_t x = 0; x < LOGIC_BASED_NUM_X_BINS; ++x)
			{
				for (size_t y = 0; y < LOGIC_BASED_NUM_Y_BINS; ++y)
				{
					const size_t bin_index = x * LOGIC_BASED_NUM_Y_BINS + y;
					const auto &nodes = index[bin_index];
					if (nodes.empty())
					{
						continue;
					}

					LOG_INFO << "位置 [" << x << "," << y << "] (Bin Index: " << bin_index
							 << ")，节点数量：" << nodes.size();
				}
			}
			LOG_INFO << "=====================================";
		}

	private:
		LogicBasedInitiator &initiator_;
	};
}

using track_project::trackinit::test_LogicBasedInitiator;

TEST_CASE("实测数据测试", "[FunctionalityCheck][real_data]")
{
	const std::filesystem::path repo_root = std::filesystem::path(__FILE__).parent_path().parent_path();
	const std::filesystem::path data_dir = repo_root / "data";

	REQUIRE(std::filesystem::exists(data_dir));

	auto frames = load_real_data_frames(data_dir);
	REQUIRE(frames.size() >= 4);

	size_t total_input_points = 0;
	for (const auto &frame : frames)
	{
		total_input_points += frame.points.size();
	}

	const auto bounds = compute_bounds(frames);
	const bool visualization_enabled = enable_visualization();

	LOG_INFO << "实测数据帧数: " << frames.size();
	LOG_INFO << "总点迹数: " << total_input_points;
	LOG_INFO << "雷达站位置(lat, lon): (" << RADAR_LATITUDE_DEG << ", " << RADAR_LONGITUDE_DEG << ")";
	LOG_INFO << "可视化开关: " << (visualization_enabled ? "开启" : "关闭");
	LOG_INFO << "可视化回放间隔: " << PLAYBACK_INTERVAL.count() << " s";

	std::unique_ptr<track_project::ManagementService> track_manager;
	if (visualization_enabled)
	{
		track_manager = std::make_unique<track_project::ManagementService>(
			bounds.lon_min, bounds.lon_max, bounds.lat_min, bounds.lat_max);
		track_manager->clear_all_command();
	}

	LogicBasedInitiator initiator;
	test_LogicBasedInitiator tester(initiator);
	std::vector<std::pair<double, double>> unified_error_table(
		track_project::trackinit::LOGIC_BASED_NUM_X_BINS * track_project::trackinit::LOGIC_BASED_NUM_Y_BINS,
		std::make_pair(UNIFIED_SIGMA_X_KM, UNIFIED_SIGMA_Y_KM));
	initiator.build_error_distribution_table(unified_error_table);
	LOG_INFO << "过滤范围: |A| <= " << MAX_ABS_ANGLE_DEG << " deg, R in ["
			 << MIN_RANGE_KM << ", " << MAX_RANGE_KM << "] km";
	LOG_INFO << "统一误差表: sigma_x=" << UNIFIED_SIGMA_X_KM << " km, sigma_y=" << UNIFIED_SIGMA_Y_KM << " km";

	initiator.set_track_callback([&track_manager](const std::vector<std::array<TrackPoint, 4>> &tracks)
								 {
		LOG_INFO << "回调函数被调用，生成了 " << tracks.size() << " 条航迹";
		if (track_manager)
		{
			track_manager->create_track_command(const_cast<std::vector<std::array<TrackPoint, 4>> &>(tracks));
		} });

	std::vector<std::array<TrackPoint, 4>> new_tracks;
	size_t frames_with_tracks = 0;
	size_t max_track_count = 0;
	size_t total_generated_tracks = 0;

    size_t frame_index_max = frames.size();
    frame_index_max=4;
	for (size_t frame_index = 0; frame_index < frame_index_max; ++frame_index)
	{
		auto &frame = frames[frame_index];
		INFO("frame=" << frame_index + 1 << ", file=" << frame.source_file.filename().string()
						<< ", points=" << frame.points.size() << ", timestamp_ms=" << frame.timestamp_ms);

		LOG_INFO << "处理第 " << frame_index + 1 << " 帧: " << frame.source_file.filename().string()
				 << "，点迹数=" << frame.points.size()
				 << "，时间戳=" << frame.timestamp_ms << " ms";

		if (track_manager)
		{
			track_manager->draw_point_command(frame.points);
		}

		const auto status = initiator.process(frame.points, new_tracks);
		REQUIRE(status == track_project::trackinit::ProcessStatus::SUCCESS);

		if (!new_tracks.empty())
		{
			++frames_with_tracks;
		}
		max_track_count = std::max(max_track_count, new_tracks.size());
		total_generated_tracks += new_tracks.size();

		LOG_INFO << "第 " << frame_index + 1 << " 帧输出航迹数: " << new_tracks.size();

		if (track_manager)
		{
			std::this_thread::sleep_for(PLAYBACK_INTERVAL);
		}
	}

	if (track_manager)
	{
		LOG_INFO << "回放结束，保留最后一帧 " << FINAL_FRAME_HOLD.count() << " s";
		std::this_thread::sleep_for(FINAL_FRAME_HOLD);
	}

	tester.printHypothesisDistribution();

	LOG_INFO << "========== 实测数据统计结果 ==========";
	LOG_INFO << "有效输入帧数: " << frames.size();
	LOG_INFO << "累计输入点迹数: " << total_input_points;
	LOG_INFO << "产生航迹的帧数: " << frames_with_tracks;
	LOG_INFO << "累计输出航迹数: " << total_generated_tracks;
	LOG_INFO << "单帧最大输出航迹数: " << max_track_count;

	CHECK(frames_with_tracks > 0);
	CHECK(max_track_count > 0);
}
