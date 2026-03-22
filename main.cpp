#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "include/defsystem.h"
#include "include/defstruct.h"
#include "utils/Logger.hpp"
#include "tests/Func_noise_pointgen.hpp"
#include "tests/Func_pointgen.hpp"

namespace
{
    using track_project::TrackPoint;

    constexpr double TIME_INTERVAL_S = 10.0;
    constexpr std::size_t LOGIC_WINDOW = 5;
    constexpr std::size_t LOGIC_CONFIRM_HITS = 4;
    constexpr std::size_t MAX_BRANCH_PER_FRAME = 6;

    struct PolarObservation
    {
        const TrackPoint *point = nullptr;
        std::size_t point_index = 0;
        double theta_deg = 0.0;
        double rho_km = 0.0;
    };

    struct TrackCandidate
    {
        std::array<const TrackPoint *, LOGIC_WINDOW> points{};
        std::array<int, LOGIC_WINDOW> point_indices{};
        std::size_t hits = 0;
        double score = 0.0;
    };

    double rad_to_deg(double angle_rad)
    {
        return angle_rad * 180.0 / M_PI;
    }

    double normalize_angle_deg(double angle_deg)
    {
        double normalized = std::fmod(angle_deg, 360.0);
        if (normalized < 0.0)
        {
            normalized += 360.0;
        }
        return normalized;
    }

    double bearing_deg(const TrackPoint &point)
    {
        return normalize_angle_deg(rad_to_deg(std::atan2(point.x, point.y)));
    }

    double angle_diff_deg(double lhs, double rhs)
    {
        double diff = std::abs(normalize_angle_deg(lhs) - normalize_angle_deg(rhs));
        return std::min(diff, 360.0 - diff);
    }

    double dynamic_theta_gate_deg(const PolarObservation &lhs, const PolarObservation &rhs)
    {
        const double max_travel_km = track_project::VELOCITY_MAX * TIME_INTERVAL_S / 1000.0;
        const double reference_rho = std::max(1.0, std::min(lhs.rho_km, rhs.rho_km));
        return track_project::BASE_THETA_SIGMA_DEG + rad_to_deg(std::atan2(max_travel_km, reference_rho));
    }

    double dynamic_rho_gate_km()
    {
        return track_project::BASE_R_SIGMA_KM + track_project::VELOCITY_MAX * TIME_INTERVAL_S / 1000.0;
    }

    double pair_cost(const PolarObservation &lhs, const PolarObservation &rhs)
    {
        const double theta_gate = dynamic_theta_gate_deg(lhs, rhs);
        const double rho_gate = dynamic_rho_gate_km();
        const double theta_cost = angle_diff_deg(lhs.theta_deg, rhs.theta_deg) / std::max(theta_gate, 1e-6);
        const double rho_cost = std::abs(lhs.rho_km - rhs.rho_km) / std::max(rho_gate, 1e-6);
        return theta_cost + rho_cost;
    }

    bool direct_association_ok(const PolarObservation &lhs, const PolarObservation &rhs)
    {
        return angle_diff_deg(lhs.theta_deg, rhs.theta_deg) <= dynamic_theta_gate_deg(lhs, rhs) &&
               std::abs(lhs.rho_km - rhs.rho_km) <= dynamic_rho_gate_km();
    }

    std::vector<PolarObservation> build_observations(const std::vector<TrackPoint> &frame)
    {
        std::vector<PolarObservation> observations;
        observations.reserve(frame.size());

        for (std::size_t index = 0; index < frame.size(); ++index)
        {
            const TrackPoint &point = frame[index];
            observations.push_back(PolarObservation{
                .point = &point,
                .point_index = index,
                .theta_deg = bearing_deg(point),
                .rho_km = std::hypot(point.x, point.y)});
        }

        return observations;
    }

    std::string candidate_key(const TrackCandidate &candidate)
    {
        std::ostringstream oss;
        for (std::size_t frame_index = 0; frame_index < LOGIC_WINDOW; ++frame_index)
        {
            if (frame_index > 0)
            {
                oss << '|';
            }
            oss << candidate.point_indices[frame_index];
        }
        return oss.str();
    }

    void enumerate_candidates(
        const std::vector<std::vector<PolarObservation>> &observations_by_frame,
        std::size_t next_frame_index,
        const PolarObservation &last_observation,
        std::size_t miss_count,
        TrackCandidate &current_candidate,
        std::vector<TrackCandidate> &raw_candidates)
    {
        if (next_frame_index == LOGIC_WINDOW)
        {
            if (current_candidate.hits >= LOGIC_CONFIRM_HITS)
            {
                raw_candidates.push_back(current_candidate);
            }
            return;
        }

        std::vector<std::pair<double, const PolarObservation *>> matches;
        matches.reserve(observations_by_frame[next_frame_index].size());
        for (const auto &observation : observations_by_frame[next_frame_index])
        {
            if (direct_association_ok(last_observation, observation))
            {
                matches.emplace_back(pair_cost(last_observation, observation), &observation);
            }
        }

        std::sort(matches.begin(), matches.end(), [](const auto &lhs, const auto &rhs)
                  { return lhs.first < rhs.first; });

        const std::size_t branch_limit = std::min(matches.size(), MAX_BRANCH_PER_FRAME);
        for (std::size_t match_index = 0; match_index < branch_limit; ++match_index)
        {
            const auto *matched_observation = matches[match_index].second;
            current_candidate.points[next_frame_index] = matched_observation->point;
            current_candidate.point_indices[next_frame_index] = static_cast<int>(matched_observation->point_index);
            current_candidate.hits += 1;
            current_candidate.score += matches[match_index].first;

            enumerate_candidates(observations_by_frame,
                                 next_frame_index + 1,
                                 *matched_observation,
                                 miss_count,
                                 current_candidate,
                                 raw_candidates);

            current_candidate.score -= matches[match_index].first;
            current_candidate.hits -= 1;
            current_candidate.points[next_frame_index] = nullptr;
            current_candidate.point_indices[next_frame_index] = -1;
        }

        if (miss_count + 1 <= (LOGIC_WINDOW - LOGIC_CONFIRM_HITS))
        {
            enumerate_candidates(observations_by_frame,
                                 next_frame_index + 1,
                                 last_observation,
                                 miss_count + 1,
                                 current_candidate,
                                 raw_candidates);
        }
    }

    std::vector<TrackCandidate> simple_logic_initiation(const std::vector<std::vector<TrackPoint>> &frames)
    {
        std::vector<std::vector<PolarObservation>> observations_by_frame;
        observations_by_frame.reserve(frames.size());
        for (const auto &frame : frames)
        {
            observations_by_frame.push_back(build_observations(frame));
        }

        std::vector<TrackCandidate> raw_candidates;
        for (const auto &seed_observation : observations_by_frame.front())
        {
            TrackCandidate candidate;
            candidate.point_indices.fill(-1);
            candidate.points.fill(nullptr);
            candidate.points[0] = seed_observation.point;
            candidate.point_indices[0] = static_cast<int>(seed_observation.point_index);
            candidate.hits = 1;
            enumerate_candidates(observations_by_frame, 1, seed_observation, 0, candidate, raw_candidates);
        }

        std::sort(raw_candidates.begin(), raw_candidates.end(), [](const TrackCandidate &lhs, const TrackCandidate &rhs)
                  {
                      if (lhs.hits != rhs.hits)
                      {
                          return lhs.hits > rhs.hits;
                      }
                      return lhs.score < rhs.score;
                  });

        std::set<std::string> seen_keys;
        std::vector<std::vector<bool>> used(frames.size());
        for (std::size_t frame_index = 0; frame_index < frames.size(); ++frame_index)
        {
            used[frame_index].assign(frames[frame_index].size(), false);
        }

        std::vector<TrackCandidate> selected_tracks;
        for (const auto &candidate : raw_candidates)
        {
            const std::string key = candidate_key(candidate);
            if (!seen_keys.insert(key).second)
            {
                continue;
            }

            bool conflict = false;
            for (std::size_t frame_index = 0; frame_index < frames.size(); ++frame_index)
            {
                const int point_index = candidate.point_indices[frame_index];
                if (point_index >= 0 && used[frame_index][static_cast<std::size_t>(point_index)])
                {
                    conflict = true;
                    break;
                }
            }

            if (conflict)
            {
                continue;
            }

            selected_tracks.push_back(candidate);
            for (std::size_t frame_index = 0; frame_index < frames.size(); ++frame_index)
            {
                const int point_index = candidate.point_indices[frame_index];
                if (point_index >= 0)
                {
                    used[frame_index][static_cast<std::size_t>(point_index)] = true;
                }
            }
        }

        return selected_tracks;
    }

    bool same_label(const TrackPoint &lhs, const TrackPoint &rhs)
    {
        return std::abs(lhs.sog - rhs.sog) < 1e-4 && std::abs(lhs.cog - rhs.cog) < 1e-4;
    }

    void print_track_summary(const std::vector<TrackCandidate> &tracks)
    {
        for (std::size_t track_index = 0; track_index < tracks.size(); ++track_index)
        {
            const auto &track = tracks[track_index];
            // LOG_INFO << "候选航迹[" << track_index << "] hits=" << track.hits << ", score=" << track.score;
            for (std::size_t frame_index = 0; frame_index < LOGIC_WINDOW; ++frame_index)
            {
                const TrackPoint *point = track.points[frame_index];
                if (point == nullptr)
                {
                    LOG_INFO << "  frame=" << frame_index << " miss";
                    continue;
                }

                // LOG_INFO << "  frame=" << frame_index
                //          << " idx=" << track.point_indices[frame_index]
                //          << " theta=" << std::fixed << std::setprecision(2) << bearing_deg(*point)
                //          << " rho=" << std::fixed << std::setprecision(2) << std::hypot(point->x, point->y)
                //          << " sog=" << point->sog
                //          << " cog=" << point->cog;
            }
        }
    }
} // namespace

int main()
{
    using namespace track_project;

    srand(time(nullptr));
    unsigned int kFixedSeed = rand();
    // kFixedSeed = 1343091126U; // 固定随机种子，保证结果可复现
    const std::array<unsigned int, LOGIC_WINDOW> cluster_seed = {
        kFixedSeed + 1,
        kFixedSeed + 2,
        kFixedSeed + 3,
        kFixedSeed + 4,
        kFixedSeed + 5};

    std::vector<int> target_num = {10, 10, 10, 10};
    const int target_num_sum = std::accumulate(target_num.begin(), target_num.end(), 0);
    const int cluster_num = 200;

    unsigned int seed = kFixedSeed;
    auto points1 = generate_gaussian_points(target_num[0], 0, 280, 10, 280, 10, 100.0, 50.0, seed++);
    auto points2 = generate_gaussian_points(target_num[1], 0, 280, 10, 150, 10, 100.0, 50.0, seed++);
    auto points3 = generate_gaussian_points(target_num[2], 0, 150, 10, 280, 10, 100.0, 50.0, seed++);
    auto points4 = generate_gaussian_points(target_num[3], 0, 150, 10, 150, 10, 100.0, 50.0, seed++);

    const std::array<double, 7> cluster_param = {
        0.0,
        140.0,
        290.0,
        140.0,
        290.0,
        0.0,
        50.0};

    std::vector<TrackPoint> targets;
    targets.reserve(static_cast<std::size_t>(target_num_sum));
    targets.insert(targets.end(), points1.begin(), points1.end());
    targets.insert(targets.end(), points2.begin(), points2.end());
    targets.insert(targets.end(), points3.begin(), points3.end());
    targets.insert(targets.end(), points4.begin(), points4.end());

    TrackExtrapolator extrapolator(seed);
    for (const auto &point : targets)
    {
        extrapolator.getErrorDistribution(point.x, point.y);
    }

    std::vector<std::vector<TrackPoint>> frames;
    frames.reserve(LOGIC_WINDOW);
    for (std::size_t frame_index = 0; frame_index < LOGIC_WINDOW; ++frame_index)
    {
        if (frame_index > 0)
        {
            extrapolator.update(targets, TIME_INTERVAL_S);
        }

        std::vector<TrackPoint> frame_points = targets;
        auto clutter = generate_uniform_points(cluster_num,
                                               static_cast<int64_t>(frame_index * TIME_INTERVAL_S * 1000.0),
                                               cluster_param[1],
                                               cluster_param[2],
                                               cluster_param[3],
                                               cluster_param[4],
                                               cluster_param[5],
                                               cluster_param[6],
                                               cluster_seed[frame_index]);
        frame_points.insert(frame_points.end(), clutter.begin(), clutter.end());
        frames.push_back(std::move(frame_points));
    }

    LOG_INFO << "========== 4/5 逻辑法起始 ==========";
    LOG_INFO << "固定随机种子: " << kFixedSeed;
    LOG_INFO << "每帧目标数: " << target_num_sum << ", 每帧杂波数: " << cluster_num;
    LOG_INFO << "theta 门限基础参数: " << BASE_THETA_SIGMA_DEG << " deg";
    LOG_INFO << "rho 门限基础参数: " << BASE_R_SIGMA_KM << " km";
    LOG_INFO << "最大速度约束: " << VELOCITY_MAX << " m/s";

    const auto selected_tracks = simple_logic_initiation(frames);
    print_track_summary(selected_tracks);

    const auto &truth_reference = frames.front();
    std::vector<bool> target_matched(static_cast<std::size_t>(target_num_sum), false);
    int true_track_count = 0;
    int false_track_count = 0;
    int missed_target_count = 0;

    for (std::size_t track_index = 0; track_index < selected_tracks.size(); ++track_index)
    {
        const auto &track = selected_tracks[track_index];
        bool valid_track = false;

        for (int target_index = 0; target_index < target_num_sum; ++target_index)
        {
            std::size_t label_match_count = 0;
            for (const TrackPoint *point : track.points)
            {
                if (point != nullptr && same_label(*point, truth_reference[static_cast<std::size_t>(target_index)]))
                {
                    label_match_count += 1;
                }
            }

            if (label_match_count >= LOGIC_CONFIRM_HITS)
            {
                valid_track = true;
                target_matched[static_cast<std::size_t>(target_index)] = true;
                break;
            }
        }

        if (valid_track)
        {
            true_track_count += 1;
        }
        else
        {
            false_track_count += 1;
        }
    }

    for (int target_index = 0; target_index < target_num_sum; ++target_index)
    {
        if (!target_matched[static_cast<std::size_t>(target_index)])
        {
            missed_target_count += 1;
        }
    }

    const int total_targets = target_num_sum;
    const int total_tracks = static_cast<int>(selected_tracks.size());
    const double detection_rate = total_targets > 0 ? (total_targets - missed_target_count) * 100.0 / total_targets : 0.0;
    const double false_alarm_rate = total_tracks > 0 ? false_track_count * 100.0 / total_tracks : 0.0;
    const double missed_rate = total_targets > 0 ? missed_target_count * 100.0 / total_targets : 0.0;

    LOG_INFO << "========== 统计结果 ==========";
    LOG_INFO << "真实目标数: " << total_targets;
    LOG_INFO << "生成航迹数: " << total_tracks;
    LOG_INFO << "正确航迹数: " << true_track_count;
    LOG_INFO << "虚警航迹数: " << false_track_count;
    LOG_INFO << "漏检目标数: " << missed_target_count;
    LOG_INFO << "检测率: " << detection_rate << " %";
    LOG_INFO << "虚警率: " << false_alarm_rate << " %";
    LOG_INFO << "漏检率: " << missed_rate << " %";

    return 0;
}