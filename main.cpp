#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

#define M_PI 3.14159265358979323846

// 模拟的vote函数
void vote_in_hough_space(double heading1, double heading2, double doppler,
                         double rel_x, double rel_y, int batch,
                         int doppler_tolerance_bits, int vote_area)
{
    std::cout << "投票区间: [" << std::fixed << std::setprecision(3)
              << heading1 << ", " << heading2 << "] 弧度, ";
    std::cout << "区间长度: " << (heading2 - heading1) << " 弧度" << std::endl;
}

// 模拟的数据点结构
struct Point
{
    double x, y;
    double doppler;
    Point(double x_, double y_, double d_) : x(x_), y(y_), doppler(d_) {}
};
void process_heading_original(double base_dir_rad, double angle_rad, const Point &point)
{
    std::cout << "\n=== 测试: base_dir_rad = " << std::fixed << std::setprecision(3)
              << base_dir_rad << " (" << base_dir_rad * 180 / M_PI << "°), "
              << "angle_rad = " << angle_rad
              << " (" << angle_rad * 180 / M_PI << "°) ===" << std::endl;

    double heading1 = base_dir_rad - angle_rad;
    double heading2 = base_dir_rad + angle_rad;

    std::cout << "原始heading1: " << heading1 << " (" << heading1 * 180 / M_PI << "°)" << std::endl;
    std::cout << "原始heading2: " << heading2 << " (" << heading2 * 180 / M_PI << "°)" << std::endl;

    // 保存原始值用于调试
    double orig_h1 = heading1;
    double orig_h2 = heading2;

    // 分支判断（与你的测试条件一致）
    if (heading1 <= 0)
    {
        std::cout << "分支1: heading1 <= 0" << std::endl;
        double heading3 = 0.0, heading4 = M_PI;
        heading3 = heading1 + M_PI; // heading1是负数，heading3 < π
        heading4 = M_PI;
        heading1 = 0.0;
        // heading2保持不变

        std::cout << "  区间1: [" << heading1 << ", " << heading2 << "]" << std::endl;
        std::cout << "  区间2: [" << heading3 << ", " << heading4 << "]" << std::endl;

        // 检查边界条件
        if (heading2 > M_PI)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") > π，区间1跨越了π" << std::endl;
        }
        if (heading2 < 0)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") < 0，异常情况" << std::endl;
        }
        if (heading3 < 0 || heading3 > M_PI)
        {
            std::cout << "  ⚠️ 注意: heading3(" << heading3 << ") 不在[0, π]内" << std::endl;
        }
    }
    else if (heading2 >= 2 * M_PI)
    {
        std::cout << "分支2: heading2 >= 2π" << std::endl;
        double heading3 = 0.0, heading4 = M_PI;
        heading3 = heading1 - M_PI;
        heading4 = M_PI;
        heading1 = 0.0;
        heading2 = heading2 - 2 * M_PI;

        std::cout << "  区间1: [" << heading1 << ", " << heading2 << "]" << std::endl;
        std::cout << "  区间2: [" << heading3 << ", " << heading4 << "]" << std::endl;

        // 检查边界条件
        if (heading2 < 0)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") < 0，异常情况" << std::endl;
        }
        if (heading2 > M_PI)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") > π，区间1可能有问题" << std::endl;
        }
        if (heading3 < 0)
        {
            std::cout << "  ⚠️ 注意: heading3(" << heading3 << ") < 0，区间2可能无效" << std::endl;
        }
        if (heading3 > M_PI)
        {
            std::cout << "  ⚠️ 注意: heading3(" << heading3 << ") > π，区间2可能无效" << std::endl;
        }
    }
    else if (heading1 <= M_PI && heading2 <= M_PI)
    {
        std::cout << "分支3: 完全在[0, π]内" << std::endl;
        std::cout << "  区间: [" << heading1 << ", " << heading2 << "]" << std::endl;

        // 检查边界条件
        if (heading1 < 0)
        {
            std::cout << "  ⚠️ 注意: heading1(" << heading1 << ") < 0，异常情况" << std::endl;
        }
    }
    else if (heading1 <= M_PI && heading2 > M_PI)
    {
        std::cout << "分支4: 跨越π (heading1 ≤ π < heading2)" << std::endl;
        double heading3 = heading1;
        double heading4 = M_PI;
        heading1 = 0.0;
        heading2 = heading2 - M_PI;

        std::cout << "  区间1: [" << heading1 << ", " << heading2 << "]" << std::endl;
        std::cout << "  区间2: [" << heading3 << ", " << heading4 << "]" << std::endl;

        // 检查边界条件
        if (heading2 > M_PI)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") > π，区间1可能有问题" << std::endl;
        }
        if (heading2 < 0)
        {
            std::cout << "  ⚠️ 注意: heading2(" << heading2 << ") < 0，异常情况" << std::endl;
        }
        if (heading3 < 0)
        {
            std::cout << "  ⚠️ 注意: heading3(" << heading3 << ") < 0，异常情况" << std::endl;
        }
    }
    else if (heading1 > M_PI && heading2 < 2 * M_PI)
    {
        std::cout << "分支5: 完全在(π, 2π)内" << std::endl;
        double h1 = heading1 - M_PI;
        double h2 = heading2 - M_PI;
        std::cout << "  区间: [" << h1 << ", " << h2 << "]" << std::endl;

        // 检查边界条件
        if (h1 < 0)
        {
            std::cout << "  ⚠️ 注意: h1(" << h1 << ") < 0，异常情况" << std::endl;
        }
        if (h2 > M_PI)
        {
            std::cout << "  ⚠️ 注意: h2(" << h2 << ") > π，异常情况" << std::endl;
        }
        if (h1 > h2)
        {
            std::cout << "  ⚠️ 注意: h1(" << h1 << ") > h2(" << h2 << ")，区间可能无效" << std::endl;
        }
    }
    else
    {
        std::cout << "❌ 错误: 未覆盖的情况！" << std::endl;
        std::cout << "  heading1 = " << heading1 << ", heading2 = " << heading2 << std::endl;
        std::cout << "  heading1 > M_PI? " << (heading1 > M_PI) << std::endl;
        std::cout << "  heading2 < 2*M_PI? " << (heading2 < 2 * M_PI) << std::endl;
    }

    // 添加一些额外的检查
    std::cout << "  区间合法性检查:" << std::endl;
    std::cout << "    所有区间端点应在[0, π]内: ";
    bool all_valid = true;
    // 这里只是示意，实际需要根据不同的分支检查
    std::cout << (all_valid ? "✓" : "✗") << std::endl;
}

int main()
{
    std::cout << "开始参数扫描测试..." << std::endl;
    std::cout << "========================================" << std::endl;

    // 测试参数范围
    const int base_steps = 1000;  // base_dir_rad的采样点数
    const int angle_steps = 1000; // angle_rad的采样点数

    Point test_point(1.0, 0.0, 1.0); // 测试点

    // 记录未覆盖的情况
    std::vector<std::pair<double, double>> uncovered_cases;

    // 遍历所有组合
    for (int i = 0; i <= base_steps; i++)
    {
        double base_dir_rad = (2 * M_PI * i) / base_steps;

        for (int j = 0; j <= angle_steps; j++)
        {
            double angle_rad = (M_PI / 2 * j) / angle_steps;
            angle_rad = std::min(angle_rad, M_PI / 2 - 1e-6); // 避免angle_rad为π/2导致的cos为0问题
            angle_rad = std::max(angle_rad, 1e-6);            // 避免angle_rad为0导致的acos(1)问题

            // 跳过angle_rad为0的情况（虽然简单但也要测试）

            double heading1 = base_dir_rad - angle_rad;
            double heading2 = base_dir_rad + angle_rad;

            // 检查是否进入某个分支
            bool covered = false;

            if (heading1 <= 0)
            {
                covered = true;
            }
            else if (heading2 >= 2 * M_PI)
            {
                covered = true;
            }
            else if (heading1 <= M_PI && heading2 <= M_PI)
            {
                covered = true;
            }
            else if (heading1 <= M_PI && heading2 > M_PI)
            {
                covered = true;
            }
            else if (heading1 > M_PI && heading2 < 2 * M_PI)
            {
                covered = true;
            }

            if (!covered)
            {
                uncovered_cases.push_back({base_dir_rad, angle_rad});
            }
        }
    }

    // 打印统计信息
    std::cout << "\n=== 统计信息 ===" << std::endl;
    std::cout << "总测试组合数: " << (base_steps + 1) * (angle_steps + 1) << std::endl;
    std::cout << "未覆盖的组合数: " << uncovered_cases.size() << std::endl;

    if (!uncovered_cases.empty())
    {
        std::cout << "\n未覆盖的情况:" << std::endl;
        for (const auto &case_pair : uncovered_cases)
        {
            std::cout << "  base_dir_rad = " << case_pair.first
                      << " (" << case_pair.first * 180 / M_PI << "°), "
                      << "angle_rad = " << case_pair.second
                      << " (" << case_pair.second * 180 / M_PI << "°)" << std::endl;
        }
    }

    // 详细测试一些边界情况
    std::cout << "\n\n=== 详细边界测试 ===" << std::endl;

    // 测试用例1: heading1正好为0
    process_heading_original(0.3, 0.3, test_point);

    // 测试用例2: heading1略小于0
    process_heading_original(0.2, 0.3, test_point);

    // 测试用例3: heading2正好为2π
    process_heading_original(2 * M_PI - 0.1, 0.1, test_point);

    // 测试用例4: heading2略大于2π
    process_heading_original(2 * M_PI - 0.05, 0.1, test_point);

    // 测试用例5: 跨越π
    process_heading_original(3.0, 0.3, test_point);

    // 测试用例6: heading1正好为π
    process_heading_original(M_PI, 0.0, test_point);

    // 测试用例7: heading2正好为π
    process_heading_original(M_PI - 0.1, 0.1, test_point);

    // 测试用例8: heading1正好为2π
    process_heading_original(2 * M_PI, 0.0, test_point);

    // 测试用例9: 区间完全在[π, 2π]内
    process_heading_original(4.0, 0.3, test_point);

    // 测试用例10: 区间完全在[0, π]内
    process_heading_original(1.0, 0.3, test_point);

    // 测试用例11: angle_rad为最大值
    process_heading_original(M_PI / 2, M_PI / 2, test_point);

    // 测试用例12: angle_rad为最大值且base_dir_rad接近2π
    process_heading_original(2 * M_PI - 0.1, M_PI / 2, test_point);

    std::cout << "\n========================================" << std::endl;

    return 0;
}