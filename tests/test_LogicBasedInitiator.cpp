#include <catch2/catch_all.hpp>
#include "../src/LogicBasedInitiator.hpp"
#include "../include/ManagementService.hpp"
#include <../include/defsystem.h>
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <memory>
#include <algorithm>
#include <cstring>
#include <set>
#include <thread>
#include <condition_variable>
#include <csignal>
#include <chrono>
#include <numeric>
#include <iomanip>

#include "../utils/Logger.hpp"

using namespace track_project;
using track_project::trackinit::LogicBasedInitiator;

// 为了匹配 LogicBasedInitiator.hpp 中的 "friend class test_LogicBasedInitiator;"
// 将测试辅助类放到同一个命名空间 track_project::trackinit 下
namespace track_project::trackinit
{
    static constexpr size_t MAX_BINS = LOGIC_BASED_NUM_X_BINS * LOGIC_BASED_NUM_Y_BINS; // 最大允许的距离门数量乘以角度门数量

    // test_LogicBasedInitiator.h
    class test_LogicBasedInitiator
    {
    public:
        explicit test_LogicBasedInitiator(LogicBasedInitiator &init) : initiator_(init) {}

        // 只保留这一个函数：保存二进制格式
        bool saveErrorDistributionToDat(const std::string &filename) const
        {
            std::ofstream file(filename, std::ios::binary);
            if (!file.is_open())
            {
                return false;
            }

            // 直接写入所有sigma_x, sigma_y交替
            for (const auto &item : initiator_.error_distribution_table_)
            {
                double sigma_x = item.first;
                double sigma_y = item.second;
                file.write(reinterpret_cast<const char *>(&sigma_x), sizeof(double));
                file.write(reinterpret_cast<const char *>(&sigma_y), sizeof(double));
            }

            file.close();
            return true;
        }

    private:
        LogicBasedInitiator &initiator_;
    };

} // namespace track_project::trackinit

// 方便在测试代码中直接使用不带命名空间前缀的名称
using track_project::trackinit::test_LogicBasedInitiator;

TEST_CASE("基础功能正确性测试", "[FunctionalityCheck]")
{
    LogicBasedInitiator initiator;
    test_LogicBasedInitiator tester(initiator);

    // 误差分布表格保存
    REQUIRE(tester.saveErrorDistributionToDat("../error_distribution.dat") == true);
}
