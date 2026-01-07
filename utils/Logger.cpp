// src/logger.cpp
#include "./Logger.hpp"

#include <iostream>
#include <chrono>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <filesystem>

class SpdlogLogger : public Logger
{
public:
    SpdlogLogger()
    {
        try
        {
            // 日志目录从cmakelists种读取 TRACKMANAGER_DEFAULT_LOG_DIR
            std::string log_dir = std::string(TRACKMANAGER_DEFAULT_LOG_DIR);
            ;

            // 生成带日期的文件名
            auto now = std::chrono::system_clock::now();
            std::time_t now_time = std::chrono::system_clock::to_time_t(now);
            std::tm *time_info = std::localtime(&now_time);
            char buffer[80];
            std::strftime(buffer, sizeof(buffer), "%Y-%m-%d", time_info);
            std::string date_str(buffer);

            // 构建日志路径
            std::string log_path = log_dir + "/kalman_" + date_str + ".log";

            // 2. 创建控制台+文件双输出
            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_path, true);

            std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};
            auto logger = std::make_shared<spdlog::logger>("kalman", sinks.begin(), sinks.end());

            // 3. 设为全局默认日志器
            spdlog::set_default_logger(logger);
            spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");
#ifdef NDEBUG
            spdlog::set_level(spdlog::level::info); // Release模式
#else
            spdlog::set_level(spdlog::level::debug); // Debug模式
#endif
        }
        catch (const spdlog::spdlog_ex &ex)
        {
            std::cout << "注意日志文件没有成功输出，检查下文件路径对不对" << std::endl;
            spdlog::set_default_logger(spdlog::stdout_color_mt("kalman_fallback"));
        }
    }

    void debug(const std::string &msg) override { spdlog::debug(msg); }
    void info(const std::string &msg) override { spdlog::info(msg); }
    void error(const std::string &msg) override { spdlog::error(msg); }
};

// 工厂方法实现
Logger *Logger::getInstance()
{
    static SpdlogLogger instance; // 线程安全的单例
    return &instance;
}