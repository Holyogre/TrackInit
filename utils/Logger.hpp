/*****************************************************************************
 * @file logger.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 通用日志库
 * 线程安全
 * 启用宏定义时数据将同步显示到日志文件中，CMAKE种进行检查 
 * 退化方式： LOG_DEBUG退化为无输出，LOG_INFO退化为cout，LOG_ERROR退化为cerr
 * @version 0.1
 * @date 2025-11-29
 *
 * @copyright Copyright (c) 2025
 *
 *****************************************************************************/
#pragma once
#include <string>
#include <sstream>
#include <functional>
#include <iostream> 

class Logger
{
public:
    virtual void debug(const std::string &msg) = 0;
    virtual void info(const std::string &msg) = 0;
    virtual void error(const std::string &msg) = 0;
    static Logger *getInstance();
};

#ifdef ENABLE_SPDLOG
class LogStream
{
public:
    explicit LogStream(std::function<void(const std::string &)> sink) : sink_(std::move(sink)) {}
    ~LogStream()
    {
        if (sink_)
            sink_(oss_.str());
    }

    template <typename T>
    LogStream &operator<<(T &&v)
    {
        oss_ << std::forward<T>(v);
        return *this;
    }

    LogStream &operator<<(std::ostream &(*manip)(std::ostream &))
    {
        manip(oss_);
        return *this;
    }

private:
    std::ostringstream oss_;
    std::function<void(const std::string &)> sink_;
};

#define LOG_DEBUG ::LogStream([](const std::string &s) { Logger::getInstance()->debug(s); })
#define LOG_INFO ::LogStream([](const std::string &s) { Logger::getInstance()->info(s); })
#define LOG_ERROR ::LogStream([](const std::string &s) { Logger::getInstance()->error(s); })

#else

#define LOG_DEBUG \
    if (false)    \
    std::cout
#define LOG_INFO std::cout
#define LOG_ERROR std::cerr

#endif