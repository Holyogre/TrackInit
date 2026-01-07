/*****************************************************************************
 * @file config.hpp
 * @author xjl (xjl20011009@126.com)
 *
 * @brief 配置文件读取工具类
 * 不具备线程安全性，严禁多线程同时操作同一实例
 * 具备事务性，加载失败时回滚旧配置
 *
 * @version 0.3
 * @date 2025-11-07
 *
 * @copyright Copyright (c) 2025
 *****************************************************************************/
#ifndef _TRACK_CONFIG_HPP_
#define _TRACK_CONFIG_HPP_

// 基础结构
#include <string>
#include <cstdint> //精确位数整数
// 文件控制
#include <fstream>
// 算法函数
#include <algorithm>
// 字符处理
#include <cctype>
#include <vector>
// 网络通信
#include <netinet/in.h>
#include <arpa/inet.h>
// 异常处理
#include <assert.h>
// 日志库
#include "../utils/Logger.hpp"

namespace track_project
{
    // 直接读取型配置项数量
    constexpr int direct_read_count = 4;

    class TrackConfig
    {
    public:
        // 直接读取项:⚠️修改此处的时候需要同步修改 applyKeyValue方法以及 direct_read_count 常量
        // std::string track_dst_ip = "127.0.0.1";
        // std::uint16_t trackmanager_dst_port = 5555;
        std::uint16_t trackmanager_recv_port = 5556;
        std::vector<std::string> trackmanager_recv_filters{};

        // 预解析项:
        sockaddr_in trackmanager_dst_sockaddr{}; // 预解析网络结构

    public:
        /*****************************************************************************
         * @brief 构造函数，加载配置文件
         * @param filepath 配置文件路径
         *****************************************************************************/
        explicit TrackConfig(const std::string &filepath)
        {
            bool success = reload(filepath);
            if (!success)
            {
                LOG_ERROR << "ATTENTION!!!CONFIG文件初始化失败";
                assert(0 && "配置文件加载失败，强制终止程序"); // 配置文件加载失败，强制终止程序
            }
        }

        /*****************************************************************************
         * @brief 重新加载配置文件
         * @param filepath 配置文件路径
         * @return true 加载成功
         * @return false 加载失败，使用旧配置
         *****************************************************************************/
        bool reload(const std::string &filepath)
        {
            std::ifstream file(filepath);
            if (!file.is_open())
            {
                LOG_ERROR << "无法打开配置文件: " << filepath << "，使用默认配置";
                return false;
            }

            // 旧配置用于回滚
            auto old = *this;

            int passed_items = 0;
            std::string line;

            while (std::getline(file, line))
            {
                trim(line);
                if (line.empty() || line[0] == '#')
                    continue;

                size_t pos = line.find('=');
                if (pos == std::string::npos)
                    continue;

                std::string key = trim_copy(line.substr(0, pos));
                std::string value = trim_copy(line.substr(pos + 1));

                if (apply_key_value(key, value))
                {
                    passed_items++;
                }
                else
                {
                    *this = old;
                    LOG_ERROR << "配置项应用失败，回滚旧配置";
                    return false;
                }
            }

            if (passed_items < direct_read_count)
            {
                LOG_ERROR << "配置项共<" << passed_items << ">项，少于要求的<" << direct_read_count << ">项，回滚旧配置";
                *this = old;
                return false;
            }

            LOG_INFO << "配置文件重载成功: " << filepath;
            return true;
        }

    private:
        //=== 去除首尾空白 ===
        static void trim(std::string &s)
        {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch)
                                            { return !std::isspace(ch); }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch)
                                 { return !std::isspace(ch); })
                        .base(),
                    s.end());
        }

        static std::string trim_copy(std::string s)
        {
            trim(s);
            return s;
        }

        /******************************************************************************
         * @brief 解析IP地址字符串到 sockaddr_in 结构，在内部处理所有异常
         * @param ipStr 待解析的IP字符串
         * @param outIpStr 输出参数，存储解析后的IP字符串（仅当返回 true 时有效）
         * @param outSockaddr 输出参数，存储解析后的sockaddr_in结构（仅当返回 true 时有效）
         * @return true 解析成功
         * @return false 解析失败
         *****************************************************************************/
        bool parse_sockaddr(
            const std::string &ipStr,
            std::string &outIpStr,   // 用于更新 track_dst_ip
            sockaddr_in &outSockaddr // 用于更新 dst_sockaddr
        )
        {
            sockaddr_in tempAddr = {};
            tempAddr.sin_family = AF_INET;

            if (inet_pton(AF_INET, ipStr.c_str(), &tempAddr.sin_addr) != 1)
            {
                LOG_ERROR << "IP 地址格式非法：" << ipStr;
                return false;
            }

            // 只有解析成功才更新输出参数
            outIpStr = ipStr;
            outSockaddr = tempAddr;
            return true;
        }

        /*****************************************************************************
         * @brief 解析端口号字符串到 uint16_t,在内部处理所有异常
         * @param portStr 待解析的端口字符串
         * @param outPort 输出参数，存储解析后的端口值（仅当返回 true 时有效）
         * @return true 解析成功
         * @return false 解析失败
         *****************************************************************************/
        bool parse_port(const std::string &portStr, std::uint16_t &outPort)
        {
            try
            {
                size_t pos = 0;
                int port_int = std::stoi(portStr, &pos); // pos 会指向解析结束的位置

                // 检查是否整个字符串都被解析（避免 "123.456" 这样的情况）
                if (pos != portStr.length())
                {
                    LOG_ERROR << "端口值无效 [port=" << portStr << "]: 包含非法字符";
                    return false;
                }

                if (port_int < 1 || port_int > 65535)
                {
                    LOG_ERROR << "端口值无效 [port=" << portStr << "]: 必须在 1–65535 范围内";
                    return false;
                }

                outPort = static_cast<std::uint16_t>(port_int);
                return true;
            }
            catch (const std::invalid_argument &)
            {
                LOG_ERROR << "端口值无效 [port=" << portStr << "]: 不是有效数字";
            }
            catch (const std::out_of_range &)
            {
                LOG_ERROR << "端口值无效 [port=" << portStr << "]: 数值超出范围";
            }
            catch (...)
            {
                LOG_ERROR << "端口值无效 [port=" << portStr << "]: 未知解析错误";
            }
            return false;
        }

        /*****************************************************************************
         * @brief 解析逗号分隔的过滤规则到 std::vector<std::string>
         * @param filtersStr 待解析的字符串（格式: "TRACK_, SYSTEM_"）
         * @param outFilters 输出参数，存储解析后的过滤规则
         * @return true 解析成功（至少有一条有效规则）
         * @return false 解析失败（无有效规则）
         *****************************************************************************/
        bool parse_filters(const std::string &filtersStr, std::vector<std::string> &outFilters)
        {
            outFilters.clear();
            size_t start = 0, end = 0;

            while (start < filtersStr.size())
            {
                end = filtersStr.find(',', start);
                if (end == std::string::npos)
                    end = filtersStr.size();

                std::string filter = filtersStr.substr(start, end - start);

                // 去除前后空格
                size_t firstNonSpace = filter.find_first_not_of(" \t\n\r");
                if (firstNonSpace != std::string::npos)
                {
                    size_t lastNonSpace = filter.find_last_not_of(" \t\n\r");
                    filter = filter.substr(firstNonSpace, lastNonSpace - firstNonSpace + 1);
                }

                if (!filter.empty())
                {
                    outFilters.push_back(filter);
                }

                start = end + 1;
            }

            if (outFilters.empty())
            {
                LOG_ERROR << "过滤规则解析失败: 未找到有效规则 [filters=" << filtersStr << "]";
                return false;
            }

            return true;
        }

        /*****************************************************************************
         * @brief 应用单个键值对到配置项
         * @param key 配置项名称
         * @param value 配置项值
         * @return true 应用成功
         * @return false 应用失败（解析错误等）
         *****************************************************************************/
        bool apply_key_value(const std::string &key, const std::string &value)
        {

            LOG_DEBUG << "应用配置项: " << key << " = " << value;

            if (false) // 方便注释
            {
            }
            // if (key == "track_dst_ip")
            // {
            //     return parseSockaddr(value, track_dst_ip, trackmanager_dst_sockaddr);
            // }
            // else if (key == "trackmanager_dst_port")
            // {
            //     return parsePort(value, trackmanager_dst_port);
            // }
            else if (key == "trackmanager_recv_port")
            {
                return parse_port(value, trackmanager_recv_port);
            }
            else if (key == "trackmanager_recv_filters")
            {
                return parse_filters(value, trackmanager_recv_filters);
            }
            else
            {
                LOG_INFO << "未知配置项: " << key << " = " << value;
                return false;
            }
        }
    };

} // namespace track_project

#endif // _TRACK_CONFIG_HPP_
