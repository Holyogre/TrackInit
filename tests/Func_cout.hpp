#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <unordered_map>

// 重载 << 打印 std::pair
std::ostream &operator<<(std::ostream &os, const std::pair<int, int> &p)
{
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

// 重载 << 打印一维向量
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        os << vec[i];
        if (i != vec.size() - 1)
        {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

// 重载 << 运算符用于 unordered_set
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::unordered_set<T> &uset)
{
    os << "{";
    // 使用标志位控制逗号打印（避免末尾多余逗号）
    bool first = true;

    // 遍历 unordered_set
    for (const auto &elem : uset)
    {
        if (!first)
        {
            os << ", ";
        }
        os << elem; // 依赖元素类型自身的 << 重载
        first = false;
    }

    os << "}";
    return os;
}

// 重载 << 操作符用于打印 unordered_map
template<typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& map) {
    os << "{";
    bool first = true;
    for (const auto& pair : map) {
        if (!first) {
            os << ", ";
        }
        os << pair.first << ": " << pair.second;
        first = false;
    }
    os << "}";
    return os;
}

// 字符串转换为二维数组：[]替换为 {}
std::vector<std::vector<int>> str2Vec(const std::string &str)
{
    std::string s = str;
    // 去除所有空格
    s.erase(std::remove_if(s.begin(), s.end(), [](char c)
                           { return std::isspace(static_cast<unsigned char>(c)); }),
            s.end());

    if (s.empty())
    {
        return {};
    }

    std::vector<std::vector<int>> res;
    std::vector<int> cur;
    int num = 0;
    int sign = 1;
    int state = 0; // 0:初始, 1:在二维数组内, 2:在子数组内, 3:解析数字
    int i = 0;
    bool break_out = false;

    while (i < s.size() && !break_out)
    {
        switch (state)
        {
        case 0: // 初始状态
            if (s[i] == '[')
            {
                state = 1;
                i++;
            }
            else
            {
                throw std::invalid_argument("Expected '[' at start");
            }
            break;

        case 1: // 在二维数组内
            if (s[i] == '[')
            {
                cur = std::vector<int>(); // 开始新子数组
                state = 2;
                i++;
            }
            else if (s[i] == ']')
            {
                // 整个数组结束
                i++;
                break_out = true;
            }
            else if (s[i] == ',')
            {
                i++; // 跳过逗号
            }
            else
            {
                throw std::invalid_argument("Unexpected char in state1");
            }
            break;

        case 2: // 在子数组内
            if (s[i] == ']')
            {
                res.push_back(cur);
                state = 1; // 返回二维数组状态
                i++;
            }
            else if (s[i] == ',')
            {
                i++; // 跳过逗号
            }
            else if (s[i] == '-')
            {
                sign = -1;
                i++;
                state = 3; // 解析数字
            }
            else if (std::isdigit(static_cast<unsigned char>(s[i])))
            {
                state = 3; // 解析数字（不移动i）
            }
            else
            {
                throw std::invalid_argument("Unexpected char in state2");
            }
            break;

        case 3: // 解析数字
            num = 0;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
            {
                num = num * 10 + (s[i] - '0');
                i++;
            }
            num *= sign;
            cur.push_back(num);
            sign = 1;  // 重置符号
            state = 2; // 返回子数组状态
            break;
        }
    }
    return res;
}
#endif // UTILS_HPP