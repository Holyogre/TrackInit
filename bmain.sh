#!/bin/bash

# 项目构建脚本
# 功能：
# 1. 生成compile_commands.json供代码分析工具使用
# 2. 可配置的编译选项（清理构建、DEBUG模式）
# 3. 可配置的测试输出级别
# 数据库生成命令：
# cd build && cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..

# ===== 用户配置区 =====
CLEAN_BUILD="modeB"   # modeA:彻底删除构建目录; modeB:用make clean; modeC:不清理
DEBUG_MODE=false       # true: Debug模式; false: Release模式
TEST_VERBOSE=true     # true: 输出详细测试报告; false: 仅输出错误
# ======================

# 进入构建目录
mkdir -p build
cd build || exit 1

# === 清理阶段 ===
case "$CLEAN_BUILD" in
    "modeA")
        echo "执行彻底清理：删除整个build目录..."
        cd ..
        rm -rf build
        mkdir build
        cd build || exit 1
        ;;
    "modeB")
        echo "执行增量清理：运行make clean..."
        make clean || { echo "警告：make clean失败，继续构建..."; }
        ;;
    "modeC")
        echo "跳过清理阶段..."
        ;;
    *)
        echo "错误：未知的CLEAN_BUILD模式 '$CLEAN_BUILD'" >&2
        exit 1
        ;;
esac

# === 生成编译命令数据库 ===
echo "生成compile_commands.json..."
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 .. || exit 1

# === WSL环境适配 ===
if grep -q "microsoft" /proc/version 2>/dev/null; then
    echo "检测到WSL环境，应用内存随机化补丁..."
    sudo sysctl -w vm.mmap_rnd_bits=28 vm.mmap_rnd_compat_bits=28 >/dev/null
fi

# === 配置构建类型 ===
BUILD_TYPE=$([ "$DEBUG_MODE" = true ] && echo "Debug" || echo "Release")
echo "配置构建类型: $BUILD_TYPE..."
cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" .. || exit 1

# === 编译项目 ===
echo "开始编译（使用$(nproc)个线程）..."
make -j$(nproc) || exit 1

# === 运行测试 ===
echo "运行测试..."
if [ "$TEST_VERBOSE" = true ]; then
    ctest -V > report.txt 2>&1  # 详细输出到文件
else
    ctest --output-on-failure > report.txt 2>&1  # 仅输出失败信息
fi
echo "测试报告已保存到 build/report.txt"