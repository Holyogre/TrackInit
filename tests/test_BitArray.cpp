#include <catch2/catch_all.hpp>
#include "../src/BitArray.hpp"
#include <cstring>
#include <random>

using namespace track_project::trackinit;

// 定义一些常用的BitArray类型别名
using BitArray32 = BitArray<32>;
using BitArray64 = BitArray<64>;
using BitArray96 = BitArray<96>;
using BitArray128 = BitArray<128>;
using BitArray256 = BitArray<256>;
using BitArray512 = BitArray<512>;
using BitArray1024 = BitArray<1024>;
using BitArray4096 = BitArray<4096>;

TEST_CASE("BitArray的功能测试", "[bitarray][basic]"){

    SECTION("构造函数和基本属性测试") {
        BitArray128 bits;
        REQUIRE(bits.size() == 128);
        REQUIRE(bits.get_data().size() == 4);
        REQUIRE(bits.byte_size() == 16);
        
        // 模板参数必须是32的倍数，否则编译失败（由static_assert保证）
        // 所以不需要运行时异常测试了
    }

    SECTION("字节写入和读取测试") {
        BitArray128 bits;
        
        // 单个字节写入读取
        bits.write_byte(0x12, 0);
        bits.write_byte(0x34, 1);
        bits.write_byte(0x56, 2);
        bits.write_byte(0x78, 3);
        
        REQUIRE(bits.read_byte(0) == 0x12);
        REQUIRE(bits.read_byte(1) == 0x34);
        REQUIRE(bits.read_byte(2) == 0x56);
        REQUIRE(bits.read_byte(3) == 0x78);
        
        // 覆盖写入
        bits.write_byte(0xFF, 1);
        REQUIRE(bits.read_byte(1) == 0xFF);
        
        // 批量写入读取
        uint8_t write_buf[8] = {0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88};
        uint8_t read_buf[8] = {0};
        
        bits.write_bytes(write_buf, 8, 8);
        bits.read_bytes(read_buf, 8, 8);
        
        for(int i = 0; i < 8; i++) {
            REQUIRE(read_buf[i] == write_buf[i]);
        }
    }

    SECTION("位设置和获取测试") {
        BitArray64 bits;
        
        bits.set_bit(0, true);
        bits.set_bit(31, true);
        bits.set_bit(32, true);
        bits.set_bit(63, true);
        
        REQUIRE(bits.get_bit(0) == true);
        REQUIRE(bits.get_bit(31) == true);
        REQUIRE(bits.get_bit(32) == true);
        REQUIRE(bits.get_bit(63) == true);
        
        // 通过字节验证
        REQUIRE(bits.read_byte(0) == 0x01);        // 第0位为1
        REQUIRE(bits.read_byte(3) == 0x80);        // 第31位为1
        REQUIRE(bits.read_byte(4) == 0x01);        // 第32位为1
        REQUIRE(bits.read_byte(7) == 0x80);        // 第63位为1
    }

    SECTION("位与操作测试") {
        BitArray64 bits1;
        BitArray64 bits2;
        
        // 设置测试数据：交替模式
        for(size_t i = 0; i < 8; i++) {
            if(i % 2 == 0) {
                bits1.write_byte(0xFF, i);  // 11111111
                bits2.write_byte(0x0F, i);  // 00001111
            } else {
                bits1.write_byte(0x00, i);
                bits2.write_byte(0xF0, i);  // 11110000
            }
        }
        
        bits1 &= bits2;
        
        // 验证结果
        REQUIRE(bits1.read_byte(0) == 0x0F);  // 0xFF & 0x0F = 0x0F
        REQUIRE(bits1.read_byte(1) == 0x00);  // 0x00 & 0xF0 = 0x00
        REQUIRE(bits1.read_byte(2) == 0x0F);
        REQUIRE(bits1.read_byte(3) == 0x00);
        
        // 不同大小的BitArray不能直接操作（编译期类型不同）
        // 所以不需要运行时异常测试
    }

    SECTION("位或操作测试") {
        BitArray64 bits1;
        BitArray64 bits2;
        
        // 设置测试数据
        bits1.write_byte(0x0F, 0);  // 00001111
        bits1.write_byte(0xF0, 2);  // 11110000
        bits2.write_byte(0xF0, 0);  // 11110000
        bits2.write_byte(0x0F, 2);  // 00001111
        
        bits1 |= bits2;
        
        // 验证结果
        REQUIRE(bits1.read_byte(0) == 0xFF);  // 0x0F | 0xF0 = 0xFF
        REQUIRE(bits1.read_byte(1) == 0x00);
        REQUIRE(bits1.read_byte(2) == 0xFF);  // 0xF0 | 0x0F = 0xFF
        REQUIRE(bits1.read_byte(3) == 0x00);
    }

    SECTION("左移操作测试") {
        BitArray128 bits;
        
        // 设置测试模式
        bits.write_byte(0x12, 0);
        bits.write_byte(0x34, 1);
        bits.write_byte(0x56, 2);
        bits.write_byte(0x78, 3);
        
        // 左移8位
        bits <<= 8;
        REQUIRE(bits.read_byte(0) == 0x00);
        REQUIRE(bits.read_byte(1) == 0x12);
        REQUIRE(bits.read_byte(2) == 0x34);
        REQUIRE(bits.read_byte(3) == 0x56);
        REQUIRE(bits.read_byte(4) == 0x78);
        
        // 左移16位
        BitArray128 bits2;
        bits2.write_byte(0x12, 0);
        bits2.write_byte(0x34, 1);
        bits2 <<= 16;
        REQUIRE(bits2.read_byte(2) == 0x12);
        REQUIRE(bits2.read_byte(3) == 0x34);
        
        // 左移超过位数
        BitArray64 bits3;
        bits3.write_byte(0xFF, 0);
        bits3 <<= 64;
        REQUIRE(bits3.read_byte(0) == 0x00);
    }

    SECTION("右移操作测试") {
        BitArray128 bits;
        
        // 设置测试模式（从中间开始）
        bits.write_byte(0x12, 4);
        bits.write_byte(0x34, 5);
        bits.write_byte(0x56, 6);
        bits.write_byte(0x78, 7);
        
        // 右移8位
        bits >>= 8;
        REQUIRE(bits.read_byte(3) == 0x12);
        REQUIRE(bits.read_byte(4) == 0x34);
        REQUIRE(bits.read_byte(5) == 0x56);
        REQUIRE(bits.read_byte(6) == 0x78);
        REQUIRE(bits.read_byte(7) == 0x00);
        
        // 右移16位
        BitArray128 bits2;
        bits2.write_byte(0x12, 4);
        bits2.write_byte(0x34, 5);
        bits2 >>= 16;
        REQUIRE(bits2.read_byte(2) == 0x12);
        REQUIRE(bits2.read_byte(3) == 0x34);
    }

    SECTION("清零操作测试") {
        BitArray64 bits;
        
        // 写入所有字节
        for(size_t i = 0; i < 8; i++) {
            bits.write_byte(0xFF, i);
        }
        
        bits.clear();
        
        // 验证所有字节清零
        for(size_t i = 0; i < 8; i++) {
            REQUIRE(bits.read_byte(i) == 0x00);
        }
        
        // 验证所有位清零
        for(size_t i = 0; i < 64; i++) {
            REQUIRE(bits.get_bit(i) == false);
        }
    }

    SECTION("拷贝构造和赋值测试") {
        BitArray96 bits1;  // 12字节
        
        // 写入测试数据
        for(size_t i = 0; i < 12; i++) {
            bits1.write_byte(static_cast<uint8_t>(i), i);
        }
        
        // 拷贝构造
        BitArray96 bits2(bits1);
        REQUIRE(bits2.size() == bits1.size());
        for(size_t i = 0; i < 12; i++) {
            REQUIRE(bits2.read_byte(i) == bits1.read_byte(i));
        }
        
        // 赋值操作
        BitArray96 bits3;
        bits3 = bits1;
        REQUIRE(bits3.size() == bits1.size());
        for(size_t i = 0; i < 12; i++) {
            REQUIRE(bits3.read_byte(i) == bits1.read_byte(i));
        }
    }
}

TEST_CASE("BitArray的边界测试", "[bitarray][edge]"){
    
    SECTION("最小尺寸测试") {
        BitArray32 bits;  // 最小合法尺寸
        REQUIRE(bits.size() == 32);
        REQUIRE(bits.byte_size() == 4);
        
        // 测试边界字节访问
        bits.write_byte(0x12, 0);
        bits.write_byte(0x34, 3);  // 最后一个字节
        
        REQUIRE(bits.read_byte(0) == 0x12);
        REQUIRE(bits.read_byte(3) == 0x34);
        
        REQUIRE_THROWS_AS(bits.write_byte(0x56, 4), std::out_of_range);
    }
    
    SECTION("最大尺寸测试（编译期大数组）") {
        constexpr size_t LARGE_BITS = 1024 * 32 * 8;  // 256K位 = 32KB
        BitArray<LARGE_BITS> bits;
        
        // 在边界位置写入
        bits.write_byte(0xAA, 0);
        bits.write_byte(0xBB, bits.byte_size() - 1);
        
        REQUIRE(bits.read_byte(0) == 0xAA);
        REQUIRE(bits.read_byte(bits.byte_size() - 1) == 0xBB);
    }
    
    SECTION("边界位访问测试") {
        BitArray64 bits;
        
        // 测试每个字节的边界位
        for(size_t byte = 0; byte < 8; byte++) {
            size_t first_bit = byte * 8;
            size_t last_bit = byte * 8 + 7;
            
            bits.set_bit(first_bit, true);
            bits.set_bit(last_bit, true);
            
            REQUIRE(bits.get_bit(first_bit) == true);
            REQUIRE(bits.get_bit(last_bit) == true);
            
            // 验证字节值
            uint8_t expected = 0x81;  // 10000001
            REQUIRE(bits.read_byte(byte) == expected);
            
            // 清除以便下次测试
            bits.write_byte(0x00, byte);
        }
    }
    
    SECTION("跨字边界移位测试") {
        BitArray96 bits;  // 3个32位字
        
        // 设置一个跨越字边界的模式
        // 让第31位和第32位（跨字边界）都为1
        bits.set_bit(31, true);
        bits.set_bit(32, true);
        
        // 左移1位
        bits <<= 1;
        
        // 验证：原来的第31位移到第32位，原来的第32位移到第33位
        REQUIRE(bits.get_bit(32) == true);
        REQUIRE(bits.get_bit(33) == true);
        
        // 右移1位回来
        bits >>= 1;
        REQUIRE(bits.get_bit(31) == true);
        REQUIRE(bits.get_bit(32) == true);
    }
    
    SECTION("全零和全一边界测试") {
        BitArray96 bits;
        
        // 全一测试
        for(size_t i = 0; i < 12; i++) {
            bits.write_byte(0xFF, i);
        }
        
        // 验证全一
        for(size_t i = 0; i < 96; i++) {
            REQUIRE(bits.get_bit(i) == true);
        }
        
        // 与操作后应为全一
        BitArray96 bits2;
        for(size_t i = 0; i < 12; i++) {
            bits2.write_byte(0xFF, i);
        }
        bits &= bits2;
        
        for(size_t i = 0; i < 96; i++) {
            REQUIRE(bits.get_bit(i) == true);
        }
        
        // 清零
        bits.clear();
        
        // 全零测试
        for(size_t i = 0; i < 96; i++) {
            REQUIRE(bits.get_bit(i) == false);
        }
    }
    
    SECTION("异常测试") {
        BitArray64 bits;
        
        // 越界位访问
        REQUIRE_THROWS_AS(bits.get_bit(64), std::out_of_range);
        REQUIRE_THROWS_AS(bits.set_bit(64, true), std::out_of_range);
        
        // 越界字节访问
        REQUIRE_THROWS_AS(bits.read_byte(8), std::out_of_range);
        REQUIRE_THROWS_AS(bits.write_byte(0xFF, 8), std::out_of_range);
        
        // 批量操作越界
        uint8_t buf[4] = {0};
        REQUIRE_THROWS_AS(bits.read_bytes(buf, 6, 4), std::out_of_range);
        REQUIRE_THROWS_AS(bits.write_bytes(buf, 6, 4), std::out_of_range);
    }
}

TEST_CASE("BitArray的性能测试", "[bitarray][benchmark]"){
    
    SECTION("批量字节操作性能") {
        constexpr size_t BYTES = 1024 * 128;  // 128KB
        constexpr size_t BITS = BYTES * 8;
        BitArray<BITS> bits;
        uint8_t* buffer = new uint8_t[BYTES];
        
        // 准备测试数据
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 255);
        
        for(size_t i = 0; i < BYTES; i++) {
            buffer[i] = static_cast<uint8_t>(dis(gen));
        }
        
        BENCHMARK("批量写入 " + std::to_string(BYTES) + " 字节") {
            bits.write_bytes(buffer, 0, BYTES);
            return bits.size();
        };
        
        BENCHMARK("批量读取 " + std::to_string(BYTES) + " 字节") {
            uint8_t* read_buf = new uint8_t[BYTES];
            bits.read_bytes(read_buf, 0, BYTES);
            delete[] read_buf;
            return bits.size();
        };
        
        delete[] buffer;
    }
    
    SECTION("位操作性能") {
        constexpr size_t BYTES = 1024 * 64;  // 64KB
        constexpr size_t BITS = BYTES * 8;
        BitArray<BITS> bits;
        
        // 填充随机数据
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> dis;
        
        // 现在get_data()返回std::array的引用
        auto& data = const_cast<std::array<uint32_t, BITS/32>&>(bits.get_data());
        for(auto& word : data) {
            word = dis(gen);
        }
        
        BENCHMARK("设置所有位") {
            for(size_t i = 0; i < BITS; i++) {
                bits.set_bit(i, true);
            }
            return bits.size();
        };
        
        BENCHMARK("读取所有位") {
            volatile bool result;
            for(size_t i = 0; i < BITS; i++) {
                result = bits.get_bit(i);
            }
            return result;
        };
    }
    
    SECTION("移位性能测试") {
        constexpr size_t BYTES = 1024 * 64;  // 64KB
        constexpr size_t BITS = BYTES * 8;
        BitArray<BITS> bits;
        
        // 填充随机数据
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> dis;
        
        auto& data = const_cast<std::array<uint32_t, BITS/32>&>(bits.get_data());
        for(auto& word : data) {
            word = dis(gen);
        }
        
        BENCHMARK("左移1位") {
            bits <<= 1;
            return bits.size();
        };
        
        BENCHMARK("左移32位") {
            bits <<= 32;
            return bits.size();
        };
        
        BENCHMARK("右移1位") {
            bits >>= 1;
            return bits.size();
        };
        
        BENCHMARK("右移32位") {
            bits >>= 32;
            return bits.size();
        };
    }
    
    SECTION("逻辑运算性能") {
        constexpr size_t BYTES = 1024 * 64;  // 64KB
        constexpr size_t BITS = BYTES * 8;
        BitArray<BITS> bits1;
        BitArray<BITS> bits2;
        
        // 填充随机数据
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint32_t> dis;
        
        auto& data1 = const_cast<std::array<uint32_t, BITS/32>&>(bits1.get_data());
        auto& data2 = const_cast<std::array<uint32_t, BITS/32>&>(bits2.get_data());
        
        for(auto& word : data1) word = dis(gen);
        for(auto& word : data2) word = dis(gen);
        
        BENCHMARK("位与操作") {
            bits1 &= bits2;
            return bits1.size();
        };
        
        BENCHMARK("位或操作") {
            bits1 |= bits2;
            return bits1.size();
        };
    }
    
    SECTION("拷贝性能") {
        constexpr size_t BYTES = 1024 * 64;  // 64KB
        constexpr size_t BITS = BYTES * 8;
        BitArray<BITS> bits;
        
        BENCHMARK("拷贝构造") {
            BitArray<BITS> copy(bits);
            return copy.size();
        };
        
        BENCHMARK("赋值操作") {
            BitArray<BITS> copy;
            copy = bits;
            return copy.size();
        };
    }
}