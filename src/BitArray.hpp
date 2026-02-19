/*****************************************************************************
 * @file BitArray.hpp
 * @author xjl (xjl20011009@126.com)
 * @brief 模板化的BitArray，位数在编译期确定
 * @version 0.3
 * @date 2026-02-19
 *
 * @copyright Copyright (c) 2026
 *
 *****************************************************************************/
#include <cstdint>
#include <array>
#include <algorithm>
#include <stdexcept>

namespace track_project::trackinit
{
    template <size_t NumBits>
    class BitArray
    {
    private:
        static constexpr size_t BITS_PER_WORD = 32;
        static constexpr size_t NUM_WORDS = (NumBits + BITS_PER_WORD - 1) / BITS_PER_WORD;
        static constexpr size_t NUM_BYTES = NumBits / 8;

        // 编译期断言：位数必须是32的倍数
        static_assert(NumBits % BITS_PER_WORD == 0, "Bits must be multiple of 32");

        std::array<uint32_t, NUM_WORDS> data{};

        // 确保位索引有效
        void check_bit_index(size_t pos) const
        {
            if (pos >= NumBits)
            {
                throw std::out_of_range("Bit position out of range");
            }
        }

        // 确保字节索引有效
        void check_byte_index(size_t pos) const
        {
            if (pos >= NUM_BYTES)
            {
                throw std::out_of_range("Byte position out of range");
            }
        }

    public:
        // 构造函数 - 默认全零
        BitArray() = default;

        // 拷贝构造
        BitArray(const BitArray &) = default;

        // 赋值运算符
        BitArray &operator=(const BitArray &) = default;

        // 位与操作
        BitArray &operator&=(const BitArray &other)
        {
            for (size_t i = 0; i < NUM_WORDS; ++i)
            {
                data[i] &= other.data[i];
            }
            return *this;
        }

        // 位或操作
        BitArray &operator|=(const BitArray &other)
        {
            for (size_t i = 0; i < NUM_WORDS; ++i)
            {
                data[i] |= other.data[i];
            }
            return *this;
        }
        // 对指定区域进行位或操作（按字节）
        void or_bytes(size_t byte_offset, const uint8_t *mask_bytes, size_t num_bytes)
        {
            std::vector<uint8_t> buffer(num_bytes);
            read_bytes(buffer.data(), byte_offset, num_bytes);
            for (size_t i = 0; i < num_bytes; i++)
            {
                buffer[i] |= mask_bytes[i];
            }
            write_bytes(buffer.data(), byte_offset, num_bytes);
        }


        // 逻辑左移
        BitArray &operator<<=(size_t shift)
        {
            if (shift == 0)
                return *this;
            if (shift >= NumBits)
            {
                data.fill(0);
                return *this;
            }

            size_t wordShift = shift / BITS_PER_WORD;
            size_t bitShift = shift % BITS_PER_WORD;

            std::array<uint32_t, NUM_WORDS> result{};

            if (bitShift == 0)
            {
                for (size_t i = wordShift; i < NUM_WORDS; ++i)
                {
                    result[i] = data[i - wordShift];
                }
            }
            else
            {
                for (size_t i = 0; i < NUM_WORDS - wordShift; ++i)
                {
                    result[i + wordShift] |= (data[i] << bitShift);
                    if (i + wordShift + 1 < NUM_WORDS)
                    {
                        result[i + wordShift + 1] |= (data[i] >> (BITS_PER_WORD - bitShift));
                    }
                }
            }

            data = std::move(result);
            return *this;
        }

        // 逻辑右移
        BitArray &operator>>=(size_t shift)
        {
            if (shift == 0)
                return *this;
            if (shift >= NumBits)
            {
                data.fill(0);
                return *this;
            }

            size_t wordShift = shift / BITS_PER_WORD;
            size_t bitShift = shift % BITS_PER_WORD;

            std::array<uint32_t, NUM_WORDS> result{};

            if (bitShift == 0)
            {
                for (size_t i = 0; i < NUM_WORDS - wordShift; ++i)
                {
                    result[i] = data[i + wordShift];
                }
            }
            else
            {
                for (size_t i = 0; i < NUM_WORDS - wordShift; ++i)
                {
                    result[i] |= (data[i + wordShift] >> bitShift);
                    if (i + wordShift + 1 < NUM_WORDS)
                    {
                        result[i] |= (data[i + wordShift + 1] << (BITS_PER_WORD - bitShift));
                    }
                }
            }

            data = std::move(result);
            return *this;
        }

        // 置零操作
        void clear()
        {
            data.fill(0);
        }

        // 写入单个字节到指定位置
        void write_byte(uint8_t value, size_t pos)
        {
            check_byte_index(pos);

            size_t bitPos = pos * 8;
            size_t wordIndex = bitPos / BITS_PER_WORD;
            size_t bitOffset = bitPos % BITS_PER_WORD;

            // 清除目标位置的8位
            uint32_t mask = ~(static_cast<uint32_t>(0xFF) << bitOffset);
            data[wordIndex] &= mask;

            // 写入新值
            data[wordIndex] |= (static_cast<uint32_t>(value) << bitOffset);
        }

        // 读取单个字节从指定位置
        uint8_t read_byte(size_t pos) const
        {
            check_byte_index(pos);

            size_t bitPos = pos * 8;
            size_t wordIndex = bitPos / BITS_PER_WORD;
            size_t bitOffset = bitPos % BITS_PER_WORD;

            return static_cast<uint8_t>((data[wordIndex] >> bitOffset) & 0xFF);
        }

        // 批量写入多个字节
        void write_bytes(const uint8_t *buffer, size_t startPos, size_t length)
        {
            if (startPos + length > NUM_BYTES)
            {
                throw std::out_of_range("Write operation would exceed bit array size");
            }

            for (size_t i = 0; i < length; ++i)
            {
                write_byte(buffer[i], startPos + i);
            }
        }

        // 批量读取多个字节
        void read_bytes(uint8_t *buffer, size_t startPos, size_t length) const
        {
            if (startPos + length > NUM_BYTES)
            {
                throw std::out_of_range("Read operation would exceed bit array size");
            }

            for (size_t i = 0; i < length; ++i)
            {
                buffer[i] = read_byte(startPos + i);
            }
        }

        // 获取指定位的值
        bool get_bit(size_t pos) const
        {
            check_bit_index(pos);
            size_t wordIndex = pos / BITS_PER_WORD;
            size_t bitIndex = pos % BITS_PER_WORD;
            return (data[wordIndex] >> bitIndex) & 1;
        }

        // 设置指定位的值
        void set_bit(size_t pos, bool value)
        {
            check_bit_index(pos);
            size_t wordIndex = pos / BITS_PER_WORD;
            size_t bitIndex = pos % BITS_PER_WORD;
            if (value)
            {
                data[wordIndex] |= (1U << bitIndex);
            }
            else
            {
                data[wordIndex] &= ~(1U << bitIndex);
            }
        }

        // 编译期常量：获取位数
        static constexpr size_t size() { return NumBits; }

        // 编译期常量：获取字节数
        static constexpr size_t byte_size() { return NUM_BYTES; }

        // 获取数据（用于调试）
        const std::array<uint32_t, NUM_WORDS> &get_data() const { return data; }
    };
}