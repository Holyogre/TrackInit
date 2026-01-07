#include <catch2/catch_all.hpp>
#include "../src/ObjectPool.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <memory>
#include <algorithm>
#include <cstring>

using namespace track_project::trackinit;

// 写个友元用来REQUIRE判断私有成员是否正确得了
//  在测试文件中定义与头文件中声明匹配的友元模板
namespace track_project::trackinit
{
    template <typename U>
    struct ObjectPoolTestFriend
    {
        static size_t allocated_count(const ObjectPool<U> &p)
        {
            size_t cnt = 0;
            for (auto v : p.alloc_flags_)
                if (v)
                    ++cnt;
            return cnt;
        }

        static size_t ptr_map_size(const ObjectPool<U> &p) { return p.ptr_to_idx_.size(); }

        static size_t free_slots_size(const ObjectPool<U> &p) { return p.free_slots_.size(); }
    };
} // namespace track_project::trackinit

constexpr size_t TEST_OBJ_SIZE_BYTES = 4;
constexpr size_t TEST_OBJPOOL_SIZE = 1024;

struct TestObj
{
    int v[TEST_OBJ_SIZE_BYTES / sizeof(int)]; // 占128字节
    void clear() { std::memset(v, 0, sizeof(v)); }
};

TEST_CASE("测试ObjectPool的基本功能", "[memory_pool][basic]")
{
    ObjectPool<TestObj> pool(4);

    SECTION("扩容后指针有效性检查")
    {
        // 初始 4 个槽位，先占满再触发扩容
        std::vector<TestObj *> first_batch;
        for (int i = 0; i < 4; ++i)
        {
            auto p = pool.acquire_target();
            p->v[0] = i + 1; // 写入标记，后面验证未被破坏
            first_batch.push_back(p);
        }

        auto before_expand_ptr0 = first_batch[0];

        // 第五次获取会触发 expand_pool，旧指针应保持有效
        auto expanded_ptr = pool.acquire_target();
        REQUIRE(expanded_ptr != nullptr);
        expanded_ptr->v[0] = 42;

        // 验证扩容后旧指针仍可读写且数据未损坏
        REQUIRE(before_expand_ptr0->v[0] == 1);
        before_expand_ptr0->v[0] = 7;
        REQUIRE(before_expand_ptr0->v[0] == 7);

        // 释放并检查内部计数状态
        for (auto p : first_batch)
        {
            REQUIRE(pool.release_target(p));
        }
        REQUIRE(pool.release_target(expanded_ptr));

        REQUIRE(ObjectPoolTestFriend<TestObj>::allocated_count(pool) == 0);
        REQUIRE(ObjectPoolTestFriend<TestObj>::ptr_map_size(pool) == 0);
        // 扩容到 8 个槽位（4->8），全部应回到空闲栈
        REQUIRE(ObjectPoolTestFriend<TestObj>::free_slots_size(pool) == 8);
    }

    SECTION("清空指令有效性检查")
    {
        std::vector<TestObj *> ptrs;
        for (int i = 0; i < 3; ++i)
        {
            auto p = pool.acquire_target();
            p->v[0] = 100 + i; // 写入非零数据
            ptrs.push_back(p);
        }

        pool.clear_all();

        // clear_all 后不应有已分配指针
        REQUIRE(ObjectPoolTestFriend<TestObj>::allocated_count(pool) == 0);
        REQUIRE(ObjectPoolTestFriend<TestObj>::ptr_map_size(pool) == 0);
        REQUIRE(ObjectPoolTestFriend<TestObj>::free_slots_size(pool) == 4);

        // 再次获取，应该复用旧指针且数据已被 clear
        std::vector<TestObj *> reacquired;
        for (int i = 0; i < 3; ++i)
        {
            auto p = pool.acquire_target();
            reacquired.push_back(p);
            REQUIRE(p->v[0] == 0);
        }

        // 至少应有一个指针与之前相同，证明复用而非新分配
        bool reused = false;
        for (auto old_p : ptrs)
        {
            for (auto new_p : reacquired)
            {
                if (old_p == new_p)
                {
                    reused = true;
                }
            }
        }
        REQUIRE(reused);
    }
}

TEST_CASE("ObjectPool - 完整基准测试", "[memory_pool][benchmark]")
{
    const size_t N = 10000;

    BENCHMARK("容量为1的对象池扩容至10k测试") // 一次出现1000个群
    {
        ObjectPool<TestObj> pool(1);
        std::vector<TestObj *> v;
        v.reserve(N);
        for (size_t i = 0; i < N; ++i)
        {
            v.push_back(pool.acquire_target());
        }
        for (auto p : v)
            pool.release_target(p);
        return v.size();
    };

    BENCHMARK("10k对象的池反复获取释放速度测试")
    {
        ObjectPool<TestObj> pool(N);
        for (size_t i = 0; i < N; ++i)
        {
            auto p = pool.acquire_target();
            pool.release_target(p);
        }
        return N;
    };

    BENCHMARK("10k对象的池获取指针速度测试")
    {
        ObjectPool<TestObj> pool(TEST_OBJPOOL_SIZE);
        for (size_t i = 0; i < N; ++i)
        {
            auto p = pool.acquire_target();
        }
        auto it = pool.get_allocated_ptrs();
        return it.size();
    };
}
