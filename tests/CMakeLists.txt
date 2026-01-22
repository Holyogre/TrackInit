find_package(Catch2 REQUIRED CONFIG HINTS "/usr/local/lib/cmake/Catch2") 


# ==================== 测试函数 ====================
# 定义一个函数来简化测试目标的创建
function(add_track_test TEST_NAME SOURCE_FILE)
    # 捕获额外参数作为编译文件名
    set(IMPL_FILES ${ARGN})

    add_executable(${TEST_NAME}
        ${SOURCE_FILE}
        ${IMPL_FILES}
        ../utils/Logger.cpp  # 添加 Logger 实现
    )

    target_include_directories(${TEST_NAME} PRIVATE 
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/utils  # 添加 utils 包含路径
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/TrackManager
        ${OpenCV_INCLUDE_DIRS}  # 添加 OpenCV 包含目录
    )

    target_link_libraries(${TEST_NAME} PRIVATE 
        Catch2::Catch2WithMain
        spdlog::spdlog  # 链接 spdlog
        "${OpenCV_LIBS}"  # 添加 OpenCV 库，推荐加引号防止变量为空时报错
    )

    # 关键：定义宏
    target_compile_definitions(${TEST_NAME} PRIVATE
        ENABLE_SPDLOG
    )

    # 启用 ASan（仅在 Debug 模式下）
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${TEST_NAME} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
        target_link_options(${TEST_NAME} PRIVATE -fsanitize=address)
    endif()

    set_target_properties(${TEST_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endfunction()


# 使用函数添加测试
# add_track_test(testObjectPool
#                 test_ObjectPool.cpp
#                 ../src/ObjectPool.hpp)
# add_track_test(test_ClustDBSCAN
#                 test_ClustDBSCAN.cpp
#                 ../src/Func_dbscan.cpp)
add_track_test(test_HoughSlice
                test_HoughSlice.cpp
                ../src/HoughSlice.cpp
                ../src/Func_dbscan.cpp
                ${MYUTILS_SOURCES}
                ${TRACKMANAGER_SOURCES})



# 注册数据层滚动
# add_test(NAME testObjectPool功能测试 COMMAND testObjectPool "[memory_pool][]")
# add_test(NAME ClustDBSCAN COMMAND test_ClustDBSCAN "[dbscan][]")
add_test(NAME HoughSlice功能测试 COMMAND test_HoughSlice "[FunctionalityCheck]")
