namespace track_project
{
    // 雷达站经纬度
    constexpr double base_longitude = 0.0;
    constexpr double base_latitude = 0.0;
    constexpr double base_normal_direction = 0.0; // 法线方向：北偏东，单位°，范围仅允许填写[0,360)，不然结果错误

    // 理论目标最大速度，单位m/s
    constexpr double velocity_max = 300.0;
}