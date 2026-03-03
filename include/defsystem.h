namespace track_project
{
    // 雷达站基础信息
    constexpr double base_longitude = 0.0;
    constexpr double base_latitude = 0.0;
    constexpr double base_normal_direction = 0.0; // 法线方向：北偏东，单位°，范围仅允许填写[0,360)，不然结果错误
    constexpr double base_theta_sigma_deg = 3.0;  // 先验航向误差标准差，单位°，
    constexpr double base_r_sigma_km = 1.0;       // 先验距离误差标准差，单位KM，

    // 目标基础信息
    constexpr double velocity_max = 250.0; // 理论目标最大速度，单位m/s
}