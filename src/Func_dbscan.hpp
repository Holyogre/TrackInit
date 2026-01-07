/*****************************************************************************
MIT License

Original work Copyright (c) 2021 Eleobert do Espírito Santo
Modified work Copyright (c) 2025 xjl (xjl20011009@126.com)

Original source: https://github.com/eleobert/dbscan
Modifications include:
- 重命名参数和常量（以符合当前项目规范）
- 修改为仅支持TrackPoint数据结构，此外float均修改成了double型

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*****************************************************************************/

#ifndef _CLUST_DBSCAN_HPP_
#define _CLUST_DBSCAN_HPP_

#include <cassert>
#include <cstddef>
#include <vector>
#include <cstdlib>
#include "defstruct.h"

namespace track_project::trackinit
{
    /*****************************************************************************
     * @brief DBSCAN聚类算法实现（基于Eleobert do Espírito Santo的原始版本）
     *
     * @param data 待聚类的TrackPoint点集
     * @param eps 邻域搜索半径
     * @param min_pts 核心点所需的最小邻域点数
     * @return std::vector<std::vector<size_t>> 聚类结果，每个子向量包含属于同一聚类的点索引
     *
     * @note 原始实现：Eleobert do Espírito Santo (2021)
     *       修改：xjl (2025-12-18)
     * @version 1.0 (修改版)
     * @author xjl (xjl20011009@126.com)
     * @date 2025-12-18
     *****************************************************************************/
    std::vector<std::vector<size_t>> dbscan(const std::vector<TrackPoint> &data, double eps, int min_pts);
}
#endif
