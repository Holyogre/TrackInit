import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_hough_debug_file(filename):
    """读取霍夫空间调试文件"""
    with open(filename, 'rb') as f:
        # 读取维度
        dims = np.fromfile(f, dtype=np.uint32, count=2)
        theta_dim, rho_dim = dims
        
        # 读取霍夫空间数据
        hough_space = np.fromfile(f, dtype=np.uint64, count=theta_dim * rho_dim)
        hough_space = hough_space.reshape(theta_dim, rho_dim)
    
    return hough_space, theta_dim, rho_dim

def visualize_hough_space(hough_space, theta_dim, rho_dim, theta_resolution=1.0, rho_resolution=0.1):
    """可视化霍夫空间"""
    fig = plt.figure(figsize=(15, 10))
    
    # 1. 2D热力图
    ax1 = fig.add_subplot(221)
    im = ax1.imshow(hough_space.T, aspect='auto', origin='lower', cmap='hot', 
                   extent=[0, theta_dim, 0, rho_dim])
    plt.colorbar(im, ax=ax1, label='投票数')
    ax1.set_xlabel(f'θ 索引 (0-{theta_dim-1})')
    ax1.set_ylabel(f'ρ 索引 (0-{rho_dim-1})')
    ax1.set_title('2D霍夫空间热力图')
    
    # 添加实际值的刻度
    ax1_secx = ax1.secondary_xaxis('top', functions=(lambda x: x*theta_resolution, lambda x: x/theta_resolution))
    ax1_secx.set_xlabel(f'角度 (°)')
    
    ax1_secy = ax1.secondary_yaxis('right', functions=(
        lambda y: (y - rho_dim/2)*rho_resolution, 
        lambda y: y/rho_resolution + rho_dim/2
    ))
    ax1_secy.set_ylabel(f'截距 (km)')
    
    # 2. 3D曲面图
    ax2 = fig.add_subplot(222, projection='3d')
    X, Y = np.meshgrid(range(theta_dim), range(rho_dim))
    surf = ax2.plot_surface(X, Y, hough_space.T, cmap='viridis', alpha=0.8)
    ax2.set_xlabel('θ 索引')
    ax2.set_ylabel('ρ 索引')
    ax2.set_zlabel('投票数')
    ax2.set_title('3D霍夫空间曲面')
    
    # 3. θ方向统计
    ax3 = fig.add_subplot(223)
    theta_sum = np.sum(hough_space, axis=1)
    ax3.bar(range(theta_dim), theta_sum, alpha=0.7)
    ax3.set_xlabel('θ 索引')
    ax3.set_ylabel('总投票数')
    ax3.set_title('各θ角的投票统计')
    ax3.grid(True, alpha=0.3)
    
    # 标记前3个峰值
    if len(theta_sum) > 0:
        top3_idx = np.argsort(theta_sum)[-3:][::-1]
        for i, idx in enumerate(top3_idx):
            ax3.plot(idx, theta_sum[idx], 'ro', markersize=8)
            ax3.text(idx, theta_sum[idx], f' θ={idx}\n({idx*theta_resolution}°)', 
                    ha='center', va='bottom', fontsize=9)
    
    # 4. ρ方向统计
    ax4 = fig.add_subplot(224)
    rho_sum = np.sum(hough_space, axis=0)
    ax4.bar(range(rho_dim), rho_sum, alpha=0.7)
    ax4.set_xlabel('ρ 索引')
    ax4.set_ylabel('总投票数')
    ax4.set_title('各ρ截距的投票统计')
    ax4.grid(True, alpha=0.3)
    
    # 标记前3个峰值
    if len(rho_sum) > 0:
        top3_idx = np.argsort(rho_sum)[-3:][::-1]
        for i, idx in enumerate(top3_idx):
            actual_rho = (idx - rho_dim/2) * rho_resolution
            ax4.plot(idx, rho_sum[idx], 'ro', markersize=8)
            ax4.text(idx, rho_sum[idx], f' ρ={idx}\n({actual_rho:.2f}km)', 
                    ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    # 打印统计信息
    total_votes = np.sum(hough_space)
    non_zero = np.count_nonzero(hough_space)
    max_votes = np.max(hough_space)
    max_idx = np.unravel_index(np.argmax(hough_space), hough_space.shape)
    
    print(f"=== 霍夫空间统计信息 ===")
    print(f"维度: {theta_dim} × {rho_dim}")
    print(f"总投票数: {total_votes}")
    print(f"非零单元格: {non_zero} ({non_zero/(theta_dim*rho_dim)*100:.1f}%)")
    print(f"最大投票数: {max_votes} 在 (θ={max_idx[0]}, ρ={max_idx[1]})")
    print(f"峰值位置对应: 角度={max_idx[0]*theta_resolution}°, 截距={(max_idx[1]-rho_dim/2)*rho_resolution:.3f}km")
    
    return fig

def print_nonzero_cells(hough_space, theta_dim, rho_dim, theta_resolution=1.0, rho_resolution=0.1, top_n=20):
    """打印非零单元格信息"""
    nonzero_indices = np.where(hough_space > 0)
    nonzero_values = hough_space[nonzero_indices]
    
    if len(nonzero_values) == 0:
        print("没有非零单元格")
        return
    
    # 按投票数排序
    sorted_indices = np.argsort(nonzero_values)[::-1]
    
    print(f"\n=== 非零单元格详情 (共{len(nonzero_values)}个) ===")
    print("排名 | θ索引 | ρ索引 | 投票数 | 角度(°) | 截距(km)")
    print("-" * 60)
    
    for i, idx in enumerate(sorted_indices[:top_n]):
        theta_idx = nonzero_indices[0][idx]
        rho_idx = nonzero_indices[1][idx]
        votes = nonzero_values[idx]
        angle = theta_idx * theta_resolution
        rho_value = (rho_idx - rho_dim/2) * rho_resolution
        
        print(f"{i+1:4d} | {theta_idx:5d} | {rho_idx:5d} | {votes:7d} | {angle:7.1f} | {rho_value:9.3f}")

# 主程序
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("使用方法: python visualize_hough.py <hough_debug_file.dat> [theta_resolution] [rho_resolution]")
        print("示例: python visualize_hough.py hough_debug_center_53.2_24.1.dat 1.0 0.1")
        sys.exit(1)
    
    filename = sys.argv[1]
    theta_resolution = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    rho_resolution = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    
    try:
        # 读取数据
        hough_space, theta_dim, rho_dim = read_hough_debug_file(filename)
        
        # 打印非零单元格信息
        print_nonzero_cells(hough_space, theta_dim, rho_dim, theta_resolution, rho_resolution)
        
        # 可视化
        fig = visualize_hough_space(hough_space, theta_dim, rho_dim, theta_resolution, rho_resolution)
        
        # 保存图像
        output_png = filename.replace('.dat', '_visualization.png')
        plt.savefig(output_png, dpi=150, bbox_inches='tight')
        print(f"\n可视化结果已保存为: {output_png}")
        
        plt.show()
        
    except Exception as e:
        print(f"错误: {e}")
        print(f"文件格式: 前8字节为2个uint32维度，后接theta_dim*rho_dim个uint64数据")