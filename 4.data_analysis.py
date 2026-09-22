import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter, minimum_filter
from scipy.integrate import cumulative_trapezoid

# 路径配置
DIR_BX = r"E:\NIFS_Upload\高温等离子体作业\result_b3\dynamic_fields_Bx"
DIR_BY = r"E:\NIFS_Upload\高温等离子体作业\result_b3\dynamic_fields_By"

START_STEP = 100        # 文件头
END_STEP = 20000         # 文件尾
STEP_INTERVAL = 100     # 文件步长
dt = 0.5e-12            # 时间步长 (s)
dx = 0.5e-3             # X方向网格间距 (m)
dy = 0.5e-3             # 假设Y方向网格间距与X相同 (m)

def compute_vector_potential(Bx, By):
    ny, nx = Bx.shape
    Az = np.zeros((ny, nx))

    # 沿底边 y=0 积分
    Az[0, :] = cumulative_trapezoid(-By[0, :], dx=dx, initial=0)

    # 沿 y 方向积分, 在 axis=0 维度上进行累积积分
    Az = Az[0, :] + cumulative_trapezoid(Bx, dx=dy, axis=0, initial=0)
    return Az

def compute_reconnection_flux_2d(Bx, By, Az):
    B_mag = np.sqrt(Bx ** 2 + By ** 2)

    # 寻找局部极小值区域
    local_min_mask = (minimum_filter(B_mag, size=5) == B_mag)

    # 设置阈值过滤本底噪声
    threshold = np.max(B_mag) * 0.1
    null_mask = local_min_mask & (B_mag < threshold)

    az_at_nulls = Az[null_mask]

    if len(az_at_nulls) >= 2:
        # 当捕捉到分离的 X 点和 O 点时
        delta_psi = np.max(az_at_nulls) - np.min(az_at_nulls)
    else:
        # 演化极早期, 拓扑尚未完全分离, 退化为在中心电流片区域附近寻找 Az 的极值差
        ny = Az.shape[0]
        mid = ny // 2

        # 取中心附近 20% 区域
        sheet_slice = Az[mid - ny // 10: mid + ny // 10, :]
        delta_psi = np.max(sheet_slice) - np.min(sheet_slice)
    return delta_psi

def plot_magnetic_reconnection(step, dir_bx, dir_by):
    file_bx = os.path.join(dir_bx, f"{step}.txt")
    file_by = os.path.join(dir_by, f"{step}.txt")
    try:
        Bx = np.loadtxt(file_bx)
        By = np.loadtxt(file_by)
        Bx = gaussian_filter(Bx, sigma=1.5)
        By = gaussian_filter(By, sigma=1.5)
    except FileNotFoundError:
        print(f"Error: 找不到第 {step} 步的磁场拓扑数据文件！")
        return

    ny, nx = Bx.shape
    x = np.linspace(-nx / 2, nx / 2, nx) * dx
    y = np.linspace(-ny / 2, ny / 2, ny) * dy
    X, Y = np.meshgrid(x, y)

    B_mag = np.sqrt(Bx ** 2 + By ** 2)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)

    # 背景着色 Bx
    im_bg = ax.pcolormesh(X, Y, Bx, cmap='seismic', shading='auto', alpha=0.5)

    # 流线图：提升 density 增加磁力线密度
    strm = ax.streamplot(X, Y, Bx, By, color=B_mag, cmap='viridis', linewidth=1.0, density=3.0)

    cbar_bg = fig.colorbar(im_bg, ax=ax, fraction=0.04, pad=0.1)
    cbar_bg.set_label('$B_x$ [T]', fontsize=12)
    cbar_line = fig.colorbar(strm.lines, ax=ax, fraction=0.04, pad=0.05)
    cbar_line.set_label('|B| [T]', fontsize=12)

    ax.set_xlabel('X (m)', fontsize=12)
    ax.set_ylabel('Y (m)', fontsize=12)
    plt.tight_layout()
    plt.show()

def compute_and_plot_reconnection_rate(start_step, end_step, step_interval, dir_bx, dir_by):
    valid_steps = []
    delta_psi_values = []

    for step in range(start_step, end_step + 1, step_interval):
        file_bx = os.path.join(dir_bx, f"{step}.txt")
        file_by = os.path.join(dir_by, f"{step}.txt")
        try:
            Bx = np.loadtxt(file_bx)
            By = np.loadtxt(file_by)

            # 平滑处理，降低网格噪声对磁零点搜索的干扰
            Bx = gaussian_filter(Bx, sigma=1.0)
            By = gaussian_filter(By, sigma=1.0)

            # 计算全场 Az
            Az = compute_vector_potential(Bx, By)

            # 全局搜索极值提取真实磁通差
            delta_psi = compute_reconnection_flux_2d(Bx, By, Az)

            valid_steps.append(step)
            delta_psi_values.append(delta_psi)
        except FileNotFoundError:
            continue

    if not valid_steps:
        print("Error: 未找到数据文件，请检查演化路径。")
        return

    time_points = np.array(valid_steps) * dt * 1e9  # 转换为 ns
    delta_psi_array = np.array(delta_psi_values)
    time_seconds = time_points * 1e-9

    # 使用二阶中心差分计算重联率
    reconnection_rate = np.gradient(delta_psi_array, time_seconds)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), dpi=150, sharex=True)
    ax1.plot(time_points, delta_psi_array, marker='.', linestyle='-', color='navy', linewidth=1.5)
    ax1.set_ylabel(r"$\Delta \psi$ (Wb/m)", fontsize=12)
    ax1.grid(True, linestyle=':', alpha=0.7)
    ax2.plot(time_points, reconnection_rate, marker='.', linestyle='-', color='firebrick', linewidth=1.2)
    ax2.set_xlabel("Time (ns)", fontsize=12)
    ax2.set_ylabel(r"$\partial(\Delta\psi)/\partial t$ (V/m)", fontsize=12)
    ax2.axhline(0, color='black', linestyle='--', linewidth=1.2)
    ax2.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_magnetic_reconnection(END_STEP, DIR_BX, DIR_BY)
    compute_and_plot_reconnection_rate(START_STEP, END_STEP, STEP_INTERVAL, DIR_BX, DIR_BY)