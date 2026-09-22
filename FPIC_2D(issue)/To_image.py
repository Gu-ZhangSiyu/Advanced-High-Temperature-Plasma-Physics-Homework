from pathlib import Path
import subprocess
import sys
import time
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import numpy as np

class Config:
    def __init__(self):
        # =====================================================================
        # 模式选择
        # 1: 单个bin文件 -> 粒子密度沉积图
        # 7: 多个bin文件 -> 粒子密度沉积 GIF
        self.function_control = 1
        # =====================================================================
        self.base_dir = r"dynamic_ions_position"    # 输入文件夹：base_dir
        self.single_step = 20           # 模式 1 文件：base_dir / f"{single_step}.bin"
        self.start_step = 10            # 起始文件步
        self.end_step = 20              # 结束文件步
        self.step_interval = 10         # 文件步长
        self.data_column_count = 4      # 粒子文件结构[id, x, y, z]
        self.bins_z_density = 200       # 密度图Z方向网格数
        self.bins_x_density = 200       # 密度图X方向网格数
        self.ticks_x = 5                # 密度图X轴主刻度数量
        self.ticks_z = 5                # 密度图Z轴主刻度数量
        self.x_min = -0.02              # X坐标范围最小值
        self.x_max = 0.02               # X坐标范围最大值
        self.z_min = -0.02              # Z坐标范围最小值
        self.z_max = 0.02               # Z坐标范围最小值
        self.fps = 20                   # 模式 7 GIF 帧数
        self.time_step = 5e-13          # 模拟时间步长 (s)，只用于标题时间
        self.confirm_before_dynamic = True

class PlotConfig:
    @staticmethod
    def apply():
        mpl.rcParams["axes.titlesize"] = 20
        mpl.rcParams["axes.labelsize"] = 22
        mpl.rcParams["xtick.labelsize"] = 18
        mpl.rcParams["ytick.labelsize"] = 18
        mpl.rcParams["figure.titlesize"] = 20
        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.serif"] = ["Times New Roman"] + mpl.rcParams.get("font.serif", [])
        mpl.rcParams["mathtext.fontset"] = "stix"
        mpl.rcParams["mathtext.rm"] = "Times New Roman"
        mpl.rcParams["mathtext.default"] = "it"

def validate_plot_config(config):
    if config.data_column_count != 4:
        raise ValueError("data_column_count 必须为 4，因为 diagnostics_solver.f90 的位置记录是 [id, x, y, z]")
    if config.bins_z_density < 1 or config.bins_x_density < 1:
        raise ValueError("bins_z_density 和 bins_x_density 必须大于 0。")
    if config.x_max <= config.x_min:
        raise ValueError("x_max 必须大于 x_min。")
    if config.z_max <= config.z_min:
        raise ValueError("z_max 必须大于 z_min。")

def get_coordinate_limits(config):
    x_limits = (config.x_min, config.x_max)
    z_limits = (config.z_min, config.z_max)
    return x_limits, z_limits

def read_xz_positions(filepath, components=4):
    filepath = Path(filepath)
    if not filepath.is_file():
        raise FileNotFoundError(f"输入文件不存在: {filepath}")
    if components != 4:
        raise ValueError("位置 bin 文件必须按 4 列 [id, x, y, z] 读取")
    file_size = filepath.stat().st_size
    if file_size == 0:
        empty = np.empty(0, dtype=np.float64)
        return None, empty, empty
    record_size = components * np.dtype(np.float64).itemsize
    if file_size % record_size != 0:
        raise ValueError(f"{filepath} 的字节数 {file_size} 不是单条记录 {record_size} 字节的整数倍")
    raw = np.memmap(filepath, dtype=np.float64, mode="r")
    table = raw.reshape(-1, components)
    return raw, table[:, 1], table[:, 3]


def expected_steps(config):
    if config.start_step < 0 or config.end_step < 0:
        raise ValueError("start_step 和 end_step 不能为负数")
    if config.step_interval <= 0:
        raise ValueError("step_interval 必须大于 0")
    if config.end_step < config.start_step:
        raise ValueError("end_step 必须大于或等于 start_step")
    steps = list(range(config.start_step, config.end_step + 1, config.step_interval))
    if not steps:
        raise ValueError("当前时间步参数没有产生任何输入帧")
    return steps

def validate_input_files(input_dir, steps):
    input_dir = Path(input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"文件夹不存在: {input_dir}")
    missing_files = [input_dir / f"{step}.bin" for step in steps if not (input_dir / f"{step}.bin").is_file()]
    if missing_files:
        preview = "\n".join(str(path) for path in missing_files[:20])
        suffix = "" if len(missing_files) <= 20 else f"\n...另有 {len(missing_files) - 20} 个缺失文件"
        raise FileNotFoundError(f"缺少输入文件：\n{preview}{suffix}")

def make_density(config, x_values, z_values):
    x_limits, z_limits = get_coordinate_limits(config)
    histogram, _, _ = np.histogram2d(z_values, x_values, bins=(config.bins_z_density, config.bins_x_density), range=(z_limits, x_limits),)
    return histogram.T

def make_xz_figure(config):
    x_limits, z_limits = get_coordinate_limits(config)
    z_edges = np.linspace(z_limits[0], z_limits[1], config.bins_z_density + 1)
    x_edges = np.linspace(x_limits[0], x_limits[1], config.bins_x_density + 1)
    initial_density = np.zeros((config.bins_x_density, config.bins_z_density), dtype=np.float64)
    cmap = mpl.colormaps.get_cmap("jet").copy()
    mesh_figure, ax = plt.subplots(figsize=(14, 8))
    mesh = ax.pcolormesh(z_edges, x_edges, initial_density, cmap=cmap, vmin=0.0, vmax=1.0, shading="flat",)
    colorbar = mesh_figure.colorbar(mesh, ax=ax, label="Macro-particle count", shrink=0.82, aspect=30,)
    ax.set_xlabel("Z (m)")
    ax.set_ylabel("X (m)")
    ax.set_xlim(z_limits)
    ax.set_ylim(x_limits)
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(LinearLocator(numticks=config.ticks_z))
    ax.yaxis.set_major_locator(LinearLocator(numticks=config.ticks_x))
    ax.grid(True, alpha=0.25)
    title = ax.set_title("X-Z particle density deposition")
    mesh_figure.tight_layout()
    return {"figure": mesh_figure, "ax": ax, "mesh": mesh, "colorbar": colorbar, "title": title, "x_edges": x_edges, "z_edges": z_edges,}

def update_xz_figure(plot, density, title_text):
    mesh = plot["mesh"]
    mesh.set_array(np.asarray(density, dtype=np.float64).ravel())
    finite_max = np.nanmax(density) if np.any(np.isfinite(density)) else 0.0
    mesh.set_clim(0.0, max(float(finite_max), 1.0))
    plot["colorbar"].update_normal(mesh)
    plot["title"].set_text(title_text)
    plot["figure"].canvas.draw()

def connect_static_click_callback(plot, density):
    def on_click(event):
        if event.inaxes != plot["ax"] or event.xdata is None or event.ydata is None:
            return
        z_index = np.searchsorted(plot["z_edges"], event.xdata, side="right") - 1
        x_index = np.searchsorted(plot["x_edges"], event.ydata, side="right") - 1
        z_index = int(np.clip(z_index, 0, density.shape[1] - 1))
        x_index = int(np.clip(x_index, 0, density.shape[0] - 1))
        count = density[x_index, z_index]
        print("-" * 40)
        print("选定 X-Z 空间网格:")
        print(f"Z 坐标 = {event.xdata:.6g} m")
        print(f"X 坐标 = {event.ydata:.6g} m")
        print(f"局部宏粒子数 = {int(count)}")
    plot["figure"].canvas.mpl_connect("button_press_event", on_click)

def create_gif_writer(config):
    if animation.writers.is_available("ffmpeg"):
        return animation.FFMpegWriter(fps=config.fps, codec="gif")
    return animation.PillowWriter(fps=config.fps)

def show_progress(message, completed, total, width=30):
    ratio = min(max(completed / total, 0.0), 1.0)
    filled = int(width * ratio)
    encoding = (getattr(sys.stdout, "encoding", "") or "").lower().replace("-", "")
    if encoding in {"utf8", "utf16", "utf32"}:
        bar = "█" * filled + "░" * (width - filled)
    else:
        bar = "#" * filled + "-" * (width - filled)
    print(f"\r{message} [{bar}] {ratio * 100:6.2f}%", end="", flush=True)

def read_density_for_step(input_dir, step, config):
    filepath = Path(input_dir) / f"{step}.bin"
    raw, x_values, z_values = read_xz_positions(filepath, components=config.data_column_count)
    try:
        particle_count = int(x_values.size)
        density = make_density(config, x_values, z_values)
        deposited_count = int(np.nansum(density))
    finally:
        del x_values, z_values, raw
    return density, particle_count, deposited_count

def run_static_mode(config):
    print("******************************")
    print(">>> 模式 1: 单张 X-Z 粒子密度沉积图")
    validate_plot_config(config)
    input_dir = Path(config.base_dir)
    input_file = input_dir / f"{config.single_step}.bin"
    validate_input_files(input_dir, [config.single_step])
    print(f"输入文件: {input_file}")
    start_time = time.perf_counter()
    density, particle_count, deposited_count = read_density_for_step(input_dir, config.single_step, config)
    plot = make_xz_figure(config)
    current_time_us = config.single_step * config.time_step * 1e6
    update_xz_figure(plot, density, f"time={current_time_us:.6f} μs, step={config.single_step}, " f"particles={deposited_count}",)
    connect_static_click_callback(plot, density)
    print(f"文件粒子数: {particle_count}")
    print(f"落入当前 X-Z 坐标范围并完成沉积的宏粒子数: {deposited_count}")
    print(f"读取和绘图耗时: {time.perf_counter() - start_time:.2f} 秒")
    print("图形已打开：如需保存，请在图形窗口中手动保存，建议 DPI 使用 200。")
    print("交互模式已开启：点击图中网格可在终端查看局部宏粒子数。")
    plt.show()
    plt.close(plot["figure"])

def write_dynamic_gif(config, input_dir, steps, output_path, plot):
    def write_with(writer):
        first_frame_count = None
        last_frame_count = None
        with writer.saving(plot["figure"], str(output_path), dpi=200):
            for frame_number, step in enumerate(steps, start=1):
                density, particle_count, deposited_count = read_density_for_step(input_dir, step, config)
                current_time_us = step * config.time_step * 1e6
                update_xz_figure(plot, density, f"time={current_time_us:.6f} μs, step={step}, " f"particles={deposited_count}",)
                writer.grab_frame(facecolor=plot["figure"].get_facecolor())
                if frame_number == 1:
                    first_frame_count = (particle_count, deposited_count)
                if frame_number == len(steps):
                    last_frame_count = (particle_count, deposited_count)
                show_progress(f"正在读取 {frame_number}/{len(steps)} 个文件", frame_number, len(steps),)
        print()
        return first_frame_count, last_frame_count
    writer = create_gif_writer(config)
    try:
        return write_with(writer)
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        if not isinstance(writer, animation.FFMpegWriter):
            raise
        print(f"ffmpeg GIF 写入失败，切换 PillowWriter：{error}")
        return write_with(animation.PillowWriter(fps=config.fps))

def run_dynamic_mode(config):
    print("******************************")
    print(">>> 模式 7: 生成一个 X-Z 粒子密度沉积 GIF")
    validate_plot_config(config)
    if config.fps <= 0:
        raise ValueError("fps 必须大于 0。")
    if config.time_step <= 0:
        raise ValueError("time_step 必须大于 0。")
    steps = expected_steps(config)
    input_dir = Path(config.base_dir)
    validate_input_files(input_dir, steps)
    output_path = Path(__file__).resolve().parent / "electrons_density_xz.gif"
    print(f"输入目录: {input_dir}")
    print(f"时间步: {steps[0]} -> {steps[-1]}，间隔 {config.step_interval}")
    print(f"帧数: {len(steps)}，GIF 帧率: {config.fps}")
    print(f"物理时间: {steps[0] * config.time_step * 1e6:.6f} -> " f"{steps[-1] * config.time_step * 1e6:.6f} μs，" f"time_step={config.time_step:.6e} s")
    print(f"输出文件: {output_path}")
    try:
        confirmation = input("确认以上输入范围并开始生成 GIF？输入 y 后开始：").strip().lower()
    except EOFError:
        print("未检测到交互输入，已取消模式 7。")
        return
    if confirmation != "y":
        print("未输入 y，已取消模式 7。")
        return
    plot = make_xz_figure(config)
    total_start = time.perf_counter()
    try:
        first_count, last_count = write_dynamic_gif(config, input_dir, steps, output_path, plot)
    finally:
        plt.close(plot["figure"])
    print(f"已生成: {output_path}")
    if first_count is not None:
        print(f"首帧文件粒子数/沉积数: {first_count[0]}/{first_count[1]}；" f"末帧文件粒子数/沉积数: {last_count[0]}/{last_count[1]}")
    print(f"GIF 生成完成，总耗时: {time.perf_counter() - total_start:.2f} 秒")

def main():
    config = Config()
    PlotConfig.apply()
    if config.function_control == 1:
        run_static_mode(config)
    elif config.function_control == 7:
        run_dynamic_mode(config)
    else:
        raise ValueError("function_control 只能设置为 1（单张图）或 7（动态图）")

if __name__ == "__main__":
    main()

