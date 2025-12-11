import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
import glob

# ================= 配置参数 =================
CONFIG = {
    'results_dir': r"./results_lgd",  # 上一步代码保存结果的文件夹
    'region_box': [88, 92, 22, 26],  # [lon_min, lon_max, lat_min, lat_max] (红色框)
    'map_extent': [70, 110, -80, 80],  # 左图地图显示的范围 [lon_min, lon_max, lat_min, lat_max]
    'lgd_scale': 1e9,  # 将LGD转换为 nm/s^2
    'x_offset_step': 5,  # 右图中每条线的间距 (可视距离)
    'track_lat_limit': (-80, 80)  # 绘图时的纬度截断范围
}


def load_data(results_dir):
    """
    读取目录下所有的 .npz 文件
    返回格式: list of dicts [{'date': '2021-07-11', 'lat':..., 'lon':..., 'lgd':...}]
    """
    files = sorted(glob.glob(os.path.join(results_dir, "LGD_Result_*.npz")))
    data_list = []

    print(f"找到 {len(files)} 个结果文件。")

    for f in files:
        # 从文件名解析日期 (假设文件名格式 LGD_Result_2021-07-11.npz)
        basename = os.path.basename(f)
        date_str = basename.replace("LGD_Result_", "").replace(".npz", "")

        try:
            raw = np.load(f)
            # 简单过滤一下纬度范围，避免画太长
            mask = (raw['lat'] >= CONFIG['track_lat_limit'][0]) & (raw['lat'] <= CONFIG['track_lat_limit'][1])

            item = {
                'date_obj': datetime.strptime(date_str, "%Y-%m-%d"),
                'date_str': date_str,
                'lat': raw['lat'][mask],
                'lon': raw['lon'][mask],
                'lgd': raw['lgd'][mask]
            }
            data_list.append(item)
        except Exception as e:
            print(f"加载 {f} 失败: {e}")

    # 按日期排序
    data_list.sort(key=lambda x: x['date_obj'])
    return data_list


def plot_combined_figure(data_list):
    if not data_list:
        print("没有数据可绘图。")
        return

    # 创建画布：1行2列，宽度比例 1:1.5
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.15)

    # ==========================================
    # Panel A: Map (Orbit Tracks)
    # ==========================================
    ax1 = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())

    # 添加地图要素
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax1.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    ax1.add_feature(cfeature.LAND, facecolor='#f0f0f0')  # 浅灰色陆地
    ax1.add_feature(cfeature.OCEAN, facecolor='#e0f3f8')  # 浅蓝色海洋

    # 设置地图范围
    ext = CONFIG['map_extent']
    ax1.set_extent(ext, crs=ccrs.PlateCarree())

    # 绘制经纬度网格
    gl = ax1.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}

    # 1. 绘制红色矩形框 (Target Area)
    box = CONFIG['region_box']  # lon_min, lon_max, lat_min, lat_max
    rect = mpatches.Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2],
                              linewidth=2, edgecolor='red', facecolor='none',
                              transform=ccrs.PlateCarree(), zorder=10)
    ax1.add_patch(rect)

    # 2. 绘制所有轨迹
    for d in data_list:
        ax1.plot(d['lon'], d['lat'], color='red', linewidth=0.8, alpha=0.4, transform=ccrs.PlateCarree())

    # 添加 "A" 标签
    ax1.text(0.02, 0.95, 'A', transform=ax1.transAxes, fontsize=20, fontweight='bold', va='top')
    ax1.set_title("Orbit Tracks Passing Target Area", fontsize=14)

    # ==========================================
    # Panel B: Waterfall Plot (LGD Profiles)
    # ==========================================
    ax2 = fig.add_subplot(gs[1])

    # 获取目标区域的纬度范围，用于绘制灰色背景带
    lat_min, lat_max = CONFIG['region_box'][2], CONFIG['region_box'][3]
    ax2.axhspan(lat_min, lat_max, color='lightgray', alpha=0.3, zorder=0)

    # 遍历数据绘图
    step = CONFIG['x_offset_step']
    colors = plt.cm.tab10.colors  # 使用一个颜色循环

    # 记录月份位置，用于绘制顶部的月份标签
    month_labels = {}

    for i, d in enumerate(data_list):
        # 计算 X 轴偏移量
        x_shift = i * step

        # 数据转换
        lgd_nm = d['lgd'] * CONFIG['lgd_scale']

        # 绘图 (X = 值 + 偏移, Y = 纬度)
        color = colors[i % len(colors)]
        ax2.plot(lgd_nm + x_shift, d['lat'], color=color, linewidth=2, alpha=0.8)

        # 在顶部绘制日期 (Day + 'th')
        day_str = d['date_obj'].strftime("%d")
        # 简单的序数词后缀处理 (1st, 2nd, 3rd, 4th...)
        if 10 <= int(day_str) % 100 <= 20:
            suffix = 'th'
        else:
            suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(int(day_str) % 10, 'th')

        ax2.text(x_shift, CONFIG['track_lat_limit'][1] + 1, f"{int(day_str)}{suffix}",
                 ha='center', va='bottom', fontsize=9, rotation=0)

        # 记录月份信息 (如果不重复则记录)
        month_str = d['date_obj'].strftime("%B %Y")
        if month_str not in month_labels:
            month_labels[month_str] = x_shift  # 记录该月份出现的第一个位置
        else:
            # 更新位置为该月份所有轨迹的中心（可选优化，这里简单记录起始点）
            pass

        # 绘制竖直分割线 (可选)
        ax2.axvline(x_shift + step / 2 + step * 0.5, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    # 绘制月份标签 (顶部第二层)
    for m_str, x_pos in month_labels.items():
        # 稍微向右偏一点作为月份的起始
        ax2.text(x_pos, CONFIG['track_lat_limit'][1] + 4, m_str,
                 ha='left', va='bottom', fontsize=12, color='darkblue', fontweight='bold')

    # 设置轴标签和范围
    ax2.set_ylabel("Latitude (deg)", fontsize=12)
    ax2.set_xlabel("LGD ($nm/s^2$)", fontsize=12)
    ax2.set_ylim(CONFIG['track_lat_limit'])

    # X轴刻度处理：因为是错位图，物理刻度意义不明显，通常隐藏或仅保留底部的相对刻度
    # 这里我们模拟参考图，隐藏顶部和右侧边框，保留底部刻度但其实际数值含义是 "Offset + Value"
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # 我们可以设置X轴显示的范围，让它看起来紧凑
    ax2.set_xlim(-step, len(data_list) * step + step)

    # 添加 "B" 标签
    ax2.text(0.02, 0.95, 'B', transform=ax2.transAxes, fontsize=20, fontweight='bold', va='top')

    start_date = data_list[0]['date_str']
    end_date = data_list[-1]['date_str']
    main_title = f"GLDAS NOAH SoilMoisture: Line-of-Sight Gravity Difference (LGD)\nAnalysis Period: {start_date} to {end_date}"
    fig.suptitle(main_title, fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout()
    output_filename = f'LGD_Plot_GLDAS_NOAH_SoilMoisture_{start_date}_{end_date}.png'
    save_path = os.path.join(CONFIG['results_dir'], output_filename)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"💾 图形已保存: {save_path}")
    plt.show()


# ================= 主程序 =================
if __name__ == "__main__":

    # 2. 加载数据
    data = load_data(CONFIG['results_dir'])

    # 3. 绘图
    plot_combined_figure(data)