from gldas2lgd import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
import glob

def get_gldas_filename(root_dir, date_str, time_str="0000"):
    """根据日期生成GLDAS文件名"""
    # 格式: GLDAS_NOAH025_3H.A20210711.0000.021.nc4
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    date_part = dt.strftime("A%Y%m%d")
    filename = f"GLDAS_NOAH025_3H.{date_part}.{time_str}.021.nc4"
    return os.path.join(root_dir, filename)


def save_results(output_dir, date_str, lat, lon, lgd, alt):
    """保存计算结果到文件"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    save_path = os.path.join(output_dir, f"LGD_Result_{date_str}.npz")
    np.savez(save_path, lat=lat, lon=lon, lgd=lgd, alt=alt)
    print(f"结果已保存至: {save_path}")


def process_single_date(calc, date_str, mean_soilmoi, config):
    """
    处理单个日期的核心逻辑
    :param calc: TrackLGDCalculator 实例
    :param date_str: 处理日期 'YYYY-MM-DD'
    :param mean_soilmoi: 预计算好的背景场土壤湿度 (2D array)
    :param config: 配置字典
    :return: (lats, lgd_values) 或 (None, None) 如果失败
    """
    print(f"\n[{date_str}] 开始处理...")

    # 1. 读取当日 GLDAS 数据
    gldas_file = get_gldas_filename(config['gldas_dir'], date_str)
    if not os.path.exists(gldas_file):
        print(f"[{date_str}] 错误: GLDAS文件不存在 -> {gldas_file}")
        return None, None

    try:
        ds = calc.read_gldas_smc(gldas_file)
    except Exception as e:
        print(f"[{date_str}] 读取GLDAS出错: {e}")
        return None, None

    # 准备网格数据
    gldas_lat_range, gldas_lon_range = np.array(ds['lat']), np.array(ds['lon'])
    lon_grid, lat_grid = np.meshgrid(gldas_lon_range, gldas_lat_range)

    # 计算异常值 (当前值 - 背景均值)
    current_soilmoi = np.squeeze(ds['SoilMoi'])

    # 确保维度一致
    if current_soilmoi.shape != mean_soilmoi.shape:
        print(f"[{date_str}] 维度不匹配: Current {current_soilmoi.shape} vs Bgd {mean_soilmoi.shape}")
        return None, None

    h_grid = (current_soilmoi - mean_soilmoi) * 1e-3  # 单位转为 m

    # 展平
    g_h = h_grid.flatten()
    g_lat = lat_grid.flatten()
    g_lon = lon_grid.flatten()

    # 2. 加载轨道数据
    try:
        # 注意：这里加载的是当天的轨道
        orbitc, orbitd = load_orbit(date_str=date_str,
                                    groops_workspace=config['groops_workspace'],
                                    coord_type='geodetic')
    except Exception as e:
        print(f"[{date_str}] 加载轨道数据出错: {e}")
        return None, None

    # 采样与提取
    interval = config['orbit_interval']
    posc_geo = np.array([obj.position for obj in orbitc])[::interval]
    posd_geo = np.array([obj.position for obj in orbitd])[::interval]

    # 3. 筛选经过研究区的轨迹
    lonlat = np.column_stack([posc_geo[:, 0], posd_geo[:, 1]])

    tracks, indices = filter_complete_tracks_passing_region(
        lonlat,
        config['region_lon'],
        config['region_lat'],
        lat_limit=config['track_lat_limit'],
        separate=False,
        direction=config['orbit_direction']
    )

    if indices is None or len(indices) == 0:
        print(f"[{date_str}] 未找到经过目标区域的有效轨迹。")
        return None, None

    # 提取筛选后的轨迹点
    track_s1 = posc_geo[np.squeeze(indices)]
    track_s2 = posd_geo[np.squeeze(indices)]

    # 更新 calc 对象中的轨迹属性
    calc.track_s1_lon, calc.track_s1_lat, calc.track_s1_h = track_s1[:, 0], track_s1[:, 1], track_s1[:, 2]
    calc.track_s2_lon, calc.track_s2_lat, calc.track_s2_h = track_s2[:, 0], track_s2[:, 1], track_s2[:, 2]

    # 4. 执行 LGD 计算
    t_start = time.time()
    lgd_series = calc.compute_track_lgd(
        calc.track_s1_lat,
        calc.track_s1_lon,
        g_h, g_lat, g_lon,
        cutoff_deg=config['cutoff_deg']
    )
    print(f"[{date_str}] 计算完成，耗时 {time.time() - t_start:.2f} 秒")

    # 保存结果
    save_results(config['output_dir'], date_str, calc.track_s1_lat, calc.track_s1_lon, lgd_series, calc.track_s1_h)

    return calc.track_s1_lat, lgd_series


def main():
    # ================= 配置区域 =================
    CONFIG = {
        # 路径设置
        'gldas_dir': r"I:\LGD\GLDAS_NOAH025_3H_2.1_2020",
        'groops_workspace': r'G:\GROOPS\PNAS2020Workspace',
        'output_dir': r"./results_lgd",

        # 背景场设置 (用于计算距平)
        'bgd_start': '2020-05-01',
        'bgd_end': '2020-05-31',

        # 区域设置
        'region_lon': (88, 92),
        'region_lat': (22, 26),
        'track_lat_limit': (-80.0, 80.0),

        # 轨道与计算参数
        'orbit_direction': 'asc',  # 'asc' or 'desc'
        'orbit_interval': 5,  # 秒
        'cutoff_deg': 20.0,  # 积分截断距离
        'n_max': 200,  # 阶数

        # 绘图参数
        'results_dir': r"./results_lgd",  # 上一步代码保存结果的文件夹
        'track_lat_limit': (-80.0, 80.0)  # 绘图时的纬度截断范围
    }

    # ================= 待处理日期列表 =================
    date_list = [
        '2020-06-04', '2020-06-10', '2020-06-15', '2020-06-21', '2020-06-26',
        '2020-07-02', '2020-07-07', '2020-07-13', '2020-07-18', '2020-07-24', '2020-07-29',
        '2020-08-04', '2020-08-09', '2020-08-15', '2020-08-20'
    ]

    print(f"计划处理日期: {date_list}")

    # ================= 初始化 =================
    calc = TrackLGDCalculator(n_max=CONFIG['n_max'])

    # 1. 预计算背景场 (只做一次)
    print(f"\n>>> 正在计算背景场 ({CONFIG['bgd_start']} 到 {CONFIG['bgd_end']})...")
    try:
        bgd_data = calc.calculate_smc_mean_by_date_range(
            CONFIG['gldas_dir'],
            CONFIG['bgd_start'],
            CONFIG['bgd_end']
        )
        mean_soilmoi = np.squeeze(bgd_data['mean_SoilMoi'])
        print("背景场计算完成。")
    except Exception as e:
        print(f"背景场计算失败，程序终止: {e}")
        return

    # ================= 批量循环 =================
    success_count = 0
    results_summary = {}  # 用于存储最后简单画图用

    for date_str in date_list:
        lats, lgds = process_single_date(calc, date_str, mean_soilmoi, CONFIG)
        if lats is not None:
            success_count += 1
            results_summary[date_str] = (lats, lgds)

    print(f"\n==========================================")
    print(f"批量处理结束。成功: {success_count}/{len(date_list)}")
    print(f"结果已保存至: {CONFIG['output_dir']}")


if __name__ == "__main__":
    main()