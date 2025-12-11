import numpy as np
from scipy.special import lpn
import time
import os
import netCDF4 as nc
from datetime import datetime, timedelta
from S02compute_grace_lgd import OrbitLoader
from S05plot_lgd_ra_cwt_filter import filter_complete_tracks_passing_region
import read_love_numbers
import matplotlib.pyplot as plt


def load_orbit(date_str="2021-08-01", groops_workspace=r'G:\GROOPS\WRR2021WorkspaceLRI', coord_type='geodetic'):
    orbit_loader = OrbitLoader(date_str=date_str, groops_workspace_dir=groops_workspace)
    orbitc = orbit_loader.load_orbit_data('gnv1b', 'C', coord_type)
    orbitd = orbit_loader.load_orbit_data('gnv1b', 'D', coord_type)
    return orbitc, orbitd

def crop_grid_safe(h_grid, lat_grid, lon_grid, lat_range, lon_range):
        """
        安全地根据给定的经纬度范围裁剪网格

        参数:
            h_grid: 等效水高网格 (2D数组)
            lat_grid: 纬度网格 (2D数组)
            lon_grid: 经度网格 (2D数组)
            lat_range: (lat_min, lat_max) 纬度范围
            lon_range: (lon_min, lon_max) 经度范围

        返回:
            h_cropped, lat_cropped, lon_cropped: 裁剪后的网格
        """
        lat_min, lat_max = lat_range
        lon_min, lon_max = lon_range

        # 创建经纬度掩码
        lat_mask = (lat_grid >= lat_min) & (lat_grid <= lat_max)
        lon_mask = (lon_grid >= lon_min) & (lon_grid <= lon_max)

        # 组合掩码
        mask = lat_mask & lon_mask

        # 找到包含所有True值的矩形区域
        rows = np.any(mask, axis=1)
        cols = np.any(mask, axis=0)

        row_indices = np.where(rows)[0]
        col_indices = np.where(cols)[0]

        if len(row_indices) == 0 or len(col_indices) == 0:
            raise ValueError("在给定的经纬度范围内没有找到数据点")

        # 裁剪网格
        row_start, row_end = row_indices[0], row_indices[-1] + 1
        col_start, col_end = col_indices[0], col_indices[-1] + 1

        h_cropped = h_grid[row_start:row_end, col_start:col_end]
        lat_cropped = lat_grid[row_start:row_end, col_start:col_end]
        lon_cropped = lon_grid[row_start:row_end, col_start:col_end]

        print(f"裁剪后网格大小: {h_cropped.shape}")
        print(f"纬度范围: [{lat_cropped.min():.2f}, {lat_cropped.max():.2f}]")
        print(f"经度范围: [{lon_cropped.min():.2f}, {lon_cropped.max():.2f}]")

        return h_cropped, lat_cropped, lon_cropped


class TrackLGDCalculator:
    def __init__(self, n_max=200):
        self.G = 6.67430e-11
        self.M_EARTH = 5.972e24
        self.R_EARTH = 6378137.0
        self.RHO_WATER = 1000.0
        self.n_max = n_max
        self.k_l = self._init_load_love_numbers(n_max)

        self.track_s1_lon = None    # 卫星C经过研究区并延伸后的经度数组（degree）
        self.track_s2_lon = None    # 卫星D经过研究区并延伸后的经度数组（degree）
        self.track_s1_lat = None    # 卫星C经过研究区并延伸后的经度数组（degree）
        self.track_s2_lat = None    # 卫星D经过研究区并延伸后的纬度数组（degree）
        self.track_s1_h = None      # 卫星C经过研究区并延伸后的大地高数组（m）
        self.track_s2_h = None      # 卫星D经过研究区并延伸后的大地高数组（m）


    def _init_load_love_numbers(self, n_max):
        """
        加载Love数

        # PREM outputs from Han and Wahr (1995)
        # https://doi.org/10.1111/j.1365-246X.1995.tb01819.x

        :param lmax: 最大Legendre阶数
        :return: Love数-kl
        """
        love_numbers = read_love_numbers.load_love_numbers(n_max, LOVE_NUMBERS=0, REFERENCE='CF', FORMAT='tuple')
        return np.array(love_numbers[1])

    def _calculate_gravity_contribution(self, sat_sph, mass_lat_rad, mass_lon_rad, dm_vals):
        """核心积分内核：计算一组质量点对单个卫星位置的 NEU 引力"""
        r_i, phi_i, lam_i = sat_sph
        phi_j, lam_j = mass_lat_rad, mass_lon_rad

        # 1. 计算 Cos(psi)
        cos_psi = np.sin(phi_i) * np.sin(phi_j) + \
                  np.cos(phi_i) * np.cos(phi_j) * np.cos(lam_i - lam_j)
        cos_psi = np.clip(cos_psi, -1.0, 1.0)

        # 2. 勒让德多项式
        P_l, P_prime_l = lpn(self.n_max, cos_psi)

        # 3. 几何偏导
        dcos_dphi = np.cos(phi_i) * np.sin(phi_j) - np.sin(phi_i) * np.cos(phi_j) * np.cos(lam_i - lam_j)
        dcos_dlam = -np.cos(phi_i) * np.cos(phi_j) * np.sin(lam_i - lam_j)

        # 4. 级数求和 (针对这一组点)
        # 为了速度，这里对 n_max 循环，利用 numpy 的广播对所有点同时计算
        sum_N = np.zeros_like(dm_vals)
        sum_E = np.zeros_like(dm_vals)
        sum_U = np.zeros_like(dm_vals)

        ratio = self.R_EARTH / r_i
        powers = np.power(ratio, np.arange(self.n_max + 3))

        for l in range(self.n_max + 1):
            factor = (1 + self.k_l[l]) * powers[l + 2]

            # P_prime_l[l] 形状是 (N_points,)
            term_n = factor * dcos_dphi * P_prime_l[l]
            sum_N += term_n

            if abs(np.cos(phi_i)) > 1e-9:
                term_e = factor * (dcos_dlam / np.cos(phi_i)) * P_prime_l[l]
                sum_E += term_e

            term_u = -factor * (l + 1) * P_l[l]
            sum_U += term_u

        # 5. 乘以各自的质量
        pre_factor = (self.G * dm_vals) / (self.R_EARTH ** 2)

        # 累加得到总 NEU
        g_N = np.sum(sum_N * pre_factor)
        g_E = np.sum(sum_E * pre_factor)
        g_U = np.sum(sum_U * pre_factor)

        return np.array([g_N, g_E, g_U])

    def neu_to_ecef(self, g_neu, lat, lon):
        sp, cp = np.sin(lat), np.cos(lat)
        sl, cl = np.sin(lon), np.cos(lon)
        gN, gE, gU = g_neu
        dx = -sp * cl * gN - sl * gE + cp * cl * gU
        dy = -sp * sl * gN + cl * gE + cp * sl * gU
        dz = cp * gN + sp * gU
        return np.array([dx, dy, dz])

    def compute_track_lgd(self, track_lats, track_lons, gldas_h, gldas_lat, gldas_lon, cutoff_deg=25.0):
        """
        计算整条轨迹的 LGD

        参数:
            track_lats: 卫星轨迹纬度数组
            track_lons: 卫星轨迹经度数组
            gldas_h:    全球水高数据 (Flattened 1D array)
            gldas_lat:  全球网格纬度 (Flattened 1D array)
            gldas_lon:  全球网格经度 (Flattened 1D array)
            cutoff_deg: 截断距离 (度)，只计算卫星周围多少度以内的质量

        return
            results_lgd: np.ndarray， unit -> m
        """
        results_lgd = []

        # 预计算全球网格的质量 (Delta Mass)
        # 假设分辨率 1度 (这里为了演示简化，实际应根据输入分辨率算)
        res_deg = 0.25
        d_rad = np.radians(res_deg)
        # 面积元 = R^2 * d_phi * d_lam * cos(lat)
        area_arr = (self.R_EARTH ** 2) * (d_rad ** 2) * np.cos(np.radians(gldas_lat))
        dm_global = gldas_h * self.RHO_WATER * area_arr

        # 只保留非零质量点 (稀疏优化)
        valid_mask = np.abs(dm_global) > 1.0
        g_lat_v = np.radians(gldas_lat[valid_mask])
        g_lon_v = np.radians(gldas_lon[valid_mask])
        g_dm_v = dm_global[valid_mask]

        print(f"全球有效(非零)网格点数: {len(g_dm_v)}")

        # 遍历轨迹上的每一个点
        for i in range(len(track_lats)):
            lat_deg1, lat_deg2 = self.track_s1_lat[i], self.track_s2_lat[i]
            lon_deg1, lon_deg2 = self.track_s1_lon[i], self.track_s2_lon[i]
            r_s1, r_s2 = self.track_s1_h[i] + self.R_EARTH, self.track_s2_h[i] + self.R_EARTH

            # --- 1. 卫星位置 (Sat1 和 Sat2) ---
            lat_s1, lon_s1 = np.radians(lat_deg1), np.radians(lon_deg1)
            lat_s2, lon_s2 = np.radians(lat_deg2), np.radians(lon_deg2)

            # --- 2. 空间索引 (Spatial Indexing / Culling) ---
            # 核心优化：只找出距离 Sat2 中心 cutoff_deg 范围内的网格点
            # 利用球面距离公式
            cos_dist = np.sin(lat_s2) * np.sin(g_lat_v) + \
                       np.cos(lat_s2) * np.cos(g_lat_v) * np.cos(lon_s2 - g_lon_v)

            # 阈值筛选 (cos(25度) approx 0.906)
            cos_cutoff = np.cos(np.radians(cutoff_deg))
            local_mask = cos_dist > cos_cutoff

            # 提取局部网格数据
            local_lats = g_lat_v[local_mask]
            local_lons = g_lon_v[local_mask]
            local_dms = g_dm_v[local_mask]
            if len(local_dms) == 0:
                results_lgd.append(0.0)
                continue

            # --- 3. 计算 Sat1 和 Sat2 的引力 ---
            s1_sph = (r_s1, lat_s1, lon_s1)
            s2_sph = (r_s2, lat_s2, lon_s2)

            g_neu_1 = self._calculate_gravity_contribution(s1_sph, local_lats, local_lons, local_dms)
            g_neu_2 = self._calculate_gravity_contribution(s2_sph, local_lats, local_lons, local_dms)

            # 转 ECEF
            g_ecef_1 = self.neu_to_ecef(g_neu_1, lat_s1, lon_s1)
            g_ecef_2 = self.neu_to_ecef(g_neu_2, lat_s2, lon_s2)

            # --- 4. 投影 LGD ---
            pos1 = self.neu_to_ecef([0, 0, r_s1], lat_s1, lon_s1)  # 简化位置向量
            pos2 = self.neu_to_ecef([0, 0, r_s2], lat_s2, lon_s2)
            e12 = (pos2 - pos1) / np.linalg.norm(pos2 - pos1)

            lgd = np.dot((g_ecef_2 - g_ecef_1), e12)
            results_lgd.append(lgd)

            if i % 10 == 0:
                print(f"进度: Lat {lat_deg1:.1f}, 局部汇聚点数: {len(local_dms)}, LGD: {lgd * 1e9:.2f} nm/s^2")

        return np.array(results_lgd)

    # 读取dirpath路径下所有文件名
    def get_file_list(self, directory, suffix='.nc4'):
        file_list = []
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith(suffix):
                    file_list.append(os.path.join(root, file))
        return file_list

    # 获取文件列表中含有给定字符串的文件索引
    def get_file_index(self, file_list, string):
        index_list = []
        for i, file in enumerate(file_list):
            if string in file:
                index_list.append(i)
        return index_list

    # 读取GLDAS NOAHT025_3H.nc文件的土壤湿度数据
    def read_gldas_smc(self, file_path):
        """
            读取GLDAS NOAHT025_3H.nc文件的土壤湿度数据。
        :param file_path: str
            GLDAS NOAHT025_3H.nc文件路径
        :return: xr.Dataset: 包含以下内容的数据集
                - SoilMoi: 各土层含水量之和 [kg/m²]
                 - 坐标维度: time, lat, lon
        """

        def safe_read_nc_variable(dataset, var_name, default_fill=-9999):
            """
            安全读取netCDF变量并清理数据

            Args:
                dataset: netCDF数据集
                var_name: 变量名
                default_fill: 默认填充值（如果没有找到属性）
            """
            var = dataset.variables[var_name]
            data = var[:]

            # 优先使用文件中的填充值属性
            fill_value = None
            if hasattr(var, '_FillValue'):
                fill_value = var._FillValue
            elif hasattr(var, 'missing_value'):
                fill_value = var.missing_value
            else:
                fill_value = default_fill

            # 处理掩码数组
            if np.ma.is_masked(data):
                data = data.filled(np.nan)

            # 替换无效值
            if fill_value is not None:
                # 浮点数比较
                if isinstance(fill_value, float):
                    data = np.where(np.isclose(data, fill_value, rtol=1e-6), 0, data)
                else:
                    data = np.where(data == fill_value, 0, data)

            # 替换NaN
            data = np.where(np.isnan(data), 0, data)

            return data
        # 检测file_path是否存在
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件{file_path}不存在。")

        try:
            # 使用netCDF4打开文件
            dataset = nc.Dataset(file_path, 'r')

            # 读取经纬度和时间
            lon = dataset.variables['lon'][:]
            lat = dataset.variables['lat'][:]
            time = dataset.variables['time'][:]

            # 读取各土层土壤湿度并求和
            smc_0_10 = safe_read_nc_variable(dataset, 'SoilMoi0_10cm_inst')
            smc_10_40 = safe_read_nc_variable(dataset, 'SoilMoi10_40cm_inst')
            smc_40_100 = safe_read_nc_variable(dataset, 'SoilMoi40_100cm_inst')
            smc_100_200 = safe_read_nc_variable(dataset, 'SoilMoi100_200cm_inst')

            # 关闭文件
            dataset.close()
            # 计算总土壤湿度
            smc = smc_0_10 + smc_10_40 + smc_40_100 + smc_100_200
            # 返回字典格式的数据
            result = {
                'SoilMoi': smc,
                'lat': lat,
                'lon': lon,
                'time': time
            }
            return result

        except Exception as e:
            print(f"读取文件时出错: {e}")
            raise

    def calculate_smc_mean_by_date_range(self, data_dir, start_date, end_date):
        """
        批量读取指定日期范围内的所有土壤湿度数据并计算平均值

        :param data_dir: str, 数据目录路径
        :param start_date: str, 开始日期 'YYYY-MM-DD'
        :param end_date: str, 结束日期 'YYYY-MM-DD'
        :return: dict, 包含平均土壤湿度和坐标信息
        """

        # 转换日期字符串为datetime对象
        start = datetime.strptime(start_date, '%Y-%m-%d')
        end = datetime.strptime(end_date, '%Y-%m-%d')

        # 存储所有读取的数据
        all_smc_data = []
        valid_files = []

        current_date = start
        while current_date <= end:
            # 生成每天的8个时间点
            for hour in [0, 3, 6, 9, 12, 15, 18, 21]:
                date_str = current_date.strftime("A%Y%m%d")
                time_str = f"{hour:02d}00"
                filename = f"GLDAS_NOAH025_3H.{date_str}.{time_str}.021.nc4"
                file_path = os.path.join(data_dir, filename)

                if os.path.exists(file_path):
                    try:
                        data = self.read_gldas_smc(file_path)
                        all_smc_data.append(data['SoilMoi'])
                        valid_files.append(filename)
                    except Exception as e:
                        print(f"读取文件 {filename} 时出错: {e}")

            current_date += timedelta(days=1)

        if not all_smc_data:
            raise ValueError(f"在目录 {data_dir} 中未找到指定日期范围内的数据文件")

        # 计算平均值
        print(f"开始计算土壤湿度格网平均值，共 {len(all_smc_data)} 个文件...")
        smc_array = np.stack(all_smc_data, axis=0)
        mean_smc = np.mean(smc_array, axis=0)

        # 获取坐标信息
        lat = data['lat']
        lon = data['lon']

        return {
            'mean_SoilMoi': mean_smc,
            'lat': lat,
            'lon': lon,
            'time_period': f"{start_date} to {end_date}",
            'num_files': len(all_smc_data),
            'file_list': valid_files
        }
