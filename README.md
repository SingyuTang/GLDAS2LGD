# GLDAS2LGD



# 1.程序概述

本程序旨在实现 Ghobadi-Far et al. (2022) 在 Journal of Geophysical Research: Solid Earth 发表的论文 "Along-Orbit Analysis of GRACE Follow-On Inter-Satellite Laser Ranging Measurements for Sub-Monthly Surface Mass Variations" 中的部分核心算法。

主要功能是利用 GLDAS (Global Land Data Assimilation System) 水文模型数据（土壤湿度），计算沿 GRACE-FO 卫星轨道的 视线方向重力差 (Line-of-Sight Gravity Difference, LGD) 模拟值。这可以用于验证 GRACE-FO LRI (Laser Ranging Interferometer) 观测到的高频（亚月级）质量变化信号。



# 2.数据准备

## 2.1 **GLDAS 水文数据**



***\*数据产品：\**** GLDAS Noah Land Surface Model L4 3 hourly 0.25 x 0.25 degree V2.1

***\*下载地址：\****[NASA Earthdata (GES DISC) ](https://search.earthdata.nasa.gov/search/granules?p=C1342986035-GES_DISC&pg[0][v]=f&pg[0][gsk]=-start_date&q=NOAH025)

***\*文件格式：\**** GLDAS_CLSM025_DA1_D.A20200501.022.nc4

***\*存放结构：\**** 建议按年或月存放，I:\LGD\GLDAS_CLSM025_D_2.2_2020

## 2.2 GRACE-FO轨道数据



程序需要 GRACE-FO 的精密轨道数据（GNV1B）来确定卫星在特定时间的位置。代码中通过 OrbitLoader 类加载 GROOPS 工作区数据，你需要根据实际情况修改 `load_orbit` 函数以适配你的轨道文件格式（如 GNV1B等）。

如果按照博客[详细介绍利用GRACE1B多日数据计算LGD工作流程二_基于LRI1B多日数据](https://singyutang.github.io/2025/11/10/%E8%AF%A6%E7%BB%86%E4%BB%8B%E7%BB%8D%E5%88%A9%E7%94%A8GRACE1B%E5%A4%9A%E6%97%A5%E6%95%B0%E6%8D%AE%E8%AE%A1%E7%AE%97LGD%E5%B7%A5%E4%BD%9C%E6%B5%81%E7%A8%8B%E4%BA%8C-%E5%9F%BA%E4%BA%8ELRI1B%E5%A4%9A%E6%97%A5%E6%95%B0%E6%8D%AE/)已经创建了工作根目录`workdir`并进行了相关处理，只需要将本项目文件全部拷贝到该目录下即可运行，因为本项目中的 `load_orbit` 函数默认读取`workdir/gracefo_dataset`路径下的GNV1B轨道文件。



# 3.用户指南



***\*第一步：配置路径与参数\****



打开 `lgd_processor.py`，修改 `main` 函数中的 `CONFIG` 字典与待处理日期列表 `date_list`：



\```

CONFIG = {

​    \# 路径设置

​    'gldas_dir': r"I:\LGD\GLDAS_NOAH025_3H_2.1_2020",

​    'groops_workspace': r'G:\GROOPS\PNAS2020Workspace',

​    'output_dir': r"./results_lgd",



​    \# 背景场设置 (用于计算距平)

​    'bgd_start': '2020-05-01',

​    'bgd_end': '2020-05-31',



​    \# 区域设置

​    'region_lon': (88, 92),

​    'region_lat': (22, 26),

​    'track_lat_limit': (-80.0, 80.0),



​    \# 轨道与计算参数

​    'orbit_direction': 'asc',  # 'asc' or 'desc'

​    'orbit_interval': 5,  # 秒

​    'cutoff_deg': 20.0,  # 积分截断距离

​    'n_max': 200,  # 阶数



​    \# 绘图参数

​    'results_dir': r"./results_lgd",  # 上一步代码保存结果的文件夹

​    'track_lat_limit': (-80.0, 80.0)  # 绘图时的纬度截断范围

  }



date_list = [

​    '2020-06-04', '2020-06-10', '2020-06-15', '2020-06-21', '2020-06-26',

​    '2020-07-02', '2020-07-07', '2020-07-13', '2020-07-18', '2020-07-24', '2020-07-29',

​    '2020-08-04', '2020-08-09', '2020-08-15', '2020-08-20'

  ]

\```



***\*第二步：运行计算程序\****



运行 lgd_processor.py。程序将执行以下操作：

\1. 读取 GLDAS 数据计算月平均背景场。

\2. 针对 date_list 中的每一天，提取轨道并计算瞬时 LGD。

\3. 打印进度并在 output_dir 生成 .npz 文件。



\```bash



python lgd_processor.py

\```



***\*第三步：运行绘图程序\****



计算完成后，打开  `lgd_plot.py`，确保 `CONFIG` 中的 `results_dir` 指向刚才的输出目录，然后运行：



\```python



python lgd_plot.py

\```



程序将生成 PNG 图片，展示目标区域内 LGD 随时间的变化。



# 4.参考文章



Ghobadi-Far, Khosro, Shin-Chan Han, Christopher M. McCullough, David N. Wiese, Richard D. Ray, Jeanne Sauber, Linus Shihora, and Henryk Dobslaw. 2022. “Along-Orbit Analysis of GRACE Follow-On Inter-Satellite Laser Ranging Measurements for Sub-Monthly Surface Mass Variations.” Journal of Geophysical Research: Solid Earth 127(2):e2021JB022983. doi:10.1029/2021JB022983.
