# AIR5023 Pure Pursuit 路径跟踪仿真

本仓库提供小组项目 M4 的 Python 仿真，用于 M5 的前视距离对照实验和 M6 的路径几何对照实验。所有工况使用同一个 Pure Pursuit 控制器、运动学自行车模型和前向欧拉积分器。

## 安装与试跑

需要 Python 3.10 或以上。进入克隆后的仓库目录，创建并激活虚拟环境。

```powershell
# Windows PowerShell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate
```

```bash
python -m pip install -r requirements.txt
python simulation.py --path s_curve --lookahead 2.0 --output outputs/demo
```

如果 PowerShell 不允许运行激活脚本，可在后续命令中用 `.\.venv\Scripts\python.exe` 代替 `python`，直接使用虚拟环境的解释器。

## 复现实验

```bash
# M5 — Effect of Lookahead Distance
python simulation.py --path s_curve --lookahead 0.5 1.0 2.0 3.0 4.0 --output outputs/m5/s_curve

# M6 — Effect of Path Geometry
python simulation.py --path straight --lookahead 0.5 2.0 4.0 --output outputs/m6/straight
python simulation.py --path circle --lookahead 0.5 2.0 4.0 --output outputs/m6/circle
python simulation.py --path s_curve --lookahead 0.5 2.0 4.0 --output outputs/m6/s_curve
```

M5 和 M6 中相同的 S 形工况会得到相同的结果。每条命令的输出放在独立目录，方便与报告中的图表对应。

| 输出 | 内容 |
| --- | --- |
| `summary.csv` | 每组工况的 RMSE、最大跟踪误差、转向变化率、完成状态和用时 |
| `trajectory_<distance>m.csv` | 逐步位姿、误差、路径进度、目标点、实际目标距离和转向角 |
| `tracking.png` | 参考与实际轨迹、跟踪误差曲线和转向角曲线 |

[示例数据](example/summary.csv)和[示例图](example/tracking.png)来自一次 S 形路径仿真。配置前视半径为 2.0 m 时，完成用时为 36.15 s，RMSE 为 0.1202 m。

## 报告用参数与指标定义

- 所有工况使用速度 1.0 m/s、轴距 1.0 m、步长 0.05 s 和 ±45° 转向限幅。初始航向沿路径切线，后轴中心从切线左侧 0.5 m 开始。
- 直线路径长 32 m，圆形路径半径为 5 m，S 形路径由两个正弦周期组成，振幅为 2 m，波长为 16 m。每条路径使用 801 个采样点和相邻点之间的线段。
- 跟踪误差是后轴中心到参考折线的最短距离，RMSE 和最大误差包含初始与终止状态。转向变化率是相邻实际转向命令的绝对变化量除以步长后取平均，单位为 °/s。逐步 CSV 的最后一行没有目标点或转向命令，因为该时刻不再积分。
- 开放路径要求到终点的欧氏距离和剩余弧长均不超过 0.25 m。圆形路径要求完整跟踪一圈，且终止时跟踪误差不超过 0.25 m。所有运行的时间上限为 60 s。

几何控制方法见 Coulter (1992), *Implementation of the Pure Pursuit Path Tracking Algorithm*, CMU-RI-TR-92-01. [技术报告](https://publications.ri.cmu.edu/storage/publications/pub_files/pub3/coulter_r_craig_1992_1/coulter_r_craig_1992_1.pdf)
