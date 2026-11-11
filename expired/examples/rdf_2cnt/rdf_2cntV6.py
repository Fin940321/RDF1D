"""
Based on rdf_2cntV7.py, but with the RDF center/target roles swapped:
  - V7: functional group atoms (near electrode) → electrode carbons
  - V8: electrode carbons (CNT or Junction)    → functional group atoms (near electrode)

Classification logic is unchanged: functional group atoms are still classified
by their proximity to CNT / Junction / Bulk before the RDF calculation.
"""
from MDAnalysis import *
import MDAnalysis.analysis.distances as distances
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from tqdm import tqdm

# --- 引入分析資料 ---
topology = "2V_start_drudes.pdb"
trajectory = "2V_FV_NVT.dcd"
electrode_pdb = "2V_start_nodrudes.pdb"  # 用於自動偵測電極位置
u = Universe(topology, trajectory)

# --- 參數設定區 ---
r_electrode_contact = 6.0   # 判斷離子是否靠近電極的距離閾值 (Å)
junction_cutoff     = 7.0  # 定義 Junction：距離石墨烯幾 Å 內的 CNT 原子 (Å)
z_max    = 65
rdf_max  = 15
n_bins   = 75
dr       = float(rdf_max) / float(n_bins)
framestart = len(u.trajectory)//2  # 從中間開始分析，避免初始幀的非平衡狀態
frameend   = len(u.trajectory)  # slice exclusive，涵蓋最後一幀

# --- 定義 5 個目標基團 (Target Functional Groups) ---
print("🎯 正在定義陰陽離子的 5 個目標基團...")

targets_tfsi = {
    'SO2_polar'     : "resname Tf2 and name Otf Otf1 Otf2 Otf3",
    'CF3_nonpolar'  : "resname Tf2 and name Ctf Ctf1"
}
targets_bmim = {
    'Im_polar'      : "resname BMI and name C1 C2 C21",
    'CH3_nonpolar'  : "resname BMI and name C3",
    'Butyl_nonpolar': "resname BMI and name C4 C5 C51 C6"
}
all_targets = {**targets_tfsi, **targets_bmim}

# --- 資料結構初始化 ---
# 3 個環境：CNT、Junction、Bulk（對照組）
environments_list = ['CNT', 'Junction', 'Bulk']

rdf_results = {
    env: {name: np.zeros(n_bins, dtype=np.float64) for name in all_targets}
    for env in environments_list
}

# 各環境的正規化分母：累積 pair 數
rdf_N_N = {
    env: {name: 0 for name in all_targets}
    for env in environments_list
}

# --- 電極環境定義（固定，只做一次）---
print("⏳ 正在定義電極環境（CNT vs Junction）...")
u.trajectory[0]

grp_graphene    = u.select_atoms("segid A")
cnt_all         = u.select_atoms("segid B or segid C")

# 1. 找出碳管上屬於 Junction 區域的原子 (靠近石墨烯的部分)
cnt_in_junction = u.select_atoms(
    f"(segid B or segid C) and around {junction_cutoff} group grp",
    grp=grp_graphene
)

# 2. 找出石墨烯上屬於 Junction 底部區域的原子 (被碳管夾在中間或壓在下方的那塊)
# 這裡使用相同的距離邏輯，找出靠近碳管的石墨烯原子
grp_in_junction = u.select_atoms(
    f"segid A and around {junction_cutoff} group cnt",
    cnt=cnt_all
)

# 3. 複合 Junction 中心：包含兩側管壁與底部地板
center_junction = cnt_in_junction + grp_in_junction

# 4. 純 CNT 表面：只包含遠離地板的碳管上半部
center_cnt = cnt_all - cnt_in_junction

n_center_cnt      = len(center_cnt)
n_center_junction = len(center_junction)

if n_center_cnt == 0:
    print("⚠️  警告：純 CNT 原子數為 0！")
if n_center_junction == 0:
    print("⚠️  警告：Junction 原子數為 0！")

print(f"✅ 電極定義完成！純 CNT 原子數: {n_center_cnt}, Junction 原子數: {n_center_junction}\n")

# --- 自動偵測電極位置 (from start_nodrudes.pdb) ---
print("=== 偵測電極位置 ===")
u_elec  = Universe(electrode_pdb)
elecgrp = u_elec.select_atoms("resname grpc")
if len(elecgrp) == 0:
    raise ValueError(f"找不到 resname 'grpc' 的電極原子，請確認 {electrode_pdb}")
z_elec        = elecgrp.positions[:, 2]
z_elec_left   = z_elec.min()   # 左側電極（正電極 / Anode）
z_elec_right  = z_elec.max()   # 右側電極（負電極 / Cathode）
print(f"  ✅ 電極原子數 (grpc): {len(elecgrp)}")
print(f"  ✅ 正電極 (Anode,  左側): Z = {z_elec_left:.2f} Å")
print(f"  ✅ 負電極 (Cathode, 右側): Z = {z_elec_right:.2f} Å")
print(f"  ✅ 電極間距: {z_elec_right - z_elec_left:.2f} Å")
print()

# 【加速】預先建立目標基團的 AtomGroup（原子數固定，只做一次）
print("⚡ 預先建立目標基團的 AtomGroup...")
atom_groups = {
    name: u.select_atoms(sel_str)
    for name, sel_str in all_targets.items()
}
for name, grp in atom_groups.items():
    print(f"   {name}: {len(grp)} 個原子")
print()

# --- 軌跡迴圈 ---
print(f"📊 開始分析軌跡: frame {framestart} → {frameend - 1}")
print(f"   總共約 {frameend - framestart} 幀\n")

n_frames = 0

for ts in tqdm(u.trajectory[framestart:frameend], desc="Processing frames", unit="frame", ncols=80):
    n_frames += 1

    for target_name, target_grp_full in atom_groups.items():

        # Z 軸截斷：只取電極附近的目標原子（加速用）
        z_coords   = target_grp_full.positions[:, 2]
        z_mask     = z_coords < (z_max + rdf_max)
        target_grp = target_grp_full[z_mask]

        if len(target_grp) == 0:
            continue

        # ==========================================
        # 【分類階段】環境分類：
        # 以官能基原子距離電極的遠近來判斷環境（與 V7 相同）
        # ==========================================
        dist_to_cnt = np.amin(
            distances.distance_array(
                target_grp.positions, center_cnt.positions,
                box=u.dimensions, backend='OpenMP'
            ), axis=1
        )
        dist_to_junc = np.amin(
            distances.distance_array(
                target_grp.positions, center_junction.positions,
                box=u.dimensions, backend='OpenMP'
            ), axis=1
        )

        # Junction 優先分類
        mask_junction = dist_to_junc <= r_electrode_contact
        mask_cnt      = (dist_to_cnt  <= r_electrode_contact) & (~mask_junction)
        mask_bulk     = (~mask_junction) & (~mask_cnt)

        env_data = {
            'CNT'     : (mask_cnt,      center_cnt),
            'Junction': (mask_junction, center_junction),
            'Bulk'    : (mask_bulk,     None)
        }

        # ==========================================
        # 【V8 核心修改】RDF 計算方向：
        #   中心 (reference)     = 電極碳原子 (CNT 或 Junction)
        #   目標 (configuration) = 官能基原子 (已分類至該環境的那些)
        #
        # V7 是反過來的：
        #   中心 = 官能基原子，目標 = 電極碳原子
        # ==========================================
        for env_name, (mask, electrode_grp) in env_data.items():
            atoms_in_env = target_grp[mask]

            if len(atoms_in_env) == 0:
                continue

            # Bulk 環境：不計算 RDF
            if env_name == 'Bulk':
                continue

            # 【V8】以電極碳為中心，官能基為目標
            _, dist_rdf = distances.capped_distance(
                electrode_grp.positions,   # reference    = 電極碳原子
                atoms_in_env.positions,    # configuration = 官能基原子（已篩選至該環境）
                max_cutoff=rdf_max,
                box=u.dimensions,
                backend='OpenMP'
            )

            if len(dist_rdf) == 0:
                continue

            counts, _ = np.histogram(dist_rdf, bins=n_bins, range=(0.0, rdf_max))
            rdf_results[env_name][target_name]  += counts
            rdf_N_N[env_name][target_name]      += len(dist_rdf)  # 累積 pair 數

print(f"\n✅ 軌跡運算完成！實際處理幀數: {n_frames}")

# ==========================================
# 正規化
#   g(r) = counts(r) / N_pairs / V_shell(r) * V_sphere
# ==========================================
print("\n⏳ 開始進行正規化與寫入 .dat 檔案...")

edges         = np.linspace(0.0, rdf_max, n_bins + 1)
r_centers     = (edges[:-1] + edges[1:]) / 2
shell_volumes = (4/3) * np.pi * (edges[1:]**3 - edges[:-1]**3)
volume_sphere = (4/3) * np.pi * rdf_max**3

for env_name in ['CNT', 'Junction']:
    for target_name in all_targets:
        filename  = f"2V_rdf_{env_name}_{target_name}_V9.dat"
        counts    = rdf_results[env_name][target_name]
        n_pairs   = rdf_N_N[env_name][target_name]

        if n_pairs > 0:
            g_r = counts / n_pairs / shell_volumes * volume_sphere
        else:
            g_r = np.zeros(n_bins)
            print(f"  ⚠️  {env_name} / {target_name} 無資料，輸出全零。")

        np.savetxt(
            filename,
            np.column_stack([r_centers, g_r]),
            header="r(Angstrom)  g(r)",
            fmt="%.6f"
        )
        print(f"  ✅ 已儲存: {filename}")

print("\n=========================================")

# ==========================================
# 繪製 2x2 圖表
# ==========================================
print("🎨 正在產生 RDF 分析圖表...")

plt.rcParams['font.family']       = 'sans-serif'
plt.rcParams['axes.linewidth']    = 1.5
plt.rcParams['xtick.direction']   = 'in'
plt.rcParams['ytick.direction']   = 'in'

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    'Electrode / Functional Groups RDF under 2V\n' +
    '(Leave a 4 Å gap on the left side of the CNT)',
    fontsize=18, fontweight='bold', y=0.96
)

plot_styles = {
    'SO2_polar'     : {'color': '#d62728', 'label': 'TFSI: $-SO_2$ (Polar)'},
    'CF3_nonpolar'  : {'color': '#1f77b4', 'label': 'TFSI: $-CF_3$ (Non-polar)'},
    'Im_polar'      : {'color': '#2ca02c', 'label': 'BMIM: $-Im$ Ring (Polar)'},
    'CH3_nonpolar'  : {'color': '#ff7f0e', 'label': 'BMIM: $-CH_3$ (Non-polar)'},
    'Butyl_nonpolar': {'color': '#9467bd', 'label': 'BMIM: $-C_4H_9$ (Non-polar)'}
}

subplots_config = [
    ('CNT',      'TFSI', ['SO2_polar', 'CF3_nonpolar'],                  axes[0, 0]),
    ('CNT',      'BMIM', ['Im_polar', 'CH3_nonpolar', 'Butyl_nonpolar'], axes[0, 1]),
    ('Junction', 'TFSI', ['SO2_polar', 'CF3_nonpolar'],                  axes[1, 0]),
    ('Junction', 'BMIM', ['Im_polar', 'CH3_nonpolar', 'Butyl_nonpolar'], axes[1, 1])
]

for env_name, ion_name, target_keys, ax in subplots_config:
    has_data = False
    for target in target_keys:
        file_name = f"2V_rdf_{env_name}_{target}_V9.dat"
        if os.path.exists(file_name):
            data = np.loadtxt(file_name, comments='#')
            r, g_r = data[:, 0], data[:, 1]
            ax.plot(r, g_r, lw=2.5,
                    color=plot_styles[target]['color'],
                    label=plot_styles[target]['label'])
            has_data = True

    if has_data:
        ax.axhline(y=1.0, color='gray', linestyle=':', lw=2)
        ax.set_xlim(0, rdf_max)
        ax.set_ylim(bottom=0)
        ax.set_xlabel('Distance from Electrode, $r$ (Å)', fontsize=14)
        ax.set_ylabel('g(r)', fontsize=14)
        ax.set_title(f"{ion_name} on {env_name}", fontsize=15, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(fontsize=12, loc='lower right')
    else:
        ax.text(0.5, 0.5, "Data Not Available",
                ha='center', va='center', fontsize=14, color='red',
                transform=ax.transAxes)
        ax.set_title(f"{ion_name} on {env_name}", fontsize=15, fontweight='bold')

plt.tight_layout(rect=[0, 0.03, 1, 0.94])
output_image = '2V_rdf_elec_ions_jc7_V9.png'
plt.savefig(output_image, dpi=600, bbox_inches='tight')
print(f"\n🎉 分析與繪圖完成！圖表已儲存為：{output_image}")
