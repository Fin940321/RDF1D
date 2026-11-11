from MDAnalysis import *
import MDAnalysis.analysis.distances as distances
import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm import tqdm

# --- 引入分析資料 ---
voltage = "6V"
topology = f"{voltage}_start_drudes.pdb"
trajectory = f"{voltage}_FV_NVT.dcd"
u = Universe(topology, trajectory)

# --- 參數設定 ---
electrode_region_cutoff = 7.0  # CNT/CNTG 電極區域切割距離 (Å)
blocked_cutoff = 5.0           # 排除遭 CNT 遮蔽的 graphene 碳 (Å)
rdf_max  = 15
n_bins   = 75
dr       = float(rdf_max) / float(n_bins)
framestart = 1000  # 後半段軌跡，避開平衡化
frameend   = len(u.trajectory)

# --- 設定分析目標 ---
all_targets = {
    # TFSI 的極性 SO2：使用四個氧原子
    'SO2_polar': (
        "resname Tf2 and "
        "name Stf Stf1"
    ),

    # TFSI 的非極性 CF3：使用兩個 CF3 碳原子
    'CF3_nonpolar': (
        "resname Tf2 and "
        "name Ctf Ctf1"
    ),

    # BMIM 的極性 imidazolium ring：使用三個環碳
    'Im_polar': (
        "resname BMI and "
        "name C1 C2 C21"
    ),

    # BMIM 的甲基
    'CH3_nonpolar': (
        "resname BMI and "
        "name C3"
    ),

    # BMIM 的 butyl chain
    'Butyl_nonpolar': (
        "resname BMI and "
        "name C4 C5 C51 C6"
    ),
}

# --- 建立電極 reference 區域 ---
# 電極是固定的，因此在第 0 幀建立靜態 AtomGroup
u.trajectory[0]

# 完整平坦 graphene 與兩根 CNT
grp_graphene = u.select_atoms("segid A")
cnt_all = u.select_atoms("segid B or segid C")

# CNT/G 中的 CNT 部分：
# 選出距離 graphene 小於 6 Å 的 CNT 碳原子
cnt_in_cntg = u.select_atoms(
    f"(segid B or segid C) and "
    f"around {electrode_region_cutoff} group grp",
    grp=grp_graphene
)

# CNT/G 附近所有可能的 graphene 碳：
# 先選出距離 CNT 小於 6 Å 的 graphene 碳原子
grp_in_cntg_all = u.select_atoms(
    f"segid A and "
    f"around {electrode_region_cutoff} group cnt",
    cnt=cnt_all
)

# 被 CNT 遮蔽的 graphene 碳：
# 距離 CNT 小於 5 Å，視為離子無法接近
grp_blocked = u.select_atoms(
    f"segid A and "
    f"around {blocked_cutoff} group cnt",
    cnt=cnt_all
)

# CNT/G 中真正可作為 RDF 中心的 graphene：
# 距 CNT 小於 6 Å，但排除小於 5 Å 的 blocked 原子
grp_in_cntg = grp_in_cntg_all - grp_blocked

# CNT/G reference：
# CNT 下緣 + 未被遮蔽的局部 graphene
center_cntg = cnt_in_cntg + grp_in_cntg

# 純 CNT reference：
# 排除靠近 graphene 6 Å 內的 CNT 下緣
center_cnt = cnt_all - cnt_in_cntg

### 原子數量檢查
print("✅ 電極 reference 區域建立完成")
print(f"   完整 graphene：{len(grp_graphene)}")
print(f"   完整 CNT：{len(cnt_all)}")
print(f"   CNT/G 中的 CNT：{len(cnt_in_cntg)}")
print(f"   CNT/G 附近所有 graphene：{len(grp_in_cntg_all)}")
print(f"   blocked graphene：{len(grp_blocked)}")
print(f"   CNT/G 中可接近的 graphene：{len(grp_in_cntg)}")
print(f"   CNT/G reference 總數：{len(center_cntg)}")
print(f"   純 CNT reference 總數：{len(center_cnt)}")

if len(center_cntg) == 0:
    raise ValueError("CNT/G reference 為空，請檢查 6/5 Å 選區。")

if len(center_cnt) == 0:
    raise ValueError("純 CNT reference 為空，請檢查 6 Å 選區。")

if len(grp_in_cntg) == 0:
    print("⚠️ 5～6 Å 之間沒有可接近的 graphene 碳，CNT/G 將只由 CNT 碳構成。")

# --- 建立離子官能基 target AtomGroup ---
# 電極最高 z 加上 RDF 最大距離。
# 超過這個 z 的離子不可能與目前工作電極形成 15 Å 內配對。
electrode_z_max = max(
    grp_graphene.positions[:, 2].max(),
    cnt_all.positions[:, 2].max()
)

target_z_cutoff = electrode_z_max + rdf_max

print(
    f"✅ 電極最高 z：{electrode_z_max:.3f} Å，"
    f"target z 截斷：{target_z_cutoff:.3f} Å"
)

# 這裡先建立完整官能基 AtomGroup。
# z 截斷會在每一幀依照當下座標套用。
atom_groups = {
    target_name: u.select_atoms(selection)
    for target_name, selection in all_targets.items()
}

print("✅ 離子官能基 target 建立完成")
for target_name, target_group in atom_groups.items():
    print(f"{target_name}：{len(target_group)} 個原子")
    if len(target_group) == 0:
        raise ValueError(
            f"找不到 target 原子：{target_name}，"
            f"請檢查 resname 與 atom name。"
        )

# --- 設定兩種電極 reference 環境 ---
electrode_centers = {
    'CNT': center_cnt,
    'CNTG': center_cntg,
}

# --- 建立 RDF 原始累積陣列 ---
# 結構：
# rdf_counts[電極環境][官能基][bin]
rdf_counts = {
    env_name: {
        target_name: np.zeros(n_bins, dtype=np.int64)
        for target_name in all_targets
    }
    for env_name in electrode_centers
}

# 累積每一種「電極環境－官能基」在 rdf_max 內的總 pair 數
rdf_pair_totals = {
    env_name: {
        target_name: 0
        for target_name in all_targets
    }
    for env_name in electrode_centers
}

# --- 開始軌跡迴圈 ---
print(
    f"📊 開始分析軌跡：frame {framestart} → {frameend - 1}"
    f"（共 {frameend - framestart} 幀）\n"
)

n_frames_processed = 0

for ts in tqdm(
    u.trajectory[framestart:frameend],
    desc="Processing frames",
    unit="frame",
    ncols=80
):
    n_frames_processed += 1
    # 依序處理五種離子官能基 target
    for target_name, target_group_full in atom_groups.items():
        # 每一幀依當下座標重新套用 z 截斷
        target_z = target_group_full.positions[:, 2]
        target_group = target_group_full[
            target_z < target_z_cutoff
        ]
        # 目前幀在工作電極側沒有這類 target
        if len(target_group) == 0:
            continue
                   
        # 分別計算純 CNT 與 CNT/G 電極碳到 target 的距離
        for env_name, electrode_group in electrode_centers.items():
            pairs, dist = distances.capped_distance(
                electrode_group.positions,
                target_group.positions,
                max_cutoff=rdf_max,
                box=u.dimensions
            )

            # 此環境與 target 在 15 Å 內沒有任何配對
            if len(dist) == 0:
                continue
            
            counts, _ = np.histogram(
                dist,
                bins=n_bins,
                range=(0.0, rdf_max)
            )
            # 單次距離配對數應與 histogram 總數一致
            if np.sum(counts) != len(dist):
                raise RuntimeError(
                    "Histogram 與距離 pair 數不一致："
                    f"frame={ts.frame}, "
                    f"environment={env_name}, "
                    f"target={target_name}, "
                    f"histogram={np.sum(counts)}, "
                    f"pairs={len(dist)}"
                )

            rdf_counts[env_name][target_name] += counts
            rdf_pair_totals[env_name][target_name] += len(dist)
print(
    f"\n✅ 軌跡距離累積完成，"
    f"實際處理 {n_frames_processed} 幀。"
)

# --- 軌跡結束，進行正規化 ---
print("\n⏳ 正在進行 RDF 正規化...")

# Histogram 邊界
edges = np.linspace(
    0.0,
    rdf_max,
    n_bins + 1
)

# .dat 第一欄預計儲存 histogram 左邊界
r_left = edges[:-1]

# 畫圖時使用 bin center
r_centers = (
    edges[:-1] + edges[1:]
) / 2.0

# 每個球殼的體積
shell_volumes = (
    (4.0 / 3.0)
    * np.pi
    * (edges[1:]**3 - edges[:-1]**3)
)

# 0～rdf_max 的完整球體積
volume_sphere = (
    (4.0 / 3.0)
    * np.pi
    * rdf_max**3
)

# 建立獨立的浮點數正規化結果
rdf_normalized = {
    env_name: {
        target_name: np.zeros(n_bins, dtype=float)
        for target_name in all_targets
    }
    for env_name in electrode_centers
}

for env_name in electrode_centers:
    for target_name in all_targets:

        counts = rdf_counts[env_name][target_name]
        pair_total = rdf_pair_totals[env_name][target_name]
        histogram_total = int(np.sum(counts))

        # 確認所有 capped_distance pair 都有進入 histogram
        if histogram_total != pair_total:
            raise RuntimeError(
                "Histogram 與累積 pair total 不一致："
                f"environment={env_name}, "
                f"target={target_name}, "
                f"histogram={histogram_total}, "
                f"pairs={pair_total}"
            )

        if pair_total == 0:
            print(
                f"⚠️ {env_name} / {target_name} 沒有有效 pair，"
                "RDF 保持全零。"
            )
            continue

        g_r = (
            counts.astype(float)
            / pair_total
            / shell_volumes
            * volume_sphere
        )

        if not np.all(np.isfinite(g_r)):
            raise RuntimeError(
                "RDF 正規化後出現 NaN 或 Inf："
                f"environment={env_name}, "
                f"target={target_name}"
            )

        rdf_normalized[env_name][target_name] = g_r

        print(
            f"✅ {env_name} / {target_name}："
            f"{pair_total} pairs"
        )

# --- 寫出電極－離子 RDF 結果 ---
output_dir = (
    f"{voltage}_rdf_cntg_ion_"
    f"cut{electrode_region_cutoff:g}_"
    f"block{blocked_cutoff:g}"
)

os.makedirs(
    output_dir,
    exist_ok=True
)

print(f"\n📁 正在寫入 RDF 資料：{output_dir}/")

for env_name in electrode_centers:
    for target_name in all_targets:
        output_filename = (
            f"{voltage}_rdf_"
            f"{env_name}_"
            f"{target_name}_"
            f"cut{electrode_region_cutoff:g}_"
            f"block{blocked_cutoff:g}.dat"
        )

        output_path = os.path.join(
            output_dir,
            output_filename
        )

        output_data = np.column_stack([
            r_left,
            rdf_normalized[env_name][target_name]
        ])

        pair_total = rdf_pair_totals[env_name][target_name]
        header = (
            "r_left(Angstrom)  g(r)\n"
            f"voltage={voltage}, "
            f"environment={env_name}, "
            f"target={target_name}, "
            f"electrode_region_cutoff={electrode_region_cutoff:g} Angstrom, "
            f"blocked_cutoff={blocked_cutoff:g} Angstrom, "
            f"rdf_max={rdf_max:g} Angstrom, "
            f"n_bins={n_bins}, "
            f"frames={n_frames_processed}, "
            f"pair_total={pair_total}"
        )

        np.savetxt(
            output_path,
            output_data,
            fmt="%.8f",
            header=header
        )

        print(
            f"✅ {env_name} / {target_name}："
            f"{output_path}"
        )

# --- 繪製並儲存分析圖表 ---
print("\n🎨 正在繪製電極－離子 RDF 摘要圖...")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True

plot_styles = {
    'SO2_polar': {
        'color': '#d62728',
        'label': 'TFSI $(-SO_2)$'
    },
    'CF3_nonpolar': {
        'color': '#1f77b4',
        'label': 'TFSI $(-CF_3)$'
    },
    'Im_polar': {
        'color': '#2ca02c',
        'label': 'BMIM $(-Im)$ Ring'
    },
    'CH3_nonpolar': {
        'color': '#ff7f0e',
        'label': 'BMIM $(-CH_3)$'
    },
    'Butyl_nonpolar': {
        'color': '#9467bd',
        'label': 'BMIM $(-C_4H_9)$'
    },
}

env_display = {
    'CNT': 'CNT',
    'CNTG': 'CNT/G',
}

fig, axes = plt.subplots(
    2,
    2,
    figsize=(14, 10)
)

fig.suptitle(
    f'Electrode–Ion Functional-Group RDF ({voltage})\n'
    f'(region cutoff = {electrode_region_cutoff:g} Å; '
    f'blocked graphene < {blocked_cutoff:g} Å excluded)',
    fontsize=17,
    fontweight='bold',
    y=0.98
)

subplots_config = [
    (
        'CNT',
        'TFSI',
        ['SO2_polar', 'CF3_nonpolar'],
        axes[0, 0]
    ),
    (
        'CNT',
        'BMIM',
        ['Im_polar', 'CH3_nonpolar', 'Butyl_nonpolar'],
        axes[0, 1]
    ),
    (
        'CNTG',
        'TFSI',
        ['SO2_polar', 'CF3_nonpolar'],
        axes[1, 0]
    ),
    (
        'CNTG',
        'BMIM',
        ['Im_polar', 'CH3_nonpolar', 'Butyl_nonpolar'],
        axes[1, 1]
    ),
]

for env_name, ion_name, target_names, ax in subplots_config:
    for target_name in target_names:
        ax.plot(
            r_centers,
            rdf_normalized[env_name][target_name],
            lw=2.5,
            color=plot_styles[target_name]['color'],
            label=plot_styles[target_name]['label']
        )

    ax.set_xlim(
        0.0,
        rdf_max
    )

    ax.set_ylim(
        bottom=0.0
    )

    ax.set_xlabel(
        'Distance from Electrode, $r$ (Å)',
        fontsize=13
    )

    ax.set_ylabel(
        '$g(r)$',
        fontsize=13
    )

    ax.set_title(
        f"{ion_name} on {env_display[env_name]}",
        fontsize=14,
        fontweight='bold',
        pad=12
    )

    ax.grid(
        True,
        linestyle='--',
        alpha=0.3
    )

    ax.legend(
        fontsize=11,
        loc='upper right',
        frameon=True
    )

fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])

output_image = os.path.join(
    output_dir,
    (
        f"{voltage}_rdf_electrode_ion_"
        f"cut{electrode_region_cutoff:g}_"
        f"block{blocked_cutoff:g}.png"
    )
)

fig.savefig(
    output_image,
    dpi=600,
    bbox_inches='tight'
)

plt.close(fig)

print(f"✅ RDF 摘要圖已儲存：{output_image}")
print("\n🎉 電極－離子 RDF 分析全部完成！")