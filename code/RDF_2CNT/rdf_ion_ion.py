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
r_electrode_contact = 6.0   # 離子族群分類閾值 (Å)，同論文 Fig. 9 的 6 Å (Å)
z_max    = 65               # 參考原子 z 截斷（工作電極側）
rdf_max  = 15
n_bins   = 75
dr       = float(rdf_max) / float(n_bins)
framestart = 1000  # 後半段軌跡，避開平衡化
frameend   = len(u.trajectory)

# --- 輸出資料夾 ---
# 將本程式產生的數據與圖表集中存放；資料夾不存在時自動建立。
output_dir = f"{voltage}_rdf_ion_ion"
os.makedirs(output_dir, exist_ok=True)

# --- 設定分析目標 ---
resname1 = 'Tf2'
SO2_atoms = 'Stf Stf1'
CF3_atoms = 'Ctf Ctf1'
resname2 = 'BMI'
Im_atoms = 'C1 C2 C21'
Butyl_atoms = 'C4 C5 C51 C6'


# --- 設定電極分類 ---
classifier=[]
classifier.append( "segid A" )
classifier.append( "segid B" )
classifier.append( "segid C" )

# make static atom groups for electrode.
u.trajectory[0]
electrode_groups=[]
for electrode in classifier:
    electrode_groups.append( u.select_atoms(electrode) )

# --- 設定 RDF 計算區域 ---
environments_counts = 2 # 0: 純碳管區, 1: 交界區
# SO2-Im RDF 獨立累積陣列
rdf_count = [ [ 0 for _ in range(n_bins) ] for _ in range(environments_counts) ]
rdf_N_N = [ 0 for _ in range(environments_counts) ] # keeps track of number of particles histogramed for normalization

# CF3-CF3 RDF 獨立累積陣列；同一顆 TFSI 的配對會在距離計算後排除
rdf_count_cf3cf3 = np.zeros((environments_counts, n_bins), dtype=np.int64)
rdf_N_N_cf3cf3 = np.zeros(environments_counts, dtype=np.int64)

# --- 開始軌跡迴圈 ---
print(f"📊 開始分析軌跡: frame {framestart} → {frameend - 1}（共 {frameend - framestart} 幀）\n")

for t0 in tqdm(range(framestart, frameend), desc="Processing frames", unit="frame", ncols=80):
    u.trajectory[t0]
    
    # 1. 直接選取 TFSI 的 SO2 氧原子群組(參考)
    group1 = u.select_atoms("name %s and resname %s and prop z < %s" % (SO2_atoms, resname1, z_max))
    # 2. 直接選取 BMIM 的 Im ring碳原子群組(目標)
    group2 = u.select_atoms("name %s and resname %s and prop z < %s" % (Im_atoms, resname2, z_max + rdf_max))
    # 3. 直接選取 TFSI 的 CF3 碳原子群組(參考)
    group3 = u.select_atoms("name %s and resname %s and prop z < %s" % (CF3_atoms, resname1, z_max))
    # 4. 直接選取 TFSI 的 CF3 碳原子群組(目標)
    group4 = u.select_atoms("name %s and resname %s and prop z < %s" % (CF3_atoms, resname1, z_max + rdf_max))

    if any(len(group) == 0 for group in (group1, group2, group3, group4)):
        continue

    # 3. 計算與三個電極的最短距離 (0: 石墨烯, 1: CNT1, 2: CNT2)
    contact_dist_electrodes = []
    close_contact_electrodes = []
    contact_dist_electrodes_cf3 = []
    close_contact_electrodes_cf3 = []

    # 參考原子: TFSI 的 SO2 氧原子群組
    for electrode in electrode_groups:
        dist = distances.distance_array(group1.positions, electrode.positions, box=u.dimensions, backend='OpenMP')
        contact_dist = np.amin(dist, axis=1)
        contact_dist_electrodes.append(contact_dist)

    # 參考原子: TFSI 的 CF3 碳原子群組
    for electrode in electrode_groups:
        dist_cf3 = distances.distance_array(group3.positions, electrode.positions, box=u.dimensions, backend='OpenMP')
        contact_dist_cf3 = np.amin(dist_cf3, axis=1)
        contact_dist_electrodes_cf3.append(contact_dist_cf3)
      
    # 4. 判斷離子是否接觸電極 (距離 < 6.0 Å)
    # 參考原子: TFSI 的 SO2 氧原子群組
    # 石墨烯板接觸
    close_contact_electrodes.append(abs(np.heaviside(contact_dist_electrodes[0] - r_electrode_contact, 1) - 1))
    # 合併兩支 CNT 的接觸
    temp1 = abs(np.heaviside(contact_dist_electrodes[1] - r_electrode_contact, 1) - 1)
    temp2 = abs(np.heaviside(contact_dist_electrodes[2] - r_electrode_contact, 1) - 1)
    temp1 += temp2
    close_contact_electrodes.append(np.ma.masked_not_equal(temp1, 0).mask.astype(int))

    # 參考原子: TFSI 的 CF3 碳原子群組
    # 石墨烯板接觸
    close_contact_electrodes_cf3.append(abs(np.heaviside(contact_dist_electrodes_cf3[0] - r_electrode_contact, 1) - 1))
    # 合併兩支 CNT 的接觸
    temp1_cf3 = abs(np.heaviside(contact_dist_electrodes_cf3[1] - r_electrode_contact, 1) - 1)
    temp2_cf3 = abs(np.heaviside(contact_dist_electrodes_cf3[2] - r_electrode_contact, 1) - 1)
    temp1_cf3 += temp2_cf3
    close_contact_electrodes_cf3.append(np.ma.masked_not_equal(temp1_cf3, 0).mask.astype(int))

    # 5. 計算環境標籤 (1: 石墨烯, 2: CNT, 3: CNT/G)
    # 建立 SO2 reference 的環境標籤
    environment = np.zeros(len(close_contact_electrodes[0]))
    for i in range(len(close_contact_electrodes)):
        environment += close_contact_electrodes[i] * (2**i)

    # 建立 CF3 reference 的環境標籤
    environment_cf3 = np.zeros(len(close_contact_electrodes_cf3[0]))
    for i in range(len(close_contact_electrodes_cf3)):
        environment_cf3 += close_contact_electrodes_cf3[i] * (2**i)

    # 6. 分類離子到對應環境
    # TFSI 的 SO2 氧原子參考群組
    environment_atom = [[] for _ in range(2)]
    for i in range(len(group1.indices)):
        if environment[i] == 2.0:
            environment_atom[0].append(group1[i])  # 存入環境 0: 純 CNT 區
        elif environment[i] == 3.0:
            environment_atom[1].append(group1[i])  # 存入環境 1: 交界區 (CNT/G)
        elif environment[i] == 1.0 or environment[i] == 0.0:
            pass  # 忽略只接觸石墨烯或都沒接觸的離子
    
    # TFSI 的 CF3 碳原子參考群組
    environment_atom_cf3 = [[] for _ in range(environments_counts)]
    for i in range(len(group3.indices)):
        if environment_cf3[i] == 2.0:
            # 純 CNT：接觸 CNT，但沒有接觸 graphene
            environment_atom_cf3[0].append(group3[i])
        elif environment_cf3[i] == 3.0:
            # CNT/G：同時接觸 CNT 與 graphene
            environment_atom_cf3[1].append(group3[i])
        elif environment_cf3[i] == 1.0 or environment_cf3[i] == 0.0:
            # graphene-only 或沒有接觸電極，不納入這次分析
            pass

    # 7. 累積距離直方圖
    # (1)計算SO2-Im RDF
    for i in range(len(environment_atom)):
        if environment_atom[i]:
            pairs, dist = distances.capped_distance(
                AtomGroup(environment_atom[i]).positions, group2.positions, 
                rdf_max, box=u.dimensions
            )
            count = np.histogram(dist, bins=n_bins, range=(0.0, rdf_max))[0]
            rdf_count[i] += count
            rdf_N_N[i] += len(pairs)

    # (2)計算CF3-CF3 RDF
    for i in range(len(environment_atom_cf3)):
        if environment_atom_cf3[i]:
            # 將目前環境的 CF3 reference 轉成 AtomGroup
            cf3_refs_in_environment = AtomGroup(environment_atom_cf3[i])

            # 計算 CF3 reference 到所有 CF3 target 的距離
            pairs_cf3, dist_cf3cf3 = distances.capped_distance(
                cf3_refs_in_environment.positions, group4.positions,
                rdf_max, box=u.dimensions
            )
            # 15 Å 內沒有任何配對
            if len(pairs_cf3) == 0:
                continue

            # 取得每一組 pair 兩端所屬的 TFSI residue index
            reference_resindices = (cf3_refs_in_environment.resindices[pairs_cf3[:, 0]])
            target_resindices = (group4.resindices[pairs_cf3[:, 1]])

            # # 排除同一顆 TFSI 分子內的 CF3-CF3 配對
            # keep_inter_residue = (reference_resindices != target_resindices)
            # dist_cf3cf3 = dist_cf3cf3[keep_inter_residue]
            # if len(dist_cf3cf3) == 0:
            #     continue

            # 將不同 TFSI 之間的 CF3-CF3 距離累積到 histogram
            count_cf3cf3 = np.histogram(dist_cf3cf3, bins=n_bins, range=(0.0, rdf_max))[0]

            rdf_count_cf3cf3[i] += count_cf3cf3
            rdf_N_N_cf3cf3[i] += len(dist_cf3cf3)

# --- 軌跡結束，進行正規化 ---
print("\n⏳ 正在進行正規化與儲存數據...")
edges = np.linspace(0.0, rdf_max, n_bins + 1)
V = (4 / 3) * np.pi * (np.power(edges[1:], 3) - np.power(edges[:-1], 3))
volume_sphere = (4 / 3) * np.pi * rdf_max**3

# 1. 執行 SO₂–Im 正規化
for i in range(environments_counts):
    if rdf_N_N[i] > 0:
        rdf_count[i] = (rdf_count[i] / rdf_N_N[i] / V * volume_sphere)

# 2-1. 建立 CF3-CF3 的浮點正規化陣列
rdf_cf3cf3 = np.zeros_like(rdf_count_cf3cf3, dtype=float)

# 2-2. 執行 CF3-CF3 正規化
for i in range(environments_counts):
    histogram_pair_total = np.sum(rdf_count_cf3cf3[i])
    if histogram_pair_total != rdf_N_N_cf3cf3[i]:
        raise RuntimeError(
            "CF3-CF3 histogram 與 pair total 不一致："
            f"environment={i}, "
            f"histogram={histogram_pair_total}, "
            f"pairs={rdf_N_N_cf3cf3[i]}"
        )
    if rdf_N_N_cf3cf3[i] > 0:
        rdf_cf3cf3[i] = (
            rdf_count_cf3cf3[i].astype(float)
            / rdf_N_N_cf3cf3[i]
            / V
            * volume_sphere
        )

# --- 寫出結果數據檔案 ---
# SO₂–Im 輸出檔案:
# 1. 寫出純碳管區 (環境 0)
output_so2im_cnt = os.path.join(output_dir, "rdf_2cnt_so2im_cnt.dat")
with open(output_so2im_cnt, "w") as rdf_ion_ion_cnt:
    for j in range(len(rdf_count[0])):
        print('{0:5.8f}  {1:5.8f}'.format(j*dr, rdf_count[0][j]), file=rdf_ion_ion_cnt)
print(f"✅ 已儲存：{output_so2im_cnt}")

# 2. 寫出交界區 (環境 1)
output_so2im_cntg = os.path.join(output_dir, "rdf_2cnt_so2im_cntg.dat")
with open(output_so2im_cntg, "w") as rdf_ion_ion_cntg:
    for j in range(len(rdf_count[1])):
        print('{0:5.8f}  {1:5.8f}'.format(j*dr, rdf_count[1][j]), file=rdf_ion_ion_cntg)
print(f"✅ 已儲存：{output_so2im_cntg}")

# CF3-CF3 輸出檔案:
# 3. 寫出純碳管區 (環境 0)
output_cf3cf3_cnt = os.path.join(output_dir, "rdf_2cnt_cf3cf3_cnt.dat")
with open(output_cf3cf3_cnt, "w") as rdf_ion_ion_cnt:
    for j in range(len(rdf_cf3cf3[0])):
        print('{0:5.8f}  {1:5.8f}'.format(j*dr, rdf_cf3cf3[0][j]), file=rdf_ion_ion_cnt)
print(f"✅ 已儲存：{output_cf3cf3_cnt}")

# 4. 寫出交界區 (環境 1)
output_cf3cf3_cntg = os.path.join(output_dir, "rdf_2cnt_cf3cf3_cntg.dat")
with open(output_cf3cf3_cntg, "w") as rdf_ion_ion_cntg:
    for j in range(len(rdf_cf3cf3[1])):
        print('{0:5.8f}  {1:5.8f}'.format(j*dr, rdf_cf3cf3[1][j]), file=rdf_ion_ion_cntg)
print(f"✅ 已儲存：{output_cf3cf3_cntg}")

# --- 繪製並儲存分析圖表 ---
print("\n🎨 正在繪製並儲存 RDF 分析圖表...")
r_centers = (edges[:-1] + edges[1:]) / 2

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 1. 繪製 SO₂–Im RDF 圖表
plt.figure(figsize=(10, 6))
# 畫第一條線：純碳管區 (環境 0)
plt.plot(r_centers, rdf_count[0], lw=2.5, color='#1f77b4', label='CNT: [TFSI-] -SO2 ⟷ [BMIM+] -Im')
# 畫第二條線：交界區 (環境 1)
plt.plot(r_centers, rdf_count[1], lw=2.5, color='#ff7f0e', label='CNT/G: [TFSI-] -SO2 ⟷ [BMIM+] -Im')

plt.xlim(0, rdf_max)
plt.ylim(bottom=0)
plt.xlabel('Distance (Å)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.title(f'[TFSI-] - [BMIM+] RDF ({voltage})', fontsize=16, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=12, loc='upper right')

output_image = os.path.join(output_dir, "rdf_so2im_combined.png")
plt.savefig(output_image, dpi=600, bbox_inches='tight')
print(f"✅ 圖表已成功儲存為：{output_image}")
plt.close()

# 2. 繪製 CF3-CF3 RDF 圖表
plt.figure(figsize=(10, 6))
# 畫第一條線：純碳管區 (環境 0)
plt.plot(r_centers, rdf_cf3cf3[0], lw=2.5, color='#1f77b4', label='CNT: [TFSI-] -CF3 ↔ [TFSI-] -CF3')
# 畫第二條線：交界區 (環境 1)
plt.plot(r_centers, rdf_cf3cf3[1], lw=2.5, color='#ff7f0e', label='CNT/G: [TFSI-] -CF3 ↔ [TFSI-] -CF3')

plt.xlim(0, rdf_max)
plt.ylim(bottom=0)
plt.xlabel('Distance (Å)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.title(f'[TFSI-] CF$_3$–CF$_3$ RDF ({voltage})', fontsize=16, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=12, loc='upper right')
output_image_cf3 = os.path.join(output_dir, "rdf_cf3cf3_combined.png")

plt.savefig(output_image_cf3, dpi=600, bbox_inches='tight')
print(f"✅ 圖表已成功儲存為：{output_image_cf3}")
print("\n🎉 全部分析與作圖完成！")
plt.close()
