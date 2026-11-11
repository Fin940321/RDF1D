# RDF1D：單一電壓的電極／離子基團分析

The original files of radial distribution function for 2CNT/G and 3CNT/G systems.

本專案保留 2CNT/G、3CNT/G 的原始 RDF 範例，並包含 BMIM/TFSI 與 2CNT／石墨烯電極的分析程式。

## 目前分析入口

| 程式 | 用途 |
| --- | --- |
| [`code/RDF_2CNT/rdf_cntg_ion.py`](code/RDF_2CNT/rdf_cntg_ion.py) | 電極與離子官能基的距離分布 |
| [`code/RDF_2CNT/rdf_ion_ion.py`](code/RDF_2CNT/rdf_ion_ion.py) | 離子基團之間的 RDF，包含 SO2–Im 與 CF3–CF3 |

兩支程式目前均設定 `voltage = "6V"`、`framestart = 1000`，分析影格索引 1000 至最後一幀；這個設定不會自動取軌跡後半段。切換電壓或取樣區間時，需核對各程式頂部設定與對應輸入。

輸入檔案從**執行目錄**讀取，目前為 `6V_start_drudes.pdb` 與 `6V_FV_NVT.dcd`。準備對應資料後執行，例如：

```bash
cd /path/to/your/6V_analysis
python /path/to/RDF1D/code/RDF_2CNT/rdf_cntg_ion.py
python /path/to/RDF1D/code/RDF_2CNT/rdf_ion_ion.py
```

將範例路徑換成自己的資料與專案位置。輸出資料夾建立於執行目錄中，重新執行會覆寫同名結果。套件版本列於 [`requirements.txt`](requirements.txt)；大型 DCD 受 `.gitignore` 排除，取得 GitHub 程式後仍需另備原始軌跡。

## 資料安排

- [`code/original_code/`](code/original_code/)：2CNT/G 與 3CNT/G 的原始程式和壓縮檔。
- [`expired/`](expired/)：歷史版本與筆記本。
- [`tab/legacy_rdf_2cnt/`](tab/legacy_rdf_2cnt/)：舊範例資料，使用前需核對是否符合目前程式的原子選取。
- [`pic/`](pic/)：保存的示意圖與歷史圖片；目前程式的輸出位置依程式設定。
- [`report/`](report/)：方法說明與整理紀錄。

<details>
<summary>歷史 V9 使用說明與分析設定</summary>

以下保留先前 V9 版本的說明。此次提交已改用上方兩支分析入口，未包含原路徑 `code/rdf_2cntV9.py`；以下 V9 的選取條件、輸出檔名與執行範例僅供歷史參考。目前程式的設定應以各自原始碼為準。

V9 原主程式路徑為 `code/rdf_2cntV9.py`。它分析 BMIM/TFSI 與 2CNT／石墨烯電極的距離分布；原本的 V1–V5、2CNT／3CNT 範例與筆記本保留作為歷史參考。

## 資料夾

```text
RDF1D/
├── code/rdf_2cntV9.py              # 目前主程式，內容未修改
├── expired/examples/              # 原 examples，保留歷史程式、筆記本與壓縮檔
├── pic/
│   ├── illustrations/             # 示意圖
│   └── legacy_overlay/            # 原 rdf_2cnt_overlay，保留內部結構
├── report/
│   ├── legacy/RDF_Code_Explanation.md
│   ├── organization.md            # 搬移範圍、驗證與已知限制
│   └── file_map.tsv               # 原路徑、新路徑、搬移前 SHA-256
└── tab/legacy_rdf_2cnt/            # 原範例的資料表、PDB 與 DCD
```

`tab/legacy_rdf_2cnt/` 是舊範例資料，不能直接視為 V9 所需的 2V BMIM/TFSI 輸入。DCD 仍受原有 `.gitignore` 排除，GitHub clone/pull 不保證取得本機大型軌跡。

## 執行 V9

需要 NumPy、MDAnalysis、Matplotlib、tqdm；本工作區沿用根目錄的 Pixi 環境。程式固定從**執行目錄**讀取：

- `2V_start_drudes.pdb`
- `2V_FV_NVT.dcd`
- `2V_start_nodrudes.pdb`

先準備含上述檔案的分析目錄，再執行以下命令；將範例路徑換成自己的輸入目錄。程式會在該目錄寫出結果，重跑會覆寫同名輸出。

```bash
cd /path/to/your/2V_analysis
/home/atuo/Academic_Workspace/.pixi/envs/default/bin/python /home/atuo/Academic_Workspace/Research/Code/RDF1D/code/rdf_2cntV9.py
```

這是目前本機路徑；未來外層資料夾改名時需更新。程式尚未提供電壓命令列參數，輸入、輸出名稱和圖標題皆固定為 2V。

## 目前計算設定

- 長度單位為 Å；分析軌跡後半段。
- 以第一幀定義 CNT／Junction 電極原子集合；Junction 判定距離為 7 Å。
- 目標是五種基團的指定原子，非基團質心；僅取 z < 80 Å 的目標原子。
- 接觸距離為 6 Å，先分配至 Junction，再分配至 CNT；Bulk 有分類，但不計算或輸出 RDF。
- 距離範圍 0–15 Å、75 bins，中心為 0.1–14.9 Å、間距 0.2 Å。
- 每個環境／基團以 `counts(r) / N_pairs / V_shell(r) * V_sphere(15 Å)` 正規化；`N_pairs` 為累積的截斷範圍內配對數。這不是以體相數密度正規化的標準 bulk RDF，解讀時需保留此差異。

輸出為 10 個 `2V_rdf_{CNT|Junction}_{group}_V9.dat` 與一張 `2V_rdf_elec_ions_jc7_V9.png`，都位於執行目錄。`pic/` 與 `tab/` 是資料整理位置，現有程式不會自動寫入這些子目錄。

舊說明 [`RDF_Code_Explanation.md`](report/legacy/RDF_Code_Explanation.md) 原文保留，部分方法描述與 V9 不同。程式圖標題中的「4 Å gap」及 V7/V8 註解也是既有文字，不能據此判定 V9 實際採用了相應設定；本次只整理檔案，未改動計算或標題。
</details>
