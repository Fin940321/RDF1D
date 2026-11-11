# 第七階段：RDF1D 整理紀錄

日期：2026-10-23。

## 範圍與搬移

專案由 `New_Analysis/Others/RDF1D/` 搬至 `New_Analysis/Code/RDF1D/`，並將原本獨立放置的 `New_Analysis/Code/rdf_2cntV9.py` 移入 `code/`。保留原 `.git`、分支、remote 與提交歷史；未 commit 或 push。

原 `examples/` 移至 `expired/examples/`。其中 2CNT 範例的原始 PDB、DCD 和四個 `.dat` 表格移至 `tab/legacy_rdf_2cnt/`；說明文件移至 `report/legacy/`；示意圖移至 `pic/illustrations/`。原 `rdf_2cnt_overlay/` 整體移至 `pic/legacy_overlay/`，保留其中空資料夾。壓縮檔與筆記本原樣保留，未刪除或解壓。

本次填寫原本空白的 README、新增本紀錄及 `file_map.tsv`，並更新工作區 `CLAUDE.md` 的目前 RDF 程式位置與 Bulk 行為說明。其餘原檔內容均保留。完整對照見 [file_map.tsv](file_map.tsv)；其中 README 的雜湊是填寫前的空檔雜湊。

## 路徑與方法限制

V9 讀取執行目錄中的三個固定 2V 檔案，程式搬移不影響這個規則。舊 2CNT 範例若以相對檔名讀取 `FV_NVT.dcd` 或 `start_drudes.pdb`，需以 `tab/legacy_rdf_2cnt/` 為工作目錄，再以完整程式路徑執行；筆記本也需設定相應工作目錄。其他 3CNT、MDTraj 範例的相依資料需逐一核對，本次未將所有舊範例改成可直接執行。

舊說明、V9 版本註解與圖標題存在歷史差異；README 依 V9 實際實作說明其條件分類及配對數正規化。本次未修正科學算法，也未重新產生正式研究結果。原 `Elec_sys/`、`Result/`、`Simulation/` 與疊圖作品集均未變更。

## 驗證方式

搬移前後使用相同 Pixi Python 與相同小型輸入執行未修改的 V9：從 `Result/BMIM_TFSI_2CNTG_100ns_CV_Done/2V/` 的軌跡中段取四幀到暫存目錄，程式分析其中後兩幀。輸入 PDB 複製至暫存區，正式輸入未寫入。

檢查 10 個數據表各為 75 × 2、數值有限、g(r) 非負，且共同距離網格為 0.1–14.9 Å、間距 0.2 Å；另比對全部表格與 PNG 的 SHA-256。這項檢查驗證搬移前後行為一致，不代表已完成全軌跡科學驗證。

驗證結果：使用原軌跡第 5000–5003 幀（零起算），實際分析第 5002–5003 幀。搬移前後 10 個表格與一張 PNG 逐位元組相同，上述數值檢查全部通過。原專案 103 個檔案中，僅原空白 README 增加內容，其餘 102 個檔案（含 Git 內部檔案）的 SHA-256、修改時間、權限及 inode 均保留；移入的 V9 內容也與原檔相同。Git HEAD、main 分支與 origin 未變更。

既有 MDAnalysis 警告包含 PDB 未知元素與 DCDReader 未來行為變更；搬移前已出現，未藉本次整理更動原始 PDB 或套件版本。

本機暫存稽核目錄：`/tmp/analysis-stage7-odeoyrn5/`，含搬移前檔案清單、驗證腳本、標準輸出、警告與搬移後檢查結果。此位置不是永久備份。

## 2026-10-29 外層整理更新

目前路徑為 `/home/atuo/Academic_Workspace/Research/Code/RDF1D/`。先前紀錄中的來源路徑保留作歷史對照。README 的目前執行範例已更新；共通手動事項見工作區根目錄 MANUAL_ADJUSTMENTS.md。
