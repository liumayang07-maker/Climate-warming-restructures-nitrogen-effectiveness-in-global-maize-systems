# Climate warming restructures nitrogen effectiveness in global maize systems

本仓库提供论文 **Climate warming restructures nitrogen effectiveness in global maize systems** 对应的代码脚本和图表层级数据导出文件。

本研究整合文献来源的玉米氮肥响应配对观测数据和小时尺度 ERA5 气候数据，用于评估物候阶段相关热暴露如何改变氮肥有效性。研究将试验点年份划分为正常环境、花前热风险环境（pre-flowering risk, PFR）和灌浆期热风险环境（grain-filling risk, GFR），并进一步分析这些热风险环境如何影响产量响应、农学氮效率、区域响应格局和候选管理措施。

## 仓库内容

仓库按论文图表组织。每个图目录下均包含对应的 `README.md`，说明该图相关脚本、公开输入文件、输出文件以及与论文 Source Data 的对应关系。

| 目录 | 内容 |
| --- | --- |
| `Figure1/` | ERA5/GDD/EDH 处理流程、Figure 1 直方图输入和热暴露相关输出。 |
| `Figure2/` | Figure 2 主图 `lnRR` 和 `AEN` 混合效应模型输入及数值输出。 |
| `Figure3/` | Figure 3 区域化分析、Supplementary Figures 3-5 输出和 Supplementary Tables 1-3。 |
| `Figure4/` | Figure 4 的种植密度分析和 Supplementary Table 4。 |
| `Figure5/` | Figure 5 的氮肥管理策略分析。 |
| `Figure6/` | Figure 6 适应性管理策略的产量比较分析。 |

## 数据说明

完整的 manuscript Source Data workbook 已作为 `Source data.xlsx` 放在仓库根目录。本仓库同时提供按图拆分的 CSV 导出文件和分析脚本，便于按图核查和复现主要数值结果。

完整的项目级整理数据表 `meta_data_v2.csv` 不包含在本仓库中。仓库中只提供各图脚本运行所需的最小输入表。

原始小时尺度 ERA5 文件不存放在 GitHub 中，因为这些文件体量较大且可由脚本重新生成。使用者可根据公开的坐标表和 ERA5 下载/转换脚本，在本地配置 Copernicus Climate Data Store 凭据后重建这些气候数据。

## 分析流程

本仓库覆盖的主要分析流程包括：

- 文献来源的玉米氮肥响应配对观测整理
- ERA5 气候数据提取和 UTC 到本地小时的转换
- 生长度日（GDD）和极端度小时（EDH）计算
- 基于物候阶段的热风险分类
- `lnRR` 和农学氮效率（AEN）效应量计算
- 基于 `nlme` 的混合效应模型
- 基于 `emmeans` 的边际均值估计和组间比较
- Wilcoxon、Dunn 和 Games-Howell 等统计检验
- 图表层级 CSV 导出和脚本化数值检查

Figure 1A 的站点分布图在 ArcGIS 中完成，Figure 6C,E 的最终绘图在 Origin 2025 中完成；这些图对应的基础数值数据已在 CSV 文件中提供。

## 说明

本仓库不包含可编辑图件、稿件文档、本地凭据、虚拟环境、受版权限制的文献文件或大型生成气候文件。

文件命名遵循当前稿件中的主图、Supplementary Figure 和 Supplementary Table 编号。若部分 Source Data block label 沿用了旧命名，对应关系会在各图目录的 README 中说明。
