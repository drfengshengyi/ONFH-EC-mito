# Figma 拼图与排版说明书

> 项目：ONFH 内皮选择性自噬受体背景图形摘要
> 版本：Genes receptor-context revision / 2026-09-03
> 本文件只用于 Figma 拼装、准确标注和投稿输出。生物学插画素材请使用独立文件 `nanobanana2_assets_receptor_context.md` 生成。

## 1. 图形契约

### 一句话结论

> Participant-aware reanalysis of public ONFH data prioritizes inflammatory programs and a heterogeneous endothelial selective-autophagy receptor context, with SQSTM1 as the primary computational case, CALCOCO2/NDP52 as a secondary candidate, and OPTN as a mechanistic-context control; experimental validation remains required.

### 视觉目的

读者应在 5 秒内理解三件事：

1. 研究对象是股骨头骨微血管内皮细胞。
2. 最稳定的是参与者层级的炎症/干扰素程序；受体证据并不一致。
3. SQSTM1、CALCOCO2/NDP52 和 OPTN 是三个不同证据等级，不是已验证的因果链。

### 必须避免

- 不能继续使用“SQSTM1 单基因机制”为视觉中心。
- 不将三个受体用实线串成通路。
- 不把虚拟敲除画成真实实验。
- 不把单张 OA 空间切片称为 ONFH 验证。
- 不把 AUC 0.870 称为诊断性能验证。
- 不使用数据库官方 logo；使用普通文字胶囊即可。
- 不把多个主文统计图缩小后拼在一起。

## 2. 整体版式

### 画布

- Frame：1920 × 1080 px，16:9。
- 背景：`#FBFAF7`。
- 四周安全边距：72 px。
- 12 列网格，gutter 24 px。
- 不写“Graphical Abstract”。

### 阅读路径

```text
公共数据 / 参与者层级分析
          ↓
股骨头 ROI → 微血管环境放大圈 → 骨微血管内皮细胞
          ↓
炎症/干扰素程序 + 异质性受体证据
          ↓
计算扩展的边界 → 下一步平行功能验证
```

### 版式线框

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ PARTICIPANT-AWARE ONFH REANALYSIS                               [HYPOTHESIS] │
│ Inflammatory programs and heterogeneous receptor context                    │
│ 3 cohorts · 18 participants · 13,426 endothelial cells                     │
├───────────────┬───────────────────────────────────┬──────────────────────────┤
│ DATA / METHOD │            BIOLOGICAL HERO        │ EVIDENCE HIERARCHY       │
│ dataset chips │ femoral head → niche → EC close-up│ ● SQSTM1 primary case    │
│ participant   │                                   │ ◆ CALCOCO2 secondary     │
│ level filter  │ inflammatory/interferon programs  │ □ OPTN context control   │
│ OA spatial:   │ heterogeneous receptor context    │                          │
│ context only  │                                   │ no causal mechanism      │
├───────────────┴───────────────────────────────────┴──────────────────────────┤
│ vKO: no replicated nuclear FDR hit | serum AUC 0.870, internal only         │
│ NEXT: matched tissue → parallel loss/rescue → flux/function → angiogenesis  │
└──────────────────────────────────────────────────────────────────────────────┘
```

采用 schematic-led asymmetric composition。中央生物主视觉最大，右侧证据等级次之，统计扩展只占底部窄条。

## 3. 区域坐标建议

### 顶部标题区

- y = 52–164 px。
- 标题 x = 72–1400 px。
- `HYPOTHESIS-GENERATING` 胶囊 x = 1570–1848 px。
- 标题与副标题左对齐，不居中。

### 主体区

主体 y = 184–814 px。

1. 左侧输入与定位：x = 72–470 px。
2. 中央主视觉：x = 440–1280 px。
3. 右侧证据等级：x = 1300–1848 px。

区域可以有极淡色块，但不要每一项都画成厚边框卡片。

### 底部条带

- y = 842–1016 px。
- 左 56%：虚拟敲除边界。
- 中 20%：血清内部评价。
- 右 24%：下一步实验。

## 4. 多尺度主视觉拼装

### 股骨头

- 建议高度 390–430 px。
- 放置在 x = 96–420 px，y = 318–760 px。
- 不要添加血管树或坏死区域。
- ROI 放在股骨头—股骨颈上方区域，为 30–36 px 实心浅青点加 2 px 深青描边。
- ROI 外圈使用 56–66 px、2 px 青绿色虚线圆。

### 微血管环境放大圈

- 直径 235–260 px。
- 放置在 x = 330–590 px，y = 245–505 px。
- 圆外加 4–6 px `#2E918A` 描边。
- 不要压住标题、关节软骨或右侧内皮细胞。

### 股骨头 ROI 到微血管放大圈

- 必须使用两条浅灰蓝虚线，形成真正的放大锥形关系。
- 颜色 `#AEBCCA`，2 px，dash 8/8。
- 两条线从 ROI 外圈上下两个切点连接到放大圈左侧两个切点。
- 不能只画一条线。
- 不能有残留的第三条线。

### 内皮细胞放大图

- 宽度 650–760 px。
- 放置在 x = 550–1260 px，y = 325–700 px。
- 保持透明背景，不叠加白色矩形。
- 细胞右侧可留出受体证据标记的视觉方向，但受体卡放在独立右栏。

### 微血管放大圈到内皮细胞

- 使用两条浅灰蓝虚线。
- 从放大圈右侧上下切点连接到细胞左侧上下边缘。
- 股骨头 ROI 不得直接连接内皮细胞。
- 所有虚线必须置于生物素材下层，不能穿过红细胞、细胞核或文字。

## 5. 文案与层级

所有英文必须在 Figma 中输入，不使用 Nano Banana 生成文字。

### 标题

```text
Participant-aware ONFH reanalysis
```

### 主结论副标题

```text
Inflammatory programs and a heterogeneous endothelial selective-autophagy receptor context
```

### 研究规模

```text
3 scRNA-seq cohorts · 18 reported participants · 13,426 endothelial cells
```

### 状态胶囊

```text
HYPOTHESIS-GENERATING
```

### 左侧数据层

```text
PUBLIC FEMORAL-HEAD DATA
Participant-aware inference
```

数据库胶囊：

```text
SRP361778
GSE169396
GSE290411
```

空间背景小标签：

```text
GSE284089 · 1 OA section
Anatomical context only
```

### 中央主要结果

```text
INFLAMMATORY / INTERFERON PROGRAMS
Directionally stable under participant omission
```

```text
HETEROGENEOUS RECEPTOR EVIDENCE
```

### 右侧三个受体等级

```text
SQSTM1
Primary computational case
```

```text
CALCOCO2 / NDP52
Secondary transcript candidate
```

```text
OPTN
Mechanistic-context control
Not prioritized by current data
```

### 计算扩展

```text
CALIBRATED COMPUTATIONAL EXTENSIONS
```

```text
SQSTM1 in silico perturbation
No replicated nuclear-gene FDR hit
Not exceptional versus matched comparators
```

```text
Serum nested cross-validation
AUC 0.870 · Internal evaluation only
```

### 证据边界

```text
RNA abundance ≠ protein activity or autophagic flux
No receptor-specific causal mechanism established
```

### 下一步实验

```text
NEXT
Matched tissue → parallel SQSTM1/CALCOCO2 loss-rescue ± OPTN comparator
→ mitochondrial and autophagic-flux assays → endothelial survival and angiogenesis
```

### 空间不足时的短版

优先使用：

```text
SQSTM1 · primary case
CALCOCO2/NDP52 · secondary candidate
OPTN · context control
```

删减文字优先于缩小字号。

## 6. 受体证据编码

同时使用形状、颜色和直接标签。

| 受体 | 形状 | 颜色 | 证据角色 |
|---|---|---|---|
| SQSTM1 | 实心圆 | `#D95F59` | primary computational case |
| CALCOCO2/NDP52 | 实心菱形 | `#C8911B` | secondary candidate |
| OPTN | 空心方形 | `#7664A3` | mechanistic-context control |

规则：

- SQSTM1 最多只比另外两个符号大 10%–15%。
- OPTN 空心方形使用 1.5–2 px 描边，不能浅到看不见。
- 三个受体之间不画因果箭头。
- 受体与线粒体之间不画实线结合关系。
- 如需表达待验证关系，只允许很浅的短虚线，并明确写 `receptor context`。

## 7. 色彩与字体

### 色彩

| 用途 | 颜色 |
|---|---|
| 背景 | `#FBFAF7` |
| 主标题 | `#17324D` |
| 正文 | `#4D6074` |
| 内皮主色 | `#2E918A` |
| 浅青区域底色 | `#EDF7F5` |
| SQSTM1 / 应激强调 | `#D95F59` |
| 选择性清除辅助 | `#3F6FA5` |
| CALCOCO2 | `#C8911B` |
| OPTN | `#7664A3` |
| 限制信息 | `#7C8796` |
| 分隔线 | `#D9E1E6` |

避免高饱和红绿对立。颜色不是唯一信息编码。

### 字体

- 首选 Arial；替代 Helvetica 或 Aptos。
- 标题 42–48 px，Semi Bold。
- 主结论 25–29 px，Medium。
- 区域标题 19–22 px，Bold，全大写但不拉宽字距。
- 正文 18–21 px。
- 边界和数据源 16–18 px。
- 最终任何文字不低于 15 px。

## 8. 数据库与计算图标

### 数据库

- 只画普通圆角胶囊和通用数据库圆柱图标。
- 不使用 NCBI、GEO、SRA 官方 logo。
- 三个 accession 胶囊大小相同。

### 虚拟敲除

- 使用简洁网络节点 + 一条被弱化的连接，旁边标 `in silico`。
- 不使用剪刀、CRISPR、培养皿、实验小鼠或删除 DNA 的图标。
- 图标只占底栏高度约 48–60 px。

### 血清模型

- 使用非常简洁的 ROC 坐标和一条曲线。
- 不画完整统计坐标、不放置信区间、不模拟主文图。
- `AUC 0.870` 旁必须同时出现 `Internal evaluation only`。

## 9. 推荐图层结构

```text
GA_1920x1080
├── 00_BG
├── 01_Header
│   ├── Title
│   ├── Take_home_line
│   ├── Study_scale
│   └── Hypothesis_badge
├── 02_Data_input
│   ├── Dataset_chips
│   ├── Participant_aware_label
│   └── Spatial_context_chip
├── 03_Multiscale_hero
│   ├── Proximal_femur_asset
│   ├── Femoral_head_ROI
│   ├── ROI_to_niche_lines
│   ├── Microvascular_niche_asset
│   ├── Niche_to_EC_lines
│   └── Endothelial_cell_asset
├── 04_Main_findings
│   ├── Inflammatory_program_label
│   └── Receptor_context_heading
├── 05_Receptor_evidence
│   ├── SQSTM1_primary
│   ├── CALCOCO2_secondary
│   ├── OPTN_context_control
│   └── Evidence_boundary
├── 06_Computational_extensions
│   ├── Virtual_KO_summary
│   └── Serum_summary
└── 07_Next_experiments
```

## 10. 拼图顺序

1. 创建 1920 × 1080 Frame、背景和 12 列网格。
2. 放置标题、副标题和研究规模。
3. 放置股骨头、微血管放大圈和内皮细胞三个素材。
4. 添加两组各两条放大虚线并检查没有残留旧线。
5. 添加炎症/干扰素主要结果文字。
6. 添加三个受体证据等级。
7. 添加底部虚拟敲除和血清模型边界。
8. 添加下一步实验条带。
9. 统一对齐、文字行数、留白和视觉重量。
10. 导出并检查 PNG/PDF。

## 11. 科学质控

- [ ] SQSTM1 写为 `primary computational case`，不是关键驱动基因。
- [ ] CALCOCO2/NDP52 写为 `secondary transcript candidate`。
- [ ] OPTN 写为 `mechanistic-context control` 且当前数据未优先支持。
- [ ] 三个受体没有被画成已确认通路。
- [ ] 虚拟敲除明确为 `in silico`。
- [ ] 写明没有 replicated nuclear-gene FDR hit。
- [ ] 写明 matched comparators 未建立 SQSTM1 特异性。
- [ ] 血清 AUC 写明 `Internal evaluation only`。
- [ ] GSE284089 写明 `1 OA section` 和 `Anatomical context only`。
- [ ] 图中没有因果、诊断、空间验证或治疗结论。

## 12. 视觉质控

- [ ] 第一眼先看到股骨头—微血管—内皮细胞三级放大。
- [ ] 第二眼看到三个受体的证据等级。
- [ ] 底部计算扩展没有抢占主视觉。
- [ ] 没有文字压框、压线或遮挡图像。
- [ ] 没有第三条残留放大线。
- [ ] 股骨头 ROI 不直接连接内皮细胞。
- [ ] 每段正文不超过 3 行，每行尽量不超过 55 个英文字符。
- [ ] 缩放至 50% 时标题、三个基因名和核心边界仍清晰。
- [ ] 色彩以外还有形状和直接标签。
- [ ] 没有大片无意义空白，也没有塞满所有空隙。

## 13. 输出规范

- Figma 主文件保留 1920 × 1080 可编辑 Frame。
- 投稿 PNG：推荐 3840 × 2160 px，sRGB。
- 同时导出 PDF 供矢量检查。
- TIFF 从最终高分辨率导出文件转换，不从截图转换。
- PNG 至少满足 1100 px 宽、560 px 高。
- 导出后检查透明素材是否出现白底、字体是否替换、线条是否裁断。

## 14. 最终文件建议

```text
graphical_abstract_genes_receptor_context_figma_v1.fig
graphical_abstract_genes_receptor_context_v1.png
graphical_abstract_genes_receptor_context_v1.pdf
graphical_abstract_genes_receptor_context_v1.tif
```

在用户确认前，不覆盖投稿目录中的旧版 `graphical_abstract.*`。
