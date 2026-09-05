# 【归档合并版，请勿直接投喂】Nano Banana 2 + Figma 图形摘要说明书

> 本文件已拆分。后续请分别使用：
> 1. `nanobanana2_assets_receptor_context.md`：只用于 Nano Banana 2 生成无文字素材。
> 2. `figma_graphical_abstract_receptor_context.md`：只用于 Figma 拼图、标注和导出。
> 本合并版仅保留作历史记录，避免将版式文字误交给图像生成模型。

> 版本：Genes 转投稿 / receptor-context revision / 2026-09-03
> 用途：先用 Nano Banana 2 生成难以手绘的**无文字生物学插画素材**，再由 Figma 完成拼图、文字、符号、箭头和出版输出。
> 重要：**不要把整份 Markdown 一次性投给 Nano Banana。** 请按“素材 1、2、3”逐个生成，每次只复制对应的英文提示词和负面提示词。

---

## 1. 这张图真正要表达什么

### 一句话结论

Participant-aware reanalysis of public ONFH data prioritizes inflammatory programs and a **heterogeneous endothelial selective-autophagy receptor context**, with **SQSTM1** as the primary computational case, **CALCOCO2/NDP52** as a secondary candidate, and **OPTN** as a mechanistic-context control; experimental validation remains required.

中文含义：

> 对公共 ONFH 数据进行参与者层级再分析后，炎症相关程序较为稳定，而线粒体选择性自噬受体证据呈异质性：SQSTM1 是主要计算案例，CALCOCO2/NDP52 是次级候选，OPTN 是有机制关联但未被当前数据优先支持的对照；这些结果仍需实验验证。

### 图形摘要不应表达的内容

- 不得把 SQSTM1 画成已被证实驱动 ONFH 的因果靶点。
- 不得把 SQSTM1、CALCOCO2/NDP52 和 OPTN 画成一条已确定的线性通路。
- 不得表现“SQSTM1 下调必然导致线粒体损伤、炎症或血管坏死”。
- 不得把一次 OA 空间切片画成 ONFH 空间验证。
- 不得把血清 AUC 0.870 画成已经验证的诊断工具。
- 不得把计算虚拟敲除画成真实 CRISPR、动物或细胞敲除实验。
- 不得出现药物、配体、分子对接或治疗箭头。
- RNA 表达不能被等同于蛋白活性、线粒体自噬通量或线粒体功能。

---

## 2. 推荐的整体视觉结构

### 版式类型

采用 **schematic-led asymmetric composition**：一个中央生物学主视觉，加少量证据分层信息，而不是多个同等大小的统计框组成的仪表盘。

### 阅读路径

```text
公共数据与参与者层级分析
          ↓
股骨头 → 骨小梁微血管环境 → 骨微血管内皮细胞
          ↓
炎症/干扰素程序较稳定 + 选择性自噬受体证据异质
          ↓
SQSTM1（主要计算案例） / CALCOCO2（次级候选） / OPTN（机制背景对照）
          ↓
计算扩展受限 → 需要平行 loss/rescue 与功能实验
```

### 推荐 Figma 画布

- 主画布：1920 × 1080 px，16:9。
- 安全边距：四周至少 72 px。
- 推荐导出：PNG 3840 × 2160 px，另存 PDF；最终再按投稿系统需要生成 TIFF。
- 背景：暖白色 `#FBFAF7`，不使用纹理背景。
- 不在画面中写 “Graphical Abstract”。
- 生物插画约占主体面积的 55%–60%，文字与计算结果约占 40%–45%。
- 避免大段文字，也避免为了“简洁”留下大片无意义空白。

### 粗略线框

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ PARTICIPANT-AWARE ONFH REANALYSIS                               [HYPOTHESIS] │
│ 3 scRNA-seq cohorts · 18 reported participants · 13,426 endothelial cells   │
├───────────────┬───────────────────────────────────┬──────────────────────────┤
│ DATA / METHOD │            BIOLOGICAL HERO        │ EVIDENCE HIERARCHY       │
│ database chips│ femoral head → niche → EC close-up│ ● SQSTM1 primary case    │
│ participant   │                                   │ ◆ CALCOCO2 secondary     │
│ level filter  │ inflammatory/interferon programs  │ □ OPTN context control   │
│ OA spatial:   │ + heterogeneous receptor context  │                          │
│ context only  │                                   │ No causal mechanism      │
├───────────────┴───────────────────────────────────┴──────────────────────────┤
│ vKO: no replicated nuclear FDR hit | serum AUC 0.870, internal only         │
│ NEXT: matched tissue → parallel loss/rescue → flux/function → angiogenesis  │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 视觉主次

1. **第一视觉中心**：股骨头—微血管环境—内皮细胞的逐级放大关系。
2. **第二视觉中心**：三个受体的证据等级，而不是单独突出 SQSTM1。
3. **第三视觉层**：虚拟敲除和血清分类结果，仅作窄条辅助证据。
4. 空间数据只保留为一个很小的“anatomical context only”提示，不另画大型空间图。

---

## 3. 统一画风与颜色

所有 Nano Banana 素材必须看起来来自同一套医学期刊插画系统。

### 画风

- Clean scientific editorial illustration.
- Flat vector-like rendering with restrained, soft dimensional shading.
- Clear anatomy, smooth contours, consistent medium line weight.
- Modern biomedical review-figure aesthetic.
- No photorealism, no glossy 3D rendering, no cartoon face, no decorative texture.
- Objects must remain legible when reduced to approximately one third of the canvas height.

### 调色板

| 用途 | 颜色 |
|---|---|
| 暖白背景 | `#FBFAF7` |
| 骨组织 | `#E8D1A6`, 阴影 `#C8A66F` |
| 软骨 | `#C7D9E7`, 描边 `#8198B5` |
| 内皮主色 | `#2E918A`, 浅填充 `#DDF1EE` |
| 细胞核 | `#B8C7CF` |
| 应激线粒体 / SQSTM1 强调 | `#D95F59` |
| 选择性清除 / 冷色线粒体 | `#3F6FA5` |
| CALCOCO2 次级候选 | `#C8911B` |
| OPTN 机制背景对照 | `#7664A3` 或中性灰紫 |
| 边界说明 | `#7C8796` |
| 主标题 | `#17324D` |

说明：受体的证据等级符号和准确颜色由 Figma 添加。Nano Banana 生成的生物素材内不要写基因名。

---

## 4. 参考素材使用方法

### 建议上传给 Nano Banana 的参考图

1. 已认可的简化股骨头：
   `figures/source/proximal_femur_simplified_no_vessels_v1.png`

2. 如果要保持上一版的整体色调，可上传当前图形摘要，仅声明“只参考色彩和线条，不参考布局与文字”：
   使用本地投稿目录中的当前图形摘要。

3. 若手头有已经认可的圆形骨小梁微血管图，请把它作为“anatomy/style reference”单独上传。

### 不要这样使用参考图

- 不要让模型照搬旧图形摘要的版式。
- 不要让模型重绘旧图中的文字、数据库标识或统计数字。
- 不要上传期刊论文截图并要求临摹，以免产生版权和相似性问题。
- NCBI、GEO、SRA 等数据库只在 Figma 中以普通文字胶囊标签表示，不让模型生成官方 logo。

---

## 5. 素材 1：简化股骨头定位图

### 是否必须重新生成

已有素材可以直接复用。只有当线条风格与新生成素材明显不一致时才重新生成。

### 推荐尺寸

- 1200 × 1600 px 竖向构图。
- 透明背景 PNG。
- 完整物体四周至少留 8% 安全边距，不能裁掉股骨头或股骨干。

### Nano Banana 英文提示词

```text
Create one isolated scientific illustration asset for a biomedical journal graphical abstract.

Subject: a simplified anatomically recognizable human proximal femur, shown in a clean three-quarter view. Include the femoral head, a thin smooth bluish articular-cartilage cap, femoral neck, greater trochanter, and a short proximal shaft. The femoral head and neck must be immediately recognizable, but the illustration should remain simplified and uncluttered.

Style: clean scientific editorial illustration, flat vector-like rendering with restrained soft shading, smooth contours, consistent medium line weight, high-end biomedical review-figure aesthetic. Use warm ivory and muted beige for bone, pale desaturated blue for cartilage, and subtle tan shadows. Preserve the visual style of the uploaded reference image without copying any text or surrounding layout.

Composition: one complete proximal femur centered on a transparent background, with generous padding around all edges. Leave an unobstructed area at the superior femoral head-neck junction where a small ROI marker can later be added in Figma.

Do not draw blood vessels, capillaries, bone-marrow cells, lesions, necrosis, labels, arrows, text, numbers, icons, logos, watermarks, frames, circles, callouts, or background scenery. Do not add a pelvis, acetabulum, knee, or whole skeleton.

Deliverable: high-resolution transparent-background PNG, one isolated object only, no baked-in text.
```

### 负面提示词

```text
photorealistic X-ray, CT scan, surgery, hip implant, pelvis, acetabulum, full skeleton, vascular tree, red arteries, necrotic lesion, fracture, tumor, inflammation glow, labels, letters, numbers, arrows, dashed lines, watermark, logo, white rectangle, cropped anatomy, exaggerated cartoon proportions, glossy 3D plastic rendering
```

### 合格标准

- 股骨头、股骨颈和大转子关系自然。
- 软骨只覆盖股骨头关节面，不覆盖整个骨。
- 不出现血管，后续由放大圈解释微血管环境。
- ROI 圆点和两条引导线由 Figma 添加，不能烘焙在图里。

---

## 6. 素材 2：圆形“骨小梁微血管环境”放大圈

这是最需要 Nano Banana 帮助的难画素材。它只负责表现解剖环境，不承担结论表达。

### 推荐尺寸

- 1400 × 1400 px 正方形。
- 圆形构图，圆外透明。
- 圆内主要结构尽量填满 80%–88%，不要空洞。

### Nano Banana 英文提示词

```text
Create a circular cutaway scientific illustration asset showing the human femoral-head trabecular-bone microvascular niche for a biomedical journal graphical abstract.

Inside the circular field, show two or three simplified warm-ivory trabecular bone struts framing a single narrow capillary that passes diagonally through the marrow space. The capillary must have a thin continuous endothelial lining made of elongated, flattened endothelial cells with small oval nuclei. Place only three or four red blood cells inside the capillary lumen. Maintain a clear lumen and a biologically plausible relationship between the endothelial wall, lumen, red blood cells, and surrounding trabecular bone.

The capillary should be a quiet anatomical context, not a large artery and not a free-floating tube. Keep the vessel wall thin. The surrounding bone struts should be smooth, porous and simplified, without excessive cellular detail. Do not show sprouting angiogenesis, immune-cell infiltration, thrombosis, bleeding, necrosis, or diseased tissue.

Style: clean scientific editorial illustration; flat vector-like forms with subtle soft shading; consistent medium line weight; warm ivory bone, pale cyan marrow background, muted teal endothelium, and restrained brick-red erythrocytes. Match the palette and line quality of the supplied proximal-femur reference. The final asset should remain clear at small size.

Composition: a balanced circular field, with the capillary running diagonally from lower left to upper right and trabecular bone occupying the outer quadrants. Keep all important structures inside the circle. Outside the circle must be transparent.

Do not add labels, text, numbers, arrows, dotted lines, ROI markers, gene symbols, mitochondria, autophagosomes, database logos, figure borders, legends, or watermarks.

Deliverable: high-resolution PNG with transparency outside the circular field, no baked-in text.
```

### 负面提示词

```text
large artery, thick muscular vessel wall, vein valves, branching vascular tree, red and blue vessel network, vessel loop, broken vessel, leaking blood, clot, thrombus, angiogenic sprout, endothelial cells outside the wall, star-shaped endothelial cells, pericytes drawn as spikes, inflammatory cells, tumor, necrosis, osteoblast labels, osteoclast labels, dense bone cells, mitochondria, DNA helix, gene names, text, arrows, watermark, photorealism, glossy 3D, excessive empty space, cropped circular field
```

### 合格标准

- 红细胞在血管腔内，不在内皮细胞内部。
- 内皮细胞是连续、扁平、沿血管排列的壁结构，不是星形细胞。
- 骨小梁与微血管同时可见，但微血管是主角。
- 圆形边缘完整，方便在 Figma 中加 4–6 px 青绿色描边。
- 不含任何文字、虚线或标签。

---

## 7. 素材 3：骨微血管内皮细胞放大剖面

该素材是整张图的视觉主角。它需要说明“内皮细胞与线粒体背景”，但不能画成已证实的自噬机制。

### 推荐尺寸

- 1800 × 1100 px 横向构图。
- 透明背景 PNG。
- 细胞完整可见，四周留 6%–8% 安全边距。

### Nano Banana 英文提示词

```text
Create one isolated bone-microvascular endothelial-cell cutaway for a high-impact biomedical journal graphical abstract.

Subject: a single elongated human bone-microvascular endothelial cell shown as a clean side-view cutaway. The cell should have a gently tapered oval-spindle outline, a thin muted-teal plasma membrane, pale aqua cytoplasm, and one flattened light blue-grey nucleus. Include exactly four clearly separated mitochondria: two coral-red mitochondria representing a stress-associated state and two desaturated blue mitochondria representing a selective-clearance hypothesis. Each mitochondrion must have anatomically plausible inner cristae and must remain fully inside the cytoplasm.

Near only one blue mitochondrion, draw a subtle incomplete double-membrane crescent or dashed contour suggesting a structure selected for future testing. It must not be shown as a completed autophagosome, lysosome, degradation event, or proven mitophagy mechanism. Keep this feature understated.

Add a few small, generic, well-spaced cytoplasmic vesicles in muted coral and blue solely for visual balance. Do not depict molecular binding, protein complexes, receptor recruitment, DNA release, inflammatory signaling, apoptosis, or causal arrows. The three receptor evidence badges and all gene labels will be added later in Figma, so leave clean open space along the right side of the cell.

Style: clean scientific editorial illustration, flat vector-like rendering with restrained soft dimensional shading, smooth contours, consistent medium line weight, minimal but anatomically coherent, suitable for a Nature-style schematic-led composite. Use a limited palette: muted teal membrane, pale aqua cytoplasm, light blue-grey nucleus, coral-red and slate-blue mitochondria. No photorealism and no glossy 3D effects.

Composition: one complete horizontally oriented cell centered on a transparent background. The cell should occupy most of the canvas while remaining uncropped. Keep organelles separated and avoid visual clutter.

Do not add red blood cells, capillary lumen, bone, extracellular matrix, immune cells, labels, text, gene symbols, numbers, arrows, dashed callout lines, receptor icons, pathway diagrams, legends, logos, watermarks, or colored background panels.

Deliverable: high-resolution transparent-background PNG, one isolated cell only, no baked-in text.
```

### 负面提示词

```text
red blood cells inside the endothelial cytoplasm, vessel tube, epithelial cell, neuron, fibroblast, round immune cell, multiple nuclei, giant nucleus, mitochondria outside the cell, tightly overlapping organelles, mitophagosome fully enclosing mitochondrion, lysosomal degradation, ubiquitin chains, protein labels, SQSTM1, CALCOCO2, NDP52, OPTN, PINK1, Parkin, DNA helix, cGAS-STING arrows, apoptosis fragments, necrosis, gene network, labels, letters, numbers, watermark, logo, photorealistic microscopy, glossy 3D render, cropped cell
```

### 合格标准

- 是单个内皮细胞，不是完整血管。
- 线粒体数量少且相互分开，不能像密集装饰图案。
- 红色/蓝色只表示示意状态，图中不写“damaged”或“cleared”。
- 只允许一个未闭合、低调的膜性轮廓，不能直接画出“线粒体被吃掉”。
- 右侧保留足够空间，供 Figma 放置三种受体证据标记。

---

## 8. 素材 4（可选）：孤立的“线粒体 + 未闭合膜结构”小图标

只有素材 3 中的局部结构不够清晰时才生成。不要用它单独宣称机制。

### 推荐尺寸

- 1000 × 800 px。
- 透明背景 PNG。

### Nano Banana 英文提示词

```text
Create a small isolated scientific illustration asset for a hypothesis-generating biomedical figure: one desaturated slate-blue mitochondrion with anatomically plausible cristae, adjacent to a subtle incomplete double-membrane crescent. The crescent must remain clearly open and must not engulf the mitochondrion. This is a neutral visual placeholder for selective-autophagy receptor testing, not a depiction of confirmed mitophagy.

Use a clean flat vector-like scientific editorial style, restrained soft shading, consistent medium line weight, and a transparent background. Leave generous empty space around the object so that three evidence-tier badges and labels can later be added in Figma.

Do not include proteins, receptors, ubiquitin chains, lysosomes, autophagosomes, degradation particles, arrows, labels, text, gene names, numbers, legends, logos, watermarks, or background panels.

Deliverable: high-resolution transparent PNG, one isolated mitochondrion and one open membrane crescent only.
```

### 负面提示词

```text
completed autophagosome, engulfed mitochondrion, lysosome fusion, degradation, mitophagy confirmed, protein complex, receptor, ubiquitin, PINK1, Parkin, SQSTM1, NDP52, OPTN, arrows, causal pathway, DNA release, text, labels, logo, watermark, photorealistic microscopy, glossy 3D render
```

---

## 9. Nano Banana 不负责生成的元素

以下内容必须在 Figma 中制作，避免 AI 产生错字、错误数字或过强的机制暗示：

- 所有英文标题、标签和说明。
- `SQSTM1`、`CALCOCO2/NDP52`、`OPTN` 三个名称。
- 三种证据等级符号。
- 数据集编号和参与者/细胞数量。
- 所有实线、虚线、箭头、括号和引导线。
- 虚拟敲除与 ROC/AUC 的简化图标。
- 数据库胶囊标签。
- 图例和证据边界声明。

---

## 10. Figma 中使用的精确文字

文字尽量短，所有结论用“prioritized / candidate / context / internal only”等限定词。

### 顶部标题与副标题

**Title**

```text
Participant-aware ONFH reanalysis
```

**Take-home line**

```text
Inflammatory programs and a heterogeneous endothelial selective-autophagy receptor context
```

**Study scale**

```text
3 scRNA-seq cohorts · 18 reported participants · 13,426 endothelial cells
```

右上角小标签：

```text
HYPOTHESIS-GENERATING
```

### 左侧数据与分析层

```text
PUBLIC FEMORAL-HEAD DATA
Participant-aware inference
```

数据库胶囊标签：

```text
SRP361778
GSE169396
GSE290411
```

空间数据小标签：

```text
GSE284089 · 1 OA section
Anatomical context only
```

### 中央主要发现

```text
INFLAMMATORY / INTERFERON PROGRAMS
Directionally stable under participant omission
```

```text
HETEROGENEOUS RECEPTOR EVIDENCE
```

受体证据等级：

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

### 底部证据边界与下一步

```text
RNA abundance ≠ protein activity or autophagic flux
No receptor-specific causal mechanism established
```

```text
NEXT
Matched tissue → parallel SQSTM1/CALCOCO2 loss-rescue ± OPTN comparator
→ mitochondrial and autophagic-flux assays → endothelial survival and angiogenesis
```

### 可选更短版本

如果 Figma 试排后文字仍拥挤，优先删去数据库编号和受体的第二行说明，保留：

```text
SQSTM1 · primary case
CALCOCO2/NDP52 · secondary candidate
OPTN · context control
```

不要通过缩小到难以阅读的字号来容纳更多文字。

---

## 11. Figma 证据等级编码

受体身份必须同时由**形状 + 颜色 + 文字**编码，不能只靠红绿颜色。

| 受体 | 视觉符号 | 颜色 | 含义 |
|---|---|---|---|
| SQSTM1 | 实心圆 | 珊瑚红 `#D95F59` | primary computational case |
| CALCOCO2/NDP52 | 实心菱形 | 金色 `#C8911B` | secondary candidate |
| OPTN | 空心方形，1.5–2 px 描边 | 灰紫 `#7664A3` | mechanistic-context control |

规则：

- 三个符号尺寸接近，SQSTM1 最多只比另外两个大 10%–15%。
- 不画从 SQSTM1 指向 CALCOCO2 或 OPTN 的因果箭头。
- 不画三个受体直接连接线粒体或自噬体的实线。
- 如需表达“待验证关系”，只用很浅的灰色短虚线，并在旁写 `receptor context`。
- `OPTN` 不要做成浅到看不见；它是明确的对照角色，不是被删除的基因。

---

## 12. Figma 拼图坐标建议（1920 × 1080）

### 全局

- 页面边距：72 px。
- 12 列网格，列间距 24 px。
- 主标题区域：y = 58–150 px。
- 主体区域：y = 180–810 px。
- 底部证据边界/下一步：y = 842–1015 px。

### 三个区域

1. **左侧输入与定位**：x = 72–470 px。
   股骨头约高 390–430 px；圆形微血管放大圈直径约 235–260 px，压在股骨头右上方，但不遮住关节面。

2. **中央主视觉**：x = 440–1280 px。
   内皮细胞约宽 660–760 px；圆形放大圈与内皮细胞之间使用两条浅灰蓝虚线，明确表达第二级放大。

3. **右侧证据等级**：x = 1300–1848 px。
   三个受体项目纵向排列；下面放一条简短证据边界。不要把右区做成三个巨大卡片。

### 放大关系

- 股骨头 ROI 圆 → 微血管环境放大圈：两条浅灰蓝虚线。
- 微血管环境放大圈 → 单个内皮细胞：两条浅灰蓝虚线。
- 股骨头 ROI **不得直接连接内皮细胞**。
- 引导线应终止在对象边缘，不能穿过文字、细胞或红细胞。

### 底部条带

- 左 56%：SQSTM1 虚拟敲除边界。
- 中 20%：血清 AUC 0.870，标注 `internal only`。
- 右 24%：下一步实验，用简单箭头串联，不使用实验照片风格图标。

---

## 13. 推荐的 Figma 图层结构

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

---

## 14. Nano Banana 生成与筛选流程

### 第一轮：只看结构

每个素材先生成 4 个候选，不追求细节，筛选标准：

- 解剖结构是否可信。
- 轮廓是否简洁。
- 缩小后是否仍可辨认。
- 是否没有文字、假字、箭头和多余器官。

### 第二轮：修解剖错误

使用局部编辑，只修一个问题，例如：

- “Move all erythrocytes into the capillary lumen.”
- “Replace star-shaped wall cells with flattened continuous endothelial cells.”
- “Reduce the number of mitochondria to exactly four and separate them.”
- “Open the double-membrane crescent; do not enclose the mitochondrion.”

不要在一次修图中同时修改结构、颜色、构图和光影。

### 第三轮：统一风格

将最终股骨头、微血管圈和内皮细胞图同时作为参考，仅要求：

```text
Match line weight, muted palette, shading restraint, and scientific editorial style across these assets. Preserve each asset's anatomy and composition. Do not add text or new objects.
```

### 第四轮：放大与透明背景

- 输出最终高分辨率 PNG。
- 确认 alpha 透明背景真实存在，不是白色棋盘或白色矩形。
- 不要使用 JPEG 作为 Figma 主素材。

### 文件命名

```text
NB2_GA_receptor_context_asset01_femur_v1.png
NB2_GA_receptor_context_asset02_microvascular_niche_v1.png
NB2_GA_receptor_context_asset03_endothelial_cell_v1.png
NB2_GA_receptor_context_asset04_mito_open_membrane_v1.png
```

每次局部修订只递增版本号，不覆盖上一版。

---

## 15. 最终科学与出版质控

### 科学边界

- [ ] 主视觉没有把三个受体画成已验证通路。
- [ ] SQSTM1 明确是 `primary computational case`，不是“关键驱动基因”。
- [ ] CALCOCO2/NDP52 明确是 `secondary candidate`，没有 FDR 支持的暗示。
- [ ] OPTN 明确是 `mechanistic-context control`，且当前数据未优先支持。
- [ ] 虚拟敲除写明 `in silico`，并写明没有跨供体重复的核编码基因 FDR 命中。
- [ ] 血清 AUC 写明 `internal evaluation only`。
- [ ] OA 空间数据只写 `anatomical context only`。
- [ ] 图中没有把 RNA 表达当作蛋白活性、线粒体自噬通量或功能验证。

### 版式与可读性

- [ ] 先看到股骨头—微血管—内皮细胞，再看到受体证据等级。
- [ ] 没有字压框、字压线、标签互相遮挡。
- [ ] 缩到单栏预览时，标题、三个基因名和核心边界仍清楚。
- [ ] 没有超过 3 行的正文块；每行尽量不超过 55 个英文字符。
- [ ] 字体统一使用 Arial、Helvetica 或 Aptos；不混用衬线体。
- [ ] 颜色之外还有形状和文字冗余编码。
- [ ] 不使用 NCBI/GEO 官方 logo 或其他可能涉及商标的图形。
- [ ] 图不是主文某一幅图的简单复制，也不是把多个主图缩小拼接。

### 导出

- [ ] PNG 至少 1100 px 宽、560 px 高；推荐 3840 × 2160 px。
- [ ] 同时保留可编辑 Figma 文件和一份 PDF。
- [ ] 检查 PNG、PDF 中没有裁切、透明度错误和字体替换。
- [ ] 提交前以 100% 和 50% 两种显示比例检查文字与引导线。

---

## 16. 给 Nano Banana 的最后一句统一提醒

每次生成素材时可在提示词末尾追加：

```text
Scientific accuracy and clean isolation are more important than decorative detail. This asset will be assembled with exact labels in Figma. Therefore, generate no text, no gene names, no arrows, no logos, and no causal pathway claims.
```
