# Nano Banana 2 素材生成提示词

> 项目：ONFH 内皮选择性自噬受体背景图形摘要
> 版本：Genes receptor-context revision / 2026-09-03
> 本文件只用于生成**无文字生物学插画素材**。图题、基因名、数字、箭头、证据等级和数据库信息全部留给 Figma。

## 使用方法

1. 不要把整份文件一次性投给 Nano Banana 2。
2. 按“素材 1、素材 2、素材 3”分别生成；素材 4 仅在需要时生成。
3. 每个素材先生成 4 个候选，再选一个做局部修正。
4. 优先保证解剖正确、对象完整和透明背景，不追求装饰细节。
5. Nano Banana 的图里不得出现文字、基因名、箭头、统计数字或数据库 logo。

## 科学边界

素材服务于以下结论，但不能直接把结论画成已证实机制：

> Participant-aware reanalysis prioritizes inflammatory programs and a heterogeneous endothelial selective-autophagy receptor context. SQSTM1 is the primary computational case, CALCOCO2/NDP52 is a secondary candidate, and OPTN is a mechanistic-context control. Experimental validation remains required.

禁止表达：

- SQSTM1 已被证实驱动 ONFH。
- SQSTM1、CALCOCO2/NDP52 和 OPTN 构成已确定的线性通路。
- SQSTM1 下调必然导致线粒体损伤、炎症或血管坏死。
- RNA 表达等同于蛋白活性、线粒体自噬通量或线粒体功能。
- 计算虚拟敲除等同于 CRISPR、动物或细胞敲除实验。
- OA 空间切片是 ONFH 空间验证。
- 血清模型是已验证诊断工具。

## 统一画风

- Clean scientific editorial illustration.
- Flat vector-like rendering with restrained soft shading.
- Smooth contours and consistent medium line weight.
- Modern biomedical review-figure aesthetic.
- No photorealism, glossy 3D rendering, decorative texture or cartoon faces.
- 缩小到图形摘要约三分之一高度时仍能辨认。

### 统一颜色

| 对象 | 建议颜色 |
|---|---|
| 骨组织 | `#E8D1A6`，阴影 `#C8A66F` |
| 软骨 | `#C7D9E7`，描边 `#8198B5` |
| 内皮膜 | `#2E918A` |
| 内皮胞质 | `#DDF1EE` |
| 细胞核 | `#B8C7CF` |
| 应激相关线粒体 | `#D95F59` |
| 选择性清除假设线粒体 | `#3F6FA5` |

受体证据等级的颜色由 Figma 添加，不要让 Nano Banana 在细胞内生成受体图标。

## 可上传的参考图

### 股骨头参考图

`E:\ONFH2026dry\figures\source\proximal_femur_simplified_no_vessels_v1.png`

### 仅供配色参考的旧图

`E:\ONFH2026dry\outputs\onfh-genes-revision-v1-20260903\submission_ready\graphical_abstract.png`

如果上传旧图，必须告诉模型：

```text
Use the uploaded graphical abstract only as a reference for muted palette and line quality. Do not reproduce its layout, labels, text, arrows, evidence hierarchy, or composition.
```

不要上传已发表论文插图让模型临摹。

---

## 素材 1：简化股骨头

已有参考素材可以直接复用。只有与新素材风格明显不一致时才重新生成。

### 输出要求

- 1200 × 1600 px，竖向。
- 透明背景 PNG。
- 完整物体四周至少保留 8% 安全边距。

### 提示词

```text
Create one isolated scientific illustration asset for a biomedical journal graphical abstract.

Subject: a simplified anatomically recognizable human proximal femur, shown in a clean three-quarter view. Include the femoral head, a thin smooth bluish articular-cartilage cap, femoral neck, greater trochanter, and a short proximal shaft. The femoral head and neck must be immediately recognizable, but the illustration should remain simplified and uncluttered.

Style: clean scientific editorial illustration, flat vector-like rendering with restrained soft shading, smooth contours, consistent medium line weight, and a high-end biomedical review-figure aesthetic. Use warm ivory and muted beige for bone, pale desaturated blue for cartilage, and subtle tan shadows. Preserve the visual style of the uploaded reference image without copying any text or surrounding layout.

Composition: one complete proximal femur centered on a transparent background, with generous padding around all edges. Leave an unobstructed area at the superior femoral head-neck junction where a small ROI marker can later be added in Figma.

Do not draw blood vessels, capillaries, bone-marrow cells, lesions, necrosis, labels, arrows, text, numbers, icons, logos, watermarks, frames, circles, callouts, or background scenery. Do not add a pelvis, acetabulum, knee, or whole skeleton.

Deliverable: high-resolution transparent-background PNG, one isolated object only, no baked-in text.

Scientific accuracy and clean isolation are more important than decorative detail. This asset will be assembled with exact labels in Figma. Therefore, generate no text, no gene names, no arrows, no logos, and no causal pathway claims.
```

### 负面提示词

```text
photorealistic X-ray, CT scan, surgery, hip implant, pelvis, acetabulum, full skeleton, vascular tree, red arteries, necrotic lesion, fracture, tumor, inflammation glow, labels, letters, numbers, arrows, dashed lines, watermark, logo, white rectangle, cropped anatomy, exaggerated cartoon proportions, glossy 3D plastic rendering
```

### 验收标准

- 股骨头、股骨颈和大转子关系自然。
- 软骨只覆盖股骨头关节面。
- 不出现血管、病灶、ROI 圆点或引导线。
- 没有裁边和伪透明白底。

---

## 素材 2：圆形骨小梁微血管环境放大圈

这是最需要 Nano Banana 生成的素材。

### 输出要求

- 1400 × 1400 px。
- 圆形构图，圆外透明。
- 圆内主要结构占 80%–88%，不要留下大面积空白。

### 提示词

```text
Create a circular cutaway scientific illustration asset showing the human femoral-head trabecular-bone microvascular niche for a biomedical journal graphical abstract.

Inside the circular field, show two or three simplified warm-ivory trabecular bone struts framing a single narrow capillary that passes diagonally through the marrow space. The capillary must have a thin continuous endothelial lining made of elongated, flattened endothelial cells with small oval nuclei. Place only three or four red blood cells inside the capillary lumen. Maintain a clear lumen and a biologically plausible relationship between the endothelial wall, lumen, red blood cells, and surrounding trabecular bone.

The capillary should be a quiet anatomical context, not a large artery and not a free-floating tube. Keep the vessel wall thin. The surrounding bone struts should be smooth, porous and simplified, without excessive cellular detail. Do not show sprouting angiogenesis, immune-cell infiltration, thrombosis, bleeding, necrosis, or diseased tissue.

Style: clean scientific editorial illustration; flat vector-like forms with subtle soft shading; consistent medium line weight; warm ivory bone, pale cyan marrow background, muted teal endothelium, and restrained brick-red erythrocytes. Match the palette and line quality of the supplied proximal-femur reference. The final asset should remain clear at small size.

Composition: a balanced circular field, with the capillary running diagonally from lower left to upper right and trabecular bone occupying the outer quadrants. Keep all important structures inside the circle. Outside the circle must be transparent.

Do not add labels, text, numbers, arrows, dotted lines, ROI markers, gene symbols, mitochondria, autophagosomes, database logos, figure borders, legends, or watermarks.

Deliverable: high-resolution PNG with transparency outside the circular field, no baked-in text.

Scientific accuracy and clean isolation are more important than decorative detail. This asset will be assembled with exact labels in Figma. Therefore, generate no text, no gene names, no arrows, no logos, and no causal pathway claims.
```

### 负面提示词

```text
large artery, thick muscular vessel wall, vein valves, branching vascular tree, red and blue vessel network, vessel loop, broken vessel, leaking blood, clot, thrombus, angiogenic sprout, endothelial cells outside the wall, star-shaped endothelial cells, pericytes drawn as spikes, inflammatory cells, tumor, necrosis, osteoblast labels, osteoclast labels, dense bone cells, mitochondria, DNA helix, gene names, text, arrows, watermark, photorealism, glossy 3D, excessive empty space, cropped circular field
```

### 验收标准

- 红细胞全部在血管腔内。
- 内皮细胞为连续、扁平、沿血管排列的壁结构。
- 骨小梁和微血管均可见，微血管仍是主角。
- 圆形边缘完整，圆外真实透明。
- 不含文字、虚线、标签或 mitochondria。

---

## 素材 3：骨微血管内皮细胞放大剖面

这是整张图的视觉主角，只表达内皮细胞和线粒体背景，不表达已证实机制。

### 输出要求

- 1800 × 1100 px，横向。
- 透明背景 PNG。
- 细胞完整可见，四周保留 6%–8% 安全边距。

### 提示词

```text
Create one isolated bone-microvascular endothelial-cell cutaway for a high-impact biomedical journal graphical abstract.

Subject: a single elongated human bone-microvascular endothelial cell shown as a clean side-view cutaway. The cell should have a gently tapered oval-spindle outline, a thin muted-teal plasma membrane, pale aqua cytoplasm, and one flattened light blue-grey nucleus. Include exactly four clearly separated mitochondria: two coral-red mitochondria representing a stress-associated state and two desaturated blue mitochondria representing a selective-clearance hypothesis. Each mitochondrion must have anatomically plausible inner cristae and must remain fully inside the cytoplasm.

Near only one blue mitochondrion, draw a subtle incomplete double-membrane crescent or dashed contour suggesting a structure selected for future testing. It must not be shown as a completed autophagosome, lysosome, degradation event, or proven mitophagy mechanism. Keep this feature understated.

Add a few small, generic, well-spaced cytoplasmic vesicles in muted coral and blue solely for visual balance. Do not depict molecular binding, protein complexes, receptor recruitment, DNA release, inflammatory signaling, apoptosis, or causal arrows. The three receptor evidence badges and all gene labels will be added later in Figma, so leave clean open space along the right side of the cell.

Style: clean scientific editorial illustration, flat vector-like rendering with restrained soft dimensional shading, smooth contours, consistent medium line weight, minimal but anatomically coherent, suitable for a schematic-led scientific composite. Use a limited palette: muted teal membrane, pale aqua cytoplasm, light blue-grey nucleus, coral-red and slate-blue mitochondria. No photorealism and no glossy 3D effects.

Composition: one complete horizontally oriented cell centered on a transparent background. The cell should occupy most of the canvas while remaining uncropped. Keep organelles separated and avoid visual clutter.

Do not add red blood cells, capillary lumen, bone, extracellular matrix, immune cells, labels, text, gene symbols, numbers, arrows, dashed callout lines, receptor icons, pathway diagrams, legends, logos, watermarks, or colored background panels.

Deliverable: high-resolution transparent-background PNG, one isolated cell only, no baked-in text.

Scientific accuracy and clean isolation are more important than decorative detail. This asset will be assembled with exact labels in Figma. Therefore, generate no text, no gene names, no arrows, no logos, and no causal pathway claims.
```

### 负面提示词

```text
red blood cells inside the endothelial cytoplasm, vessel tube, epithelial cell, neuron, fibroblast, round immune cell, multiple nuclei, giant nucleus, mitochondria outside the cell, tightly overlapping organelles, mitophagosome fully enclosing mitochondrion, lysosomal degradation, ubiquitin chains, protein labels, SQSTM1, CALCOCO2, NDP52, OPTN, PINK1, Parkin, DNA helix, cGAS-STING arrows, apoptosis fragments, necrosis, gene network, labels, letters, numbers, watermark, logo, photorealistic microscopy, glossy 3D render, cropped cell
```

### 验收标准

- 是一个单独的内皮细胞，不是血管。
- 只有一个核，四个线粒体均在胞质内并彼此分开。
- 不出现红细胞、骨小梁或免疫细胞。
- 只允许一个未闭合、低调的膜性轮廓。
- 右侧留有叠加受体标记的空间。
- 不含任何文字或箭头。

---

## 素材 4（可选）：线粒体与未闭合膜结构

只有素材 3 的局部不够清楚时才生成。

### 输出要求

- 1000 × 800 px。
- 透明背景 PNG。

### 提示词

```text
Create a small isolated scientific illustration asset for a hypothesis-generating biomedical figure: one desaturated slate-blue mitochondrion with anatomically plausible cristae, adjacent to a subtle incomplete double-membrane crescent. The crescent must remain clearly open and must not engulf the mitochondrion. This is a neutral visual placeholder for selective-autophagy receptor testing, not a depiction of confirmed mitophagy.

Use a clean flat vector-like scientific editorial style, restrained soft shading, consistent medium line weight, and a transparent background. Leave generous empty space around the object so that three evidence-tier badges and labels can later be added in Figma.

Do not include proteins, receptors, ubiquitin chains, lysosomes, autophagosomes, degradation particles, arrows, labels, text, gene names, numbers, legends, logos, watermarks, or background panels.

Deliverable: high-resolution transparent PNG, one isolated mitochondrion and one open membrane crescent only.

Scientific accuracy and clean isolation are more important than decorative detail. This asset will be assembled with exact labels in Figma. Therefore, generate no text, no gene names, no arrows, no logos, and no causal pathway claims.
```

### 负面提示词

```text
completed autophagosome, engulfed mitochondrion, lysosome fusion, degradation, mitophagy confirmed, protein complex, receptor, ubiquitin, PINK1, Parkin, SQSTM1, NDP52, OPTN, arrows, causal pathway, DNA release, text, labels, logo, watermark, photorealistic microscopy, glossy 3D render
```

---

## 局部修图指令

一次只修一个问题：

```text
Move all erythrocytes into the capillary lumen. Preserve every other element.
```

```text
Replace the star-shaped wall cells with flattened continuous endothelial cells. Preserve the capillary path, bone struts, palette, and composition.
```

```text
Reduce the mitochondria to exactly four and separate them clearly. Preserve the cell outline, nucleus, palette, and transparent background.
```

```text
Open the double-membrane crescent so that it does not enclose the mitochondrion. Do not add any other organelle or annotation.
```

## 跨素材风格统一指令

将最终候选同时作为参考图上传，仅执行风格统一：

```text
Match line weight, muted palette, shading restraint, and scientific editorial style across these assets. Preserve each asset's anatomy and composition. Do not add text, arrows, labels, new organelles, or background panels.
```

## 文件命名

```text
NB2_GA_receptor_context_asset01_femur_v1.png
NB2_GA_receptor_context_asset02_microvascular_niche_v1.png
NB2_GA_receptor_context_asset03_endothelial_cell_v1.png
NB2_GA_receptor_context_asset04_mito_open_membrane_v1.png
```

不要覆盖上一版本，只递增 `v2`、`v3`。

## 最终素材检查

- [ ] 无文字、伪文字、数字、箭头、图例和 logo。
- [ ] 无裁边，主体四周有安全边距。
- [ ] 透明背景不是白色矩形或棋盘图案。
- [ ] 解剖结构可信，红细胞没有进入内皮细胞胞质。
- [ ] 没有把受体关系或线粒体自噬画成已证实机制。
- [ ] 各素材线宽、阴影和颜色一致。
- [ ] 缩小后仍能辨认股骨头、毛细血管、内皮细胞和线粒体。
