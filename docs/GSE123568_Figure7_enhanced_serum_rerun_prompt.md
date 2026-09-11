# GSE123568 血清分类器增强重跑与 Figure 7 更新提示词

你是一个熟悉微阵列表达数据、嵌套交叉验证、小样本预测建模、bootstrap、模型校准和生物信息学可复现性审计的助手，精通 Python、pandas、NumPy、scikit-learn、SciPy，以及 R、pROC、ggplot2 和 patchwork。

请使用 `academic-research-suite` 的实验规划与验证流程，在本仓库根目录内，基于公开数据 GSE123568 和现有血清分类器代码，编写并实际运行增强版血清分析。目标是消除训练折外探针选择、增加包含模型重新拟合的内部不确定性评估、量化基因选择方向稳定性、检验线粒体候选空间相对于匹配随机基因集的特异性，并补充 Brier score 与探索性校准图。

本任务只增强 Figure 7 的内部血清分析。不得把结果表述为外部验证、临床诊断模型、锁定签名或因果证据。不得为了获得更高 AUC 而调整候选基因、分组、交叉验证种子、模型类型、阈值或图形展示范围。

## Material Passport

- 项目：ONFH participant-aware public-data reanalysis
- 数据集：GSE123568
- 数据类型：外周血清微阵列表达数据
- 分析性质：探索性内部预测评估
- 当前投稿目标：Scientific Reports
- 当前任务：增强版 Figure 7 分析、审计和候选图形
- 允许写入：新的分析脚本、工作流、独立结果目录和审阅版图形目录
- 禁止写入：原始公共数据、`results/figure_inputs`、`figures/final`、已冻结投稿包
- 停止门：完成分析、质控、图形和中文报告后停止；未经作者审核不得替换主稿或投稿包

## 一、目录、脚本与版本边界

1. 项目根目录固定为当前仓库根目录。
2. 原始公开数据目录通过 `--data-dir` 指定；默认使用仓库内的 `data/`。
3. 新建独立工作目录：

   `results/serum_enhanced_20260910`

4. 在该目录下建立：

   - `00_frozen_plan`：冻结方案、输入清单、哈希和参数；
   - `01_smoke`：小规模烟雾测试结果；
   - `02_formal`：正式分析结果；
   - `03_figures`：审阅版 Figure 7 和补充图；
   - `04_report`：中文报告、运行日志和输出清单；
   - `checkpoints`：置换、bootstrap 和随机基因集断点文件。

5. 新增而不是覆盖以下代码文件：

   - `analysis/serum_classifier_enhanced.py`；
   - `plotting/make_figure7_enhanced.R`；
   - `workflow/run_serum_enhanced.ps1`；
   - `qa/check_serum_enhanced.py`。

6. 只把以下现有文件作为只读参考：

   - `analysis/serum_classifier.py`；
   - `analysis/serum_paired_comparison.R`；
   - `analysis/genesets_final.json`；
   - `plotting/make_genes_revision_figures.R`；
   - `results/figure_inputs/diag_summary_v4.json`；
   - `results/figure_inputs/diag_oof_predictions_aggregated_v8.csv`；
   - `results/figure_inputs/diag_delong_model_comparison_v8.csv`。

7. 不得覆盖既有 Figure 7、既有血清结果或任何已冻结投稿目录。不得执行 Git commit、tag、push 或 release。

## 二、核对原始输入文件

必须读取并核对：

- `data/GSE123568_series_matrix.txt.gz`；
- `data/GSE123568_family.soft.gz`；
- `analysis/genesets_final.json`。

运行前计算并保存三个文件的 SHA-256。将实际路径、文件大小、修改时间、SHA-256 和解析状态写入：

`results/serum_enhanced_20260910/00_frozen_plan/input_manifest.csv`

不得在本地文件可用时重新下载 GEO 数据，不得静默更换平台注释版本。不得修改或重新压缩原始文件。

检查 series matrix 是否已经是 GEO 存档的处理后表达矩阵，记录其最小值、最大值、分位数及是否存在缺失值。不得在没有证据的情况下再次 log2 转换或重新标准化整张矩阵；若无法确认表达尺度，停止正式分析并报告。

## 三、冻结样本名单与标签

必须从 GEO 元数据解析样本和疾病标签，并与以下预期完全核对：

- 总样本数：40；
- steroid-exposed non-SONFH：10；
- SONFH：30；
- 对照：`GSM3507251` 至 `GSM3507260`，标签为 0；
- SONFH：`GSM3507261` 至 `GSM3507290`，标签为 1。

输出 `sample_manifest.csv`，至少包含：

- `sample_id`；
- `label_numeric`；
- `label_text`；
- `source_metadata_text`；
- `included`；
- `exclusion_reason`。

检查：

- 样本 ID 唯一；
- 表达矩阵列和元数据样本一一对应；
- 样本顺序重排后完全一致；
- 标签只能为 0/1；
- 不存在同一样本重复进入训练集和测试集；
- 不根据表达值、预测结果或 AUC 排除样本。

如样本数、标签数、样本 ID 或顺序不符合预期，使用显式异常停止，不得静默修正后继续。

## 四、冻结候选基因空间和比较器

候选空间来源固定为 `genesets_final.json` 中以下六个预设模块的并集：

- `Mito_fission`；
- `Mito_fusion`；
- `Mitophagy_core`；
- `mtDNA_release`；
- `cGAS_STING`；
- `YAP_mTOR`。

芯片上预期可测的77个候选基因固定为：

```text
AKT1; AMBRA1; ATG5; ATG7; BAK1; BAX; BBC3; BCL2; BCL2L1; BCL2L13; BID; BNIP3; BNIP3L; CALCOCO2; CCL5; CHUK; CXCL10; DNM1; DNM1L; DNM2; EIF4EBP1; FIS1; FKBP8; FUNDC1; GABARAP; GABARAPL2; GDAP1; GIGYF1; GIGYF2; HIF1A; IKBKB; IL6; INF2; IRF3; IRF7; LATS1; LATS2; MAP1LC3B; MCL1; MFF; MFN1; MFN2; MIEF1; MIEF2; MST1; MTCH2; MTFR1; MTFR1L; MTOR; NBR1; NFKB1; OPA1; OPTN; PHB2; PINK1; PLD6; PMAIP1; RELA; RHEB; RPS6KB1; RPTOR; SQSTM1; TAX1BP1; TBK1; TEAD1; TEAD4; TMEM173; TNF; TSC1; TSC2; ULK1; VCP; VDAC1; VDAC2; VDAC3; WWTR1; YAP1
```

Ma 等人原候选集固定为：

```text
BID; FTH1; LACTB; PDK3; RAB5IF; SOD2; SQOR
```

本平台预期可测的比较器基因固定为：

```text
BID; FTH1; LACTB; PDK3
```

必须输出 `candidate_gene_manifest.csv` 和 `ma_comparator_gene_manifest.csv`，逐基因记录：来源模块、平台是否有探针、探针数、最终是否可用及不可用原因。

如候选基因可用数不是77，或Ma比较器可用数不是4，应停止正式分析并报告，不得通过替换基因、同义词猜测或放宽规则凑足数量。

## 五、在看到新结果前冻结分析方案

在任何正式模型拟合前创建：

`00_frozen_plan/analysis_plan_frozen.json`

至少固定：

- 40个样本及标签；
- 77个候选基因；
- Ma 4基因比较器；
- 主模型、比较器模型和参数网格；
- 外层和内层交叉验证方案；
- 所有随机种子；
- 标签置换次数；
- refit-aware bootstrap次数；
- 匹配随机基因集次数；
- 主指标、次要指标和敏感性指标；
- 所有经验P值计算公式；
- 缺失、失败和丢弃重采样的处理规则；
- 预定图形结构；
- 停止条件和解释边界。

冻结参数如下：

- 基础随机种子：`20260820`；
- 外层：分层5折，重复5次；
- 外层随机状态：`0, 1, 2, 3, 4`；
- 内层：分层4折；
- C网格：`numpy.logspace(-3, 2, 16)`；
- 候选模型：L1 logistic regression，`solver="liblinear"`；
- Ma比较器：L2 logistic regression，`solver="liblinear"`；
- `class_weight="balanced"`；
- `max_iter=5000`；
- 正式标签置换：1,000次；
- 正式refit-aware OOB bootstrap：1,000个有效重采样；
- 正式表达/方差匹配随机基因空间：1,000个；
- 并行任务数默认：4；
- checkpoint间隔：每20个任务保存一次。

对 `analysis_plan_frozen.json` 计算 SHA-256。脚本重新运行时，如计划文件存在且参数不一致，应停止并要求建立新的版本目录，不得覆盖原计划。

## 六、重构探针级数据读取

不要直接调用现有 `load_gse123568()` 返回的全数据基因级矩阵作为增强版主分析输入，因为该函数在全部40个样本上按平均表达选择探针。

新代码必须先保留：

- 完整探针 × 样本表达矩阵；
- probe ID；
- Gene Symbol；
- probe到gene的原始映射；
- 无符号、重复符号和多符号记录的处理审计。

Gene Symbol处理沿用原分析的明确规则：优先使用注释中的第一个有效符号；不得根据疾病标签或差异表达选择符号。对 ` /// `、`//`、空符号和重复探针进行显式记录。

输出 `probe_annotation_audit.csv`，至少包含：

- `probe_id`；
- `raw_gene_symbol`；
- `parsed_gene_symbol`；
- `mapping_status`；
- `n_probes_for_gene`。

## 七、训练折内探针折叠

在每一个外层训练折内，独立执行探针折叠：

1. 对同一gene的全部候选探针，仅使用外层训练样本计算平均表达。
2. 选择训练折平均表达最高的探针。
3. 如平均表达完全相同，按probe ID字典序选择，保证确定性。
4. 将训练折选定的同一个probe应用到该外层测试折。
5. 外层测试样本不得参与探针选择、缺失处理、方差判断或任何参数决定。
6. 内层调参也必须通过scikit-learn `Pipeline` 实现，使 `StandardScaler` 仅在每个内层训练子折拟合。
7. 外层最终模型在全部外层训练样本上重新拟合pipeline后，才能预测外层测试样本。

保存 `probe_selection_by_outer_fold.csv`，至少包含：

- `repeat`；
- `outer_fold`；
- `gene`；
- `selected_probe`；
- `n_candidate_probes`；
- `training_mean_expression`；
- `tie_count`；
- `selection_rule`。

必须编写并运行泄漏单元测试：人为改变外层测试折的表达值后，该折训练阶段选择的probe ID必须完全不变。测试失败即停止正式分析。

## 八、增强版重复嵌套交叉验证

为候选模型和Ma比较器创建完全相同的外层分组。每次重复中，每个样本必须恰好进入一次测试折。

每个外层折执行：

1. 仅用训练折选择gene对应probe；
2. 仅用训练数据建立预处理和模型pipeline；
3. 在内层4折中分别拟合 `StandardScaler` 和logistic regression；
4. 用ROC-AUC选择C；
5. 候选模型使用L1，Ma比较器使用L2；
6. 在相同外层测试样本上输出两个模型的概率；
7. 保存所选C、收敛状态、迭代次数、系数和测试概率。

主要统计量固定为：每个样本5次外层OOF概率的平均值所计算的aggregate ROC-AUC。

次要指标固定为：

- aggregate average precision；
- aggregate Brier score；
- 5次repeat-specific AUC；
- 5次repeat-specific AP；
- 5次repeat-specific Brier score；
- 候选模型减Ma比较器的AUC差；
- 候选模型减Ma比较器的AP差；
- 候选模型减Ma比较器的Brier差，注意Brier越低越好。

不得使用外层测试标签选择模型、特征、C值、探针、概率阈值或校准方法。

输出：

- `outer_fold_assignments.csv`；
- `oof_predictions_long.csv`；
- `oof_predictions_aggregated.csv`；
- `repeat_performance.csv`；
- `aggregate_performance.csv`；
- `selected_hyperparameters.csv`；
- `model_comparison_fixed_oof.csv`。

## 九、基因选择频率与方向稳定性

候选模型共有25个外层最终拟合。对77个候选基因逐一汇总标准化尺度上的logistic回归系数。

至少计算：

- `n_outer_fits=25`；
- `n_nonzero`；
- `nonzero_frequency`；
- `n_positive`；
- `n_negative`；
- `positive_frequency_all_fits`；
- `negative_frequency_all_fits`；
- `dominant_direction`；
- `direction_consistency_given_selected = max(n_positive,n_negative)/n_nonzero`；
- 包含0值的系数中位数、Q1和Q3；
- 仅非零系数的中位数、Q1和Q3；
- 各重复中进入模型的次数；
- 是否在至少50%的外层拟合中非零；
- 是否同时满足非零频率≥50%且方向一致率≥80%。

方向定义固定为：

- 正系数：表达升高对应更高SONFH预测概率；
- 负系数：表达升高对应更低SONFH预测概率。

若基因从未被选择，方向一致率记为NA，不得记为0或100%。不得把选择频率解释为统计显著性、因果效应或独立生物标志物证据。

输出：

- `feature_coefficients_long.csv`；
- `feature_stability_summary.csv`；
- `features_per_outer_fit.csv`。

## 十、完整流程标签置换检验

保留现有标签置换思想，但增强版必须重新执行完整流程：

1. 每次在40个样本层面置换标签；
2. 按置换标签重新建立分层外层折；
3. 每个外层折重新进行训练折内探针选择；
4. 每个内层折重新拟合scaler并重新选择C；
5. 重新拟合L1模型并生成OOF概率；
6. 计算与观察分析定义完全相同的aggregate OOF AUC。

正式置换次数为1,000。经验P值固定为：

```text
p = (1 + 置换AUC >= 观察AUC的次数) / (1 + 有效置换次数)
```

不得用repeat-level AUC平均值替代观察分析的aggregate OOF AUC。不得只置换最终预测或跳过模型重新调参。

输出 `label_permutation_null.csv` 和 `label_permutation_summary.csv`。即使结果不显著，也必须保留完整输出，不得追加置换直至显著。

## 十一、包含模型重新拟合的不确定性分析

增加分层out-of-bag bootstrap，目标是评估重新抽样、探针选择、标准化、调参和变量选择共同造成的模型不确定性。

每个bootstrap执行：

1. 从30个SONFH样本中有放回抽取30次，从10个对照中有放回抽取10次。
2. 未被抽中的原始样本组成该次唯一OOB测试集。
3. 候选模型和Ma比较器使用完全相同的bootstrap训练样本和OOB测试样本。
4. 训练样本中的重复GSM ID代表bootstrap权重，但在内层交叉验证中，同一原始GSM ID的全部重复副本必须位于同一个子折。
5. 使用 `StratifiedGroupKFold` 或等价的显式group约束，group固定为原始GSM ID。
6. 内层折数取 `min(4, 训练集中两类的唯一GSM数)`；少于2折则该次bootstrap无效。
7. 探针选择仅使用bootstrap训练样本；相同原始样本的重复副本按抽样次数参与均值计算。
8. scaler、C选择和模型拟合全部在bootstrap训练流程内完成。
9. OOB测试集中必须至少包含2个SONFH和2个对照，否则该次记为无效并按预设规则重新抽取。
10. 记录无效原因；最多尝试5,000次以获得1,000个有效bootstrap。若无法获得1,000个有效重采样，停止并报告，不得降低有效性标准。

每个有效bootstrap同时记录：

- 候选模型OOB AUC、AP、Brier；
- Ma比较器OOB AUC、AP、Brier；
- 配对AUC差、AP差和Brier差；
- OOB两类样本数；
- 训练集唯一GSM数；
- 候选模型选择的基因数；
- 两个模型的最佳C。

用2.5%和97.5%分位数汇总refit-aware OOB分布。必须称为“refit-aware OOB bootstrap sensitivity interval”或“重新拟合敏感性区间”，不得无条件称为外部泛化95%置信区间。不得把1,000次bootstrap当作1,000个独立患者。

输出：

- `refit_oob_bootstrap_all.csv`；
- `refit_oob_bootstrap_invalid.csv`；
- `refit_oob_bootstrap_summary.csv`。

## 十二、匹配随机基因空间特异性检验

检验77基因线粒体候选空间是否优于同样大小、表达和方差相近的非候选基因空间。

背景基因池必须：

- 来自同一GPL15207平台；
- 有有效Gene Symbol；
- 在当前外层训练折中有有限表达和非零方差；
- 排除77个线粒体候选基因；
- 排除Ma 7基因原候选集；
- 不根据疾病标签、差异表达P值或全数据AUC筛选。

匹配必须在每个外层训练折内进行：

1. 先对平台上所有可用gene按本提示词第七节的规则，仅用当前外层训练样本选择代表probe。
2. 使用各gene在训练折选定probe上的表达值计算训练均值和log1p方差；测试折不得参与代表probe或匹配层的确定。
3. 按训练折内gene表达均值和log1p方差分别划分十分位数。
4. 形成“均值十分位 × 方差十分位”二维匹配层。
5. 对77个候选基因逐一在相同二维层中无放回抽取一个背景基因。
6. 若原层无可用基因，按预先规定的曼哈顿距离逐层扩展；优先距离较小者，距离相同按固定随机种子抽取。
7. 每个随机空间在每个外层折均保持77个不重复基因。
8. 将每个随机空间视为一个预设的匹配分析规则；由于匹配统计量来自训练数据，不同外层折允许选中不同具体背景基因，但必须保存完整映射。
9. 对每个随机空间运行与候选模型完全相同的scaler、内层调参、L1模型和外层预测流程；代表probe沿用该折训练数据确定的映射。
10. 所有随机空间使用与观察候选模型相同的外层分组和种子。

正式生成1,000个匹配随机基因空间。每个空间计算aggregate OOF AUC、AP和Brier。AUC经验P值固定为：

```text
p_specificity = (1 + 匹配随机空间AUC >= 候选空间AUC的次数) / 1001
```

必须输出匹配质量：

- 每个候选基因匹配到的背景基因；
- 原始二维层；
- 是否扩层；
- 扩展距离；
- 训练均值差；
- 训练log1p方差差；
- 每个随机空间是否完整包含77个唯一基因。

输出：

- `matched_random_space_auc_null.csv`；
- `matched_random_space_gene_map.csv.gz`；
- `matched_random_space_qc.csv`；
- `matched_random_space_summary.csv`。

如果背景池不足以按预设规则形成77基因空间，应停止并报告，不得降低基因数或允许重复基因。

## 十三、固定OOF比较与配对统计

基于两个模型在相同40个样本上的aggregate OOF概率，计算：

- 候选模型和Ma比较器的ROC-AUC、AP和Brier；
- 配对DeLong检验；
- 对固定OOF预测进行20,000次分层配对bootstrap，得到AUC差区间；
- 固定预测bootstrap和refit-aware OOB bootstrap必须分开报告。

固定预测bootstrap只反映这40个既有cross-fitted预测的样本重采样，不包含重新拟合；refit-aware OOB bootstrap包含模型重建，但OOB样本较少且仍属于同一数据集内部评估。不得把两者混成同一个“95% CI”。

## 十四、Brier score与探索性校准

使用aggregate OOF概率评估校准：

1. 计算候选模型和Ma比较器的Brier score。
2. 计算常数概率0.75对应的同数据Brier基线。
3. 可计算Brier skill score，但必须明确参照的是样本内患病率0.75。
4. 生成4个等频概率组的探索性校准图；每组标注样本数、平均预测概率和观察SONFH比例。
5. 计算探索性calibration-in-the-large和calibration slope；预测概率先截断到 `[1e-6, 1-1e-6]`。
6. 如logistic校准模型不收敛或发生完全分离，应输出NA和失败原因，不得用任意惩罚参数替代后继续宣称校准良好。

由于总样本只有40例且类别为30:10：

- 不进行Hosmer–Lemeshow显著性检验；
- 不根据本数据选择临床阈值；
- 不报告敏感度、特异度或最佳Youden阈值作为主要结果；
- 不把校准曲线解释为临床校准验证。

输出：

- `calibration_summary.csv`；
- `calibration_bins.csv`；
- `brier_comparison.csv`。

## 十五、烟雾测试与正式运行

先运行独立烟雾测试：

- 标签置换20次；
- refit-aware OOB bootstrap20个有效重采样；
- 匹配随机基因空间20个；
- 只输出到 `01_smoke`；
- 不允许把烟雾测试结果并入正式统计。

烟雾测试必须检查：

- 所有40个样本均产生OOF预测；
- 每个repeat内每个样本恰好测试一次；
- 训练集和测试集无重叠；
- 测试折表达变化不会改变训练折probe选择；
- candidate和Ma模型的外层折完全一致；
- 概率均在0至1之间且无缺失；
- 77个候选基因顺序稳定；
- 随机空间每折均含77个唯一非候选基因；
- bootstrap重复GSM不会跨内层训练/验证子折；
- 相同种子重复运行得到完全相同的烟雾测试摘要。

只有全部烟雾测试通过后，才能运行正式分析。正式运行写入 `02_formal`，不得复制或沿用烟雾测试统计量。

三个高计算量模块必须支持：

- `--resume`；
- 按任务ID去重；
- 原子化checkpoint写入；
- 中断后从最后一个完整ID继续；
- 不因重新启动而重复计数；
- 每20个任务更新进度和日志；
- 任务失败时保留错误信息，不得静默重试或丢弃。

## 十六、与旧版Figure 7结果对照

正式分析完成后，只读比较以下旧版数值：

- 候选空间aggregate OOF AUC 0.870；
- 候选空间aggregate AP 0.927；
- Ma四基因aggregate OOF AUC 0.813；
- AUC差0.057；
- paired DeLong P值0.192；
- 固定预测AUC差区间约为−0.027至0.153；
- 旧版候选空间有77个可测基因。

这些数值只是复核参照，不是新分析必须复现的目标。由于探针折叠和内层scaler拟合方式更严格，新结果允许变化。不得修改代码以逼近旧AUC。

生成 `legacy_vs_enhanced_comparison.csv`，逐项说明：

- 旧值；
- 新值；
- 绝对差；
- 相对差；
- 变化可能来自哪个预处理或重采样环节；
- 是否改变原有审慎结论。

## 十七、结果判定规则

严格采用以下解释规则：

- 标签置换经验P值小于0.05：只说明候选模型在该数据集内部的判别信息超过标签随机化流程；
- 匹配随机基因空间经验P值小于0.05：只说明预设候选空间在当前内部分析中优于表达/方差匹配的同规模随机空间；
- 两者均不显著：不得宣称候选空间具有特异性；
- candidate减Ma的配对AUC区间包含0或DeLong P值不小于0.05：不得宣称优于Ma比较器；
- refit-aware OOB区间很宽或跨越0.5：明确报告模型重新拟合不稳定；
- 基因选择频率高但方向不一致：不得列为稳定候选；
- Brier低于患病率常数模型：可描述为内部概率误差较低，但不是临床校准验证；
- 校准斜率或截距不稳定：如实报告，不得隐藏；
- 即使全部内部结果阳性，也必须保留“无独立外部验证、无锁定签名、无临床阈值”的限制。

不得根据名义P值、单次最佳折、单个高频基因或漂亮的校准曲线宣称诊断价值。

## 十八、Figure 7审阅版设计

生成新的审阅版而不覆盖正式图：

- `03_figures/Figure7_enhanced_review.pdf`；
- `03_figures/Figure7_enhanced_review.png`，300 dpi且宽度至少2200像素；
- 所有底层绘图表格同时保存为CSV。

推荐五面板结构：

### A. 完整流程标签置换

- 展示1,000次aggregate OOF AUC空分布；
- 标出观察AUC；
- 显示经验P值和有效置换数；
- 横轴必须写明aggregate repeated-OOF AUC。

### B. 匹配随机基因空间特异性

- 展示1,000个表达/方差匹配77基因空间的AUC分布；
- 标出候选空间AUC；
- 显示特异性经验P值；
- 图注说明匹配只使用外层训练数据。

### C. 基因选择与方向稳定性

- 展示非零频率最高的15个基因；
- 点或条长度表示非零频率；
- 颜色表示主导系数方向；
- 同时显示非零系数中位数和IQR，必要时拆成上下对齐的两个小图；
- 不让基因名称互相遮挡。

### D. 候选模型与Ma比较器

- 在完全相同的5次外层重复中连接两个模型的AUC、AP和Brier；
- Brier面板明确标注“lower is better”；
- 旁边用简洁文字给出aggregate AUC差、paired DeLong P值、固定预测区间和refit-aware OOB敏感性区间；
- 不把5个重复当作5个独立样本进行显著性检验。

### E. ROC与探索性校准

- 左侧为两个模型aggregate OOF ROC；
- 右侧为4等频组校准图；
- 显示Brier score；
- 标注“internal cross-fitted predictions; n=40”；
- 不绘制或推荐最佳临床阈值。

图形要求：

- 统一使用项目现有Figure 7配色；
- 面板标签A–E大小一致；
- 图例符号必须与图内形状一致；
- 不叠字、不截字、不跨面板；
- 不在图中使用“validated diagnostic model”“biomarker signature”或“superior”一类措辞；
- PDF为矢量文本；
- PNG为至少300 dpi；
- 对最终PNG和PDF分别执行视觉QA并保存检查记录。

另外生成补充诊断图：

- refit-aware OOB候选AUC及AUC差分布；
- 每个外层拟合选择基因数分布；
- 最佳C分布；
- 匹配扩层距离和表达/方差匹配误差；
- 每个校准分箱的样本数。

## 十九、中文分析报告

生成：

`04_report/GSE123568_Figure7_enhanced_analysis_report_zh.md`

报告至少包括：

- 分析目的；
- 实际输入文件、路径和SHA-256；
- 样本数、分组、顺序和标签核对；
- 77候选基因和Ma 4基因核对；
- 训练折内探针折叠实现；
- 内层pipeline如何避免scaler泄漏；
- 外层、内层、C网格和随机种子；
- aggregate与repeat-level AUC、AP和Brier；
- paired DeLong和固定预测bootstrap；
- refit-aware OOB bootstrap结果及其限制；
- 标签置换结果；
- 匹配随机基因空间结果及匹配质量；
- 基因非零频率、方向一致性和系数分布；
- 探索性校准结果；
- 与旧版Figure 7数值比较；
- 所有失败、警告、无效bootstrap和模型不收敛情况；
- 结果判定；
- 不能支持的结论；
- 所有输出文件清单。

报告必须区分：

- 固定cross-fitted预测的不确定性；
- 包含模型重新拟合的OOB敏感性；
- 标签置换检验；
- 候选空间特异性检验；
- 外部验证。

前四项均不是外部验证。

## 二十、可复现性、软件环境和质量控制

使用项目指定运行时：

- Python：使用满足 `analysis/requirements_serum_enhanced.txt` 的 Python 3.12 环境；
- Rscript：使用安装了绘图依赖的 R 4.6.1 环境。

如需独立虚拟环境，创建在：

仓库外或被 `.gitignore` 排除的独立 Python 3.12 虚拟环境

不得修改系统Python环境。记录：

- Python版本；
- R版本；
- 所有直接依赖包版本；
- `pip freeze`；
- R `sessionInfo()`；
- CPU核数和实际jobs数；
- 开始、结束和总运行时间；
- 完整命令行；
- 脚本SHA-256；
- 输出文件SHA-256。

正式完成前运行：

- Python语法检查；
- 新增单元测试；
- `qa/check_serum_enhanced.py`；
- 项目级 `qa/check_repository.py`；
- CSV/JSON字段和行数检查；
- OOF覆盖检查；
- checkpoint去重检查；
- PDF/PNG存在性、尺寸和分辨率检查；
- Figure 7视觉检查。

任何关键QC失败时，不得把分析标记为完成。

## 二十一、需要输出的最低文件集合

正式结果至少包括：

```text
00_frozen_plan/analysis_plan_frozen.json
00_frozen_plan/input_manifest.csv
00_frozen_plan/sample_manifest.csv
00_frozen_plan/candidate_gene_manifest.csv
00_frozen_plan/ma_comparator_gene_manifest.csv
00_frozen_plan/probe_annotation_audit.csv
00_frozen_plan/file_checksums_sha256.txt
02_formal/probe_selection_by_outer_fold.csv
02_formal/outer_fold_assignments.csv
02_formal/oof_predictions_long.csv
02_formal/oof_predictions_aggregated.csv
02_formal/repeat_performance.csv
02_formal/aggregate_performance.csv
02_formal/selected_hyperparameters.csv
02_formal/model_comparison_fixed_oof.csv
02_formal/feature_coefficients_long.csv
02_formal/feature_stability_summary.csv
02_formal/features_per_outer_fit.csv
02_formal/label_permutation_null.csv
02_formal/label_permutation_summary.csv
02_formal/refit_oob_bootstrap_all.csv
02_formal/refit_oob_bootstrap_invalid.csv
02_formal/refit_oob_bootstrap_summary.csv
02_formal/matched_random_space_auc_null.csv
02_formal/matched_random_space_gene_map.csv.gz
02_formal/matched_random_space_qc.csv
02_formal/matched_random_space_summary.csv
02_formal/calibration_summary.csv
02_formal/calibration_bins.csv
02_formal/brier_comparison.csv
02_formal/legacy_vs_enhanced_comparison.csv
02_formal/serum_enhanced_summary.json
03_figures/Figure7_enhanced_review.pdf
03_figures/Figure7_enhanced_review.png
04_report/GSE123568_Figure7_enhanced_analysis_report_zh.md
04_report/run_log.txt
04_report/python_environment.txt
04_report/R_sessionInfo.txt
04_report/output_manifest.csv
```

若某项结果为空，也应输出带完整表头的空文件，并在报告中说明原因；不得省略失败或阴性结果文件。

## 二十二、停止条件

完成以下事项后停止：

1. 冻结并哈希分析方案；
2. 通过全部输入、样本、基因和泄漏QC；
3. 通过烟雾测试；
4. 完成正式嵌套CV、1,000次标签置换、1,000个有效refit-aware OOB bootstrap和1,000个匹配随机基因空间；
5. 完成固定预测配对统计、Brier和探索性校准；
6. 生成审阅版Figure 7及补充诊断图；
7. 完成中文报告和输出清单；
8. 运行代码、结果和图形QA；
9. 汇报新旧结果差异及是否改变原有审慎结论。

然后先向作者汇报并等待审核。不得自行：

- 替换 `figures/final/Figure7.*`；
- 覆盖 `results/figure_inputs`；
- 修改主稿、摘要、图注或补充材料；
- 创建新投稿包；
- 上传外部服务；
- 提交或推送Git；
- 根据结果追加新模型、放宽标准、增加置换次数或选择更有利的分析。
