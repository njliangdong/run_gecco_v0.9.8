# GECCO 使用说明（BGC预测工具）

## 📌 简介

GECCO（Genome-scale Enzyme Cluster Calling and cOntext）是一个基于条件随机场（CRF, Conditional Random Field）的次级代谢基因簇（BGC, Biosynthetic Gene Cluster）预测工具。

相比传统工具，GECCO 的特点是：

* 🚀 速度快
* 🧠 利用基因上下文信息
* 🧬 能发现非典型/未知 BGC

---

## 🧠 工作原理

GECCO 采用两步流程：

```
基因组 → 蛋白序列 → Pfam注释 → CRF模型 → BGC预测
```

核心思想：

* 使用 Pfam domain 作为特征
* 使用 CRF 模型识别基因是否属于 BGC
* 根据上下文动态确定 BGC 边界

---

## 🚀 安装（推荐）

```bash
mamba create -n gecco_env \
  python=3.9 \
  gecco=0.9.8 \
  biopython=1.79 \
  -c conda-forge -c bioconda \
  --strict-channel-priority -y
```

如果服务器较老（glibc < 2.28），建议使用：

* 旧版本 Python（3.9）
* 或 Singularity / Docker

---

## 📦 命令总览

```bash
gecco <command> [options]
```

### 可用命令：

| 命令       | 功能              |
| -------- | --------------- |
| run      | 一键完成 BGC 预测（推荐） |
| annotate | 蛋白功能注释（Pfam）    |
| predict  | 基于特征进行 BGC 预测   |
| convert  | 输出格式转换          |
| train    | 训练 CRF 模型       |
| cv       | 模型交叉验证          |

---

## 🔥 常用用法

### 1️⃣ 一键运行（推荐）

```bash
gecco run \
  --genome genome.fna \
  --proteins proteins.faa \
  --output gecco_out \
  -v
```

适用于：

* 大多数研究场景
* 快速 BGC 挖掘

---

### 2️⃣ 分步运行（高级用法）

#### Step 1：注释蛋白

```bash
gecco annotate \
  --proteins proteins.faa \
  -o features.tsv
```

#### Step 2：预测 BGC

```bash
gecco predict \
  --features features.tsv \
  -o clusters.tsv
```

适用于：

* HPC 并行计算
* 重复使用注释结果

---

### 3️⃣ 输出转换

```bash
gecco convert --input clusters.tsv
```

用途：

* 与其他工具兼容
* 下游可视化

---

## 🧬 输出结果说明

GECCO 输出通常包括：

* BGC 区域（cluster）
* 每个基因的概率评分
* Pfam 功能注释

---

## ⚖️ 与 antiSMASH 对比

| 特点      | GECCO   | antiSMASH |
| ------- | ------- | --------- |
| 方法      | CRF机器学习 | 规则+HMM    |
| 速度      | ⭐⭐⭐⭐⭐   | ⭐⭐        |
| 新型BGC发现 | ⭐⭐⭐⭐    | ⭐⭐        |
| 类型注释    | ⭐⭐      | ⭐⭐⭐⭐⭐     |

👉 推荐策略：

```
GECCO → 快速筛选
antiSMASH → 精细注释
```

---

## 🧪 科研应用建议

### 1️⃣ BGC 挖掘

重点关注：

* PKS / NRPS
* 毒素相关基因簇
* 黑色素（melanin）相关 cluster

---

### 2️⃣ 与转录组结合

```
BGC 区域 ∩ 差异表达基因（DEGs）
```

用于识别：

* 感染期激活的基因簇
* 潜伏 → 侵染转换相关 cluster

---

### 3️⃣ 新型 BGC 发现

GECCO 的优势：

* 不依赖已知规则
* 可发现 antiSMASH 未识别的 cluster

👉 适合：

* 新毒素
* 新代谢产物

---

## ⚠️ 常见问题

### ❌ 输入文件不匹配

* genome.fna 与 proteins.faa 必须来自同一注释

### ❌ contig 太短

* 会影响预测准确性

### ❌ 环境问题（glibc）

* 老服务器建议使用容器（Singularity）

---

## 📌 总结

GECCO 是一个基于机器学习的高速 BGC 预测工具，
能够利用基因上下文信息识别潜在的次级代谢基因簇，
特别适合大规模基因组和真菌研究。

---

## 🔗 项目地址

[https://github.com/zellerlab/GECCO](https://github.com/zellerlab/GECCO)

