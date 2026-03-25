import argparse

def get_args():
    parser = argparse.ArgumentParser(
        prog="gecco_pipeline.py",
        description="""
GECCO BGC analysis pipeline

功能：
1. 解析 GECCO 输出（clusters / genes）
2. 筛选高置信 BGC
3. 整合 DESeq2 差异表达分析（可选）
4. 识别候选毒力相关 BGC
5. 输出 Excel + 可视化图

适用于：真菌/细菌次级代谢簇分析
""",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ==============================
    # 必选参数
    # ==============================
    required = parser.add_argument_group("🔥 必选参数")

    required.add_argument(
        "--clusters",
        required=True,
        help="GECCO 输出的 clusters.tsv 文件（BGC区域信息）"
    )

    required.add_argument(
        "--genes",
        required=True,
        help="GECCO 输出的 genes.tsv 文件（基因级信息）"
    )

    # ==============================
    # 可选参数
    # ==============================
    optional = parser.add_argument_group("⚙️ 可选参数")

    optional.add_argument(
        "--deseq",
        default=None,
        help="DESeq2 差异表达结果文件（CSV/TSV），用于整合表达分析"
    )

    optional.add_argument(
        "--outdir",
        default="gecco_analysis_out",
        help="输出目录（默认: gecco_analysis_out）"
    )

    optional.add_argument(
        "--prob_threshold",
        type=float,
        default=0.7,
        help="BGC置信度筛选阈值（默认: 0.7）"
    )

    optional.add_argument(
        "--log2fc",
        type=float,
        default=1.0,
        help="DESeq2 log2FoldChange 阈值（默认: 1.0）"
    )

    optional.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="DESeq2 padj/FDR 阈值（默认: 0.05）"
    )

    optional.add_argument(
        "--top_n",
        type=int,
        default=20,
        help="输出候选BGC数量（默认: 20）"
    )

    return parser.parse_args()


args = get_args()
