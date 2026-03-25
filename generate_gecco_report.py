#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path

# ---- helpers ----

def html_escape(s):
    if s is None:
        return ""
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def read_gene_scores(path: Path):
    scores = {}
    with path.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pid = row.get("protein_id")
            if not pid:
                continue
            try:
                scores[pid] = float(row.get("average_p", ""))
            except ValueError:
                pass
    return scores


def parse_location(loc: str):
    strand = "+"
    loc = loc.strip()
    if loc.startswith("complement("):
        strand = "-"
        loc = loc[len("complement("):-1]
    if loc.startswith("join(") and loc.endswith(")"):
        loc = loc[len("join("):-1]
    nums = [int(n) for n in re.findall(r"\d+", loc)]
    if not nums:
        return 0, 0, strand
    return min(nums), max(nums), strand


def _assign_qualifier(current, key, value):
    if key == "locus_tag":
        current["locus_tag"] = value
    elif key == "gene":
        current["gene"] = value
    elif key == "function":
        current["function"] = value
    elif key == "product":
        current["product"] = value
    elif key == "note":
        current.setdefault("notes", []).append(value)
    elif key == "db_xref":
        current.setdefault("db_xref", []).append(value)
    elif key == "translation":
        current["translation"] = value


def parse_gbk(path: Path, gene_scores):
    length = None
    genes_raw = []
    in_features = False
    current = None
    pending_key = None
    pending_val = ""

    with path.open() as f:
        for line in f:
            if line.startswith("LOCUS"):
                m = re.search(r"\s(\d+)\s+bp", line)
                if m:
                    length = int(m.group(1))
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            if line.startswith("ORIGIN"):
                break
            if line.startswith("     CDS"):
                if pending_key and current is not None:
                    _assign_qualifier(current, pending_key, pending_val.strip())
                pending_key = None
                pending_val = ""

                loc = line[21:].strip()
                current = {
                    "loc": loc,
                    "locus_tag": None,
                    "gene": None,
                    "function": None,
                    "product": None,
                    "notes": [],
                    "db_xref": [],
                    "translation": "",
                }
                genes_raw.append(current)
                continue
            if current is None:
                continue

            s = line.rstrip("\n")
            if pending_key:
                piece = s.strip()
                if piece.endswith('"'):
                    pending_val += " " + piece[:-1]
                    _assign_qualifier(current, pending_key, pending_val.strip())
                    pending_key = None
                    pending_val = ""
                else:
                    pending_val += " " + piece
                continue

            s_strip = s.strip()
            if not s_strip.startswith("/"):
                continue
            if "=" not in s_strip:
                continue
            key, val = s_strip[1:].split("=", 1)
            key = key.strip()
            val = val.strip()
            if val.startswith('"'):
                val = val[1:]
                if val.endswith('"'):
                    val = val[:-1]
                    _assign_qualifier(current, key, val)
                else:
                    pending_key = key
                    pending_val = val
            else:
                _assign_qualifier(current, key, val)

    if pending_key and current is not None:
        _assign_qualifier(current, pending_key, pending_val.strip())

    genes = []
    for g in genes_raw:
        start, end, strand = parse_location(g["loc"])
        locus = g.get("locus_tag") or ""
        genes.append({
            "start": start,
            "end": end,
            "strand": strand,
            "locus_tag": locus,
            "gene": g.get("gene") or "",
            "function": g.get("function") or "",
            "product": g.get("product") or "",
            "notes": g.get("notes") or [],
            "db_xref": g.get("db_xref") or [],
            "avg_p": gene_scores.get(locus),
            "gff_id": "",
            "protein_seq": g.get("translation") or "",
        })

    return length, genes


def parse_gff3(path: Path, feature_type: str, id_attr: str):
    records = []
    if not path or not path.exists():
        return records
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != feature_type:
                continue
            attr_map = {}
            for item in attrs.split(";"):
                if not item:
                    continue
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_map[k] = v
            gff_id = attr_map.get(id_attr) or attr_map.get("ID") or attr_map.get("Name") or ""
            if "," in gff_id:
                gff_id = gff_id.split(",", 1)[0]
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            records.append({
                "seqid": seqid,
                "start": s,
                "end": e,
                "strand": strand,
                "id": gff_id,
            })
    return records


def read_expr_table(path: Path, sep: str, id_col: str, cond_col: str, log2fc_col: str, padj_col: str):
    expr = {}
    if not path or not path.exists():
        return expr
    with path.open() as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            gid = (row.get(id_col) or "").strip()
            if not gid:
                continue
            cond = (row.get(cond_col) or "").strip()
            if not cond:
                continue
            try:
                log2fc = float(row.get(log2fc_col, ""))
            except ValueError:
                log2fc = None
            try:
                padj = float(row.get(padj_col, ""))
            except ValueError:
                padj = None
            expr.setdefault(gid, {})
            # keep the most significant entry per condition (lowest padj)
            prev = expr[gid].get(cond)
            if prev is None or (padj is not None and (prev.get("padj") is None or padj < prev.get("padj"))):
                expr[gid][cond] = {"log2fc": log2fc, "padj": padj}
    return expr


def classify_expr(log2fc, padj, fc_th, padj_th):
    if log2fc is None or padj is None:
        return "na"
    if padj <= padj_th and log2fc >= fc_th:
        return "up"
    if padj <= padj_th and log2fc <= -fc_th:
        return "down"
    return "ns"


def lerp(a, b, t):
    return int(a + (b - a) * t)


def color_for(avg_p):
    if avg_p is None:
        return "#8b8b8b"
    # high contrast gradient: blue (low) -> amber (high), avoid red/green clash
    t = max(0.0, min(1.0, avg_p))
    r = lerp(37, 245, t)
    g = lerp(99, 158, t)
    b = lerp(235, 11, t)
    return f"#{r:02x}{g:02x}{b:02x}"


def svg_gene_track(length, genes):
    if not length or not genes:
        return "<div class=\"muted\">No CDS features found in GBK.</div>"

    width = 1000
    height = 64
    y = 20
    h = 18
    parts = [f"<svg class=\"gene-track\" viewBox=\"0 0 {width} {height}\" preserveAspectRatio=\"none\">"]
    parts.append(f"<rect x=\"0\" y=\"{y + h/2 - 1}\" width=\"{width}\" height=\"2\" fill=\"#cbd5e1\" />")

    for g in genes:
        start = g["start"]
        end = g["end"]
        if end <= start or start <= 0:
            continue
        x = (start / length) * width
        w = max(2.0, (end - start) / length * width)
        arrow = min(10.0, w * 0.35)
        color = color_for(g.get("avg_p"))
        title = f"{g.get('locus_tag','')}  {g.get('function','') or g.get('product','')}  avg_p={g.get('avg_p') if g.get('avg_p') is not None else 'NA'}"

        expr_json = json.dumps(g.get("expr", []), ensure_ascii=False)
        data_attrs = (
            f"data-gene=\"1\" "
            f"data-locus=\"{html_escape(g.get('locus_tag',''))}\" "
            f"data-gene-name=\"{html_escape(g.get('gene',''))}\" "
            f"data-function=\"{html_escape(g.get('function',''))}\" "
            f"data-product=\"{html_escape(g.get('product',''))}\" "
            f"data-start=\"{g.get('start')}\" "
            f"data-end=\"{g.get('end')}\" "
            f"data-strand=\"{g.get('strand')}\" "
            f"data-avgp=\"{g.get('avg_p') if g.get('avg_p') is not None else 'NA'}\" "
            f"data-notes=\"{html_escape('; '.join(g.get('notes') or []))}\" "
            f"data-dbx=\"{html_escape('; '.join(g.get('db_xref') or []))}\" "
            f"data-gffid=\"{html_escape(g.get('gff_id',''))}\" "
            f"data-protein=\"{html_escape(g.get('protein_seq',''))}\" "
            f"data-protlen=\"{len(g.get('protein_seq',''))}\" "
            f"data-expr=\"{html_escape(expr_json)}\""
        )

        if w < 8:
            parts.append(
                f"<rect x=\"{x:.2f}\" y=\"{y}\" width=\"{w:.2f}\" height=\"{h}\" fill=\"{color}\" opacity=\"0.9\" {data_attrs}><title>{title}</title></rect>"
            )
            continue

        if g["strand"] == "+":
            points = [
                (x, y),
                (x + w - arrow, y),
                (x + w, y + h / 2),
                (x + w - arrow, y + h),
                (x, y + h),
            ]
        else:
            points = [
                (x + arrow, y),
                (x + w, y),
                (x + w, y + h),
                (x + arrow, y + h),
                (x, y + h / 2),
            ]
        pts = " ".join(f"{px:.2f},{py:.2f}" for px, py in points)
        parts.append(
            f"<polygon points=\"{pts}\" fill=\"{color}\" opacity=\"0.9\" {data_attrs}><title>{title}</title></polygon>"
        )

    parts.append("</svg>")
    return "\n".join(parts)


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate an antiSMASH-like HTML report from GECCO output",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input-dir", type=Path, default=Path(__file__).resolve().parent,
                   help="GECCO 输出目录")
    p.add_argument("-c", "--clusters", type=Path, default=None,
                   help="clusters TSV 文件路径")
    p.add_argument("-g", "--genes", type=Path, default=None,
                   help="genes TSV 文件路径")
    p.add_argument("--gbk-dir", type=Path, default=None,
                   help="GBK 文件目录（默认与输入目录相同）")
    p.add_argument("--gbk-suffix", type=str, default=".gbk",
                   help="GBK 文件后缀")
    p.add_argument("--gff", type=Path, default=None,
                   help="GFF3 注释文件路径（可选）")
    p.add_argument("--gff-feature", type=str, default="CDS",
                   help="GFF3 里用于匹配的 feature 类型")
    p.add_argument("--gff-id-attr", type=str, default="ID",
                   help="GFF3 里作为基因 ID 的属性键名（如 ID/Name/gene_id）")
    p.add_argument("--gff-min-overlap", type=float, default=0.5,
                   help="GFF3 匹配的最小覆盖率（基于 GECCO 基因长度）")
    p.add_argument("--gff-ignore-strand", action="store_true", default=False,
                   help="匹配时忽略链方向（默认要求一致）")
    p.add_argument("--expr", type=Path, default=None,
                   help="差异表达结果表（长表格式，可选）")
    p.add_argument("--expr-sep", type=str, default="\\t",
                   help="表达表分隔符（默认TAB）")
    p.add_argument("--expr-id-col", type=str, default="gene_id",
                   help="表达表中基因ID列名")
    p.add_argument("--expr-cond-col", type=str, default="condition",
                   help="表达表中条件列名")
    p.add_argument("--expr-log2fc-col", type=str, default="log2fc",
                   help="表达表中log2FC列名")
    p.add_argument("--expr-padj-col", type=str, default="padj",
                   help="表达表中校正P值列名")
    p.add_argument("--expr-log2fc-threshold", type=float, default=1.0,
                   help="log2FC阈值（>=上调，<=-阈值为下调）")
    p.add_argument("--expr-padj-threshold", type=float, default=0.05,
                   help="显著性阈值（padj）")
    p.add_argument("--expr-id-source", type=str, default="auto",
                   choices=["auto", "locus_tag", "gff_id"],
                   help="表达表中的基因ID对应来源（auto优先gff_id）")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="输出 HTML 路径")
    p.add_argument("--title", type=str, default="GECCO 次生代谢簇可视化报告",
                   help="报告标题")
    p.add_argument("--only-type", action="append", default=None,
                   help="只保留指定类型（可多次使用）")
    p.add_argument("--min-avg-p", type=float, default=None,
                   help="只保留 average_p >= 阈值的簇")
    p.add_argument("--max-clusters", type=int, default=0,
                   help="最多输出多少个簇（0 表示不限制）")
    return p.parse_args()


def main():
    args = parse_args()

    base = args.input_dir.resolve()
    clusters_tsv = args.clusters or (base / "Cs_JDA12.genome.clusters.tsv")
    genes_tsv = args.genes or (base / "Cs_JDA12.genome.genes.tsv")
    gbk_dir = (args.gbk_dir or base).resolve()
    gff_path = args.gff
    output_html = args.output or (base / "gecco_report.html")

    if not clusters_tsv.exists():
        raise SystemExit(f"Missing clusters file: {clusters_tsv}")

    scores = read_gene_scores(genes_tsv) if genes_tsv.exists() else {}
    gff_records = parse_gff3(gff_path, args.gff_feature, args.gff_id_attr) if gff_path else []
    gff_by_seqid = {}
    for rec in gff_records:
        gff_by_seqid.setdefault(rec["seqid"], []).append(rec)

    expr_by_gene = read_expr_table(
        args.expr, args.expr_sep, args.expr_id_col, args.expr_cond_col,
        args.expr_log2fc_col, args.expr_padj_col
    ) if args.expr else {}
    expr_id_source = args.expr_id_source
    all_conditions = sorted({c for gene in expr_by_gene.values() for c in gene})

    clusters = []
    counts = {}
    only_types = set(args.only_type) if args.only_type else None

    with clusters_tsv.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ctype = row.get("type") or "Unknown"
            if only_types and ctype not in only_types:
                continue

            avg_p = None
            try:
                avg_p = float(row.get("average_p", ""))
            except ValueError:
                avg_p = None
            if args.min_avg_p is not None:
                if avg_p is None or avg_p < args.min_avg_p:
                    continue

            counts[ctype] = counts.get(ctype, 0) + 1
            cluster_id = row.get("cluster_id")
            gbk_path = gbk_dir / f"{cluster_id}{args.gbk_suffix}" if cluster_id else None
            length = None
            genes = []
            if gbk_path and gbk_path.exists():
                length, genes = parse_gbk(gbk_path, scores)
                # map gene positions to genome coordinates for GFF3 matching
                try:
                    cluster_start = int(row.get("start", "0"))
                except ValueError:
                    cluster_start = 0
                seqid = row.get("sequence_id") or ""
                if cluster_start > 0:
                    for g in genes:
                        g["genome_start"] = cluster_start + g["start"] - 1
                        g["genome_end"] = cluster_start + g["end"] - 1
                        g["seqid"] = seqid
                        best_id = ""
                        best_cov = 0.0
                        gene_len = max(1, g["genome_end"] - g["genome_start"] + 1)
                        for rec in gff_by_seqid.get(seqid, []):
                            if (not args.gff_ignore_strand) and rec["strand"] not in (".", g["strand"]):
                                continue
                            ov_start = max(g["genome_start"], rec["start"])
                            ov_end = min(g["genome_end"], rec["end"])
                            if ov_end < ov_start:
                                continue
                            overlap = ov_end - ov_start + 1
                            cov = overlap / gene_len
                            if cov > best_cov:
                                best_cov = cov
                                best_id = rec["id"]
                        g["gff_id"] = best_id if best_cov >= args.gff_min_overlap else ""
                        # expression mapping
                        if expr_id_source == "gff_id":
                            expr_id = g["gff_id"]
                        elif expr_id_source == "locus_tag":
                            expr_id = g.get("locus_tag", "")
                        else:
                            expr_id = g["gff_id"] or g.get("locus_tag", "")
                        g["expr_id"] = expr_id
                        expr_entries = expr_by_gene.get(expr_id, {})
                        expr_list = []
                        expr_sig = False
                        for cond, vals in sorted(expr_entries.items()):
                            status = classify_expr(
                                vals.get("log2fc"), vals.get("padj"),
                                args.expr_log2fc_threshold, args.expr_padj_threshold
                            )
                            if status in ("up", "down"):
                                expr_sig = True
                            expr_list.append({
                                "condition": cond,
                                "log2fc": vals.get("log2fc"),
                                "padj": vals.get("padj"),
                                "status": status,
                            })
                        g["expr"] = expr_list
                        g["expr_sig"] = expr_sig
            clusters.append({
                "sequence_id": row.get("sequence_id"),
                "cluster_id": cluster_id,
                "start": row.get("start"),
                "end": row.get("end"),
                "average_p": row.get("average_p"),
                "max_p": row.get("max_p"),
                "type": ctype,
                "alkaloid_probability": row.get("alkaloid_probability"),
                "nrp_probability": row.get("nrp_probability"),
                "polyketide_probability": row.get("polyketide_probability"),
                "ripp_probability": row.get("ripp_probability"),
                "saccharide_probability": row.get("saccharide_probability"),
                "terpene_probability": row.get("terpene_probability"),
                "genes": genes,
                "length": length,
            })

            if args.max_clusters and len(clusters) >= args.max_clusters:
                break

    total_clusters = len(clusters)
    unique_types = len(counts)
    max_count = max(counts.values()) if counts else 1

    summary_rows = []
    for t, c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        bar_w = (c / max_count) * 100
        summary_rows.append(
            f"<div class=\"summary-row\">"
            f"<div class=\"summary-type\">{t}</div>"
            f"<div class=\"summary-bar\"><div class=\"summary-fill\" style=\"width:{bar_w:.1f}%\"></div></div>"
            f"<div class=\"summary-count\">{c}</div>"
            f"</div>"
        )

    cluster_cards = []
    for c in clusters:
        probs = [
            ("Alkaloid", c.get("alkaloid_probability")),
            ("NRP", c.get("nrp_probability")),
            ("Polyketide", c.get("polyketide_probability")),
            ("RiPP", c.get("ripp_probability")),
            ("Saccharide", c.get("saccharide_probability")),
            ("Terpene", c.get("terpene_probability")),
        ]
        probs_html = "".join(
            f"<div class=\"prob\"><span>{name}</span><span>{val}</span></div>" for name, val in probs
        )

        expr_section = ""
        if all_conditions:
            cond_cols = len(all_conditions)
            expr_rows = []
            for g in c.get("genes", []):
                expr_list = g.get("expr", [])
                if not expr_list:
                    continue
                status_by_cond = {e["condition"]: e["status"] for e in expr_list}
                if not any(s in ("up", "down") for s in status_by_cond.values()):
                    continue
                gene_label = g.get("gff_id") or g.get("locus_tag") or "-"
                badges = []
                for cond in all_conditions:
                    status = status_by_cond.get(cond, "na")
                    badges.append(f"<span class=\"badge {status}\">{status.upper()}</span>")
                expr_rows.append(
                    f"<div class=\"expr-row\">"
                    f"<span class=\"expr-gene\">{html_escape(gene_label)}</span>"
                    f"<div class=\"expr-badges\" style=\"grid-template-columns: repeat({cond_cols}, minmax(0,1fr));\">{''.join(badges)}</div>"
                    f"</div>"
                )
            if expr_rows:
                header = (
                    f"<div class=\"expr-head\">"
                    f"<span class=\"expr-gene\">Gene</span>"
                    f"<div class=\"expr-badges\" style=\"grid-template-columns: repeat({cond_cols}, minmax(0,1fr));\">"
                    + "".join(f"<span class=\"cond\">{html_escape(cn)}</span>" for cn in all_conditions)
                    + "</div></div>"
                )
                expr_section = (
                    f"<div class=\"expr-summary\">"
                    f"<div class=\"expr-title\">Expression summary (DE)</div>"
                    f"{header}"
                    f"{''.join(expr_rows)}"
                    f"</div>"
                )
            else:
                expr_section = (
                    f"<div class=\"expr-summary\">"
                    f"<div class=\"expr-title\">Expression summary (DE)</div>"
                    f"<div class=\"muted\">未检测到显著上/下调基因（或表达表未覆盖该簇）。</div>"
                    f"</div>"
                )

        gene_svg = svg_gene_track(c.get("length"), c.get("genes"))

        cluster_cards.append(
            f"<section class=\"cluster-card\" id=\"{c.get('cluster_id')}\">"
            f"<div class=\"cluster-header\">"
            f"<h3>{c.get('cluster_id')}</h3>"
            f"<div class=\"cluster-type\">{c.get('type')}</div>"
            f"</div>"
            f"<div class=\"cluster-meta\">"
            f"<div><strong>Contig</strong> {c.get('sequence_id')}</div>"
            f"<div><strong>Region</strong> {c.get('start')}..{c.get('end')}</div>"
            f"<div><strong>Average_p</strong> {c.get('average_p')}</div>"
            f"<div><strong>Max_p</strong> {c.get('max_p')}</div>"
            f"</div>"
            f"<div class=\"probs\">{probs_html}</div>"
            f"<div class=\"gene-panel\">{gene_svg}</div>"
            f"{expr_section}"
            f"<div class=\"gene-info\">"
            f"  <div class=\"gene-info-title\">Gene details</div>"
            f"  <div class=\"gene-info-body muted\">点击上方基因箭头查看信息</div>"
            f"</div>"
            f"</section>"
        )

    filter_notes = []
    if only_types:
        filter_notes.append(f"types={', '.join(sorted(only_types))}")
    if args.min_avg_p is not None:
        filter_notes.append(f"min_avg_p={args.min_avg_p}")
    if args.max_clusters:
        filter_notes.append(f"max_clusters={args.max_clusters}")
    filter_note_html = "" if not filter_notes else f"<div class=\"sub\">筛选: {'; '.join(filter_notes)}</div>"

    html = f"""
<!doctype html>
<html lang=\"zh\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>GECCO BGC Report</title>
  <style>
    :root {{
      --bg: #0b1220;
      --panel: #10192c;
      --panel-2: #111f39;
      --text: #e5e7eb;
      --muted: #a0aec0;
      --accent: #22c55e;
      --border: #1f2a44;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Fira Sans", "Avenir Next", "Helvetica Neue", Arial, sans-serif;
      background: radial-gradient(1200px 700px at 20% -10%, #1f2a44 0%, #0b1220 60%), #0b1220;
      color: var(--text);
    }}
    header {{
      padding: 28px 28px 18px;
      border-bottom: 1px solid var(--border);
      background: linear-gradient(120deg, rgba(34,197,94,0.15), rgba(14,165,233,0.05));
    }}
    h1 {{ margin: 0 0 6px; font-size: 24px; }}
    .sub {{ color: var(--muted); font-size: 14px; }}
    main {{ padding: 24px 28px 48px; }}
    .summary {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 18px;
      margin-bottom: 22px;
    }}
    .summary-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
    .summary-kpis {{ display: grid; grid-template-columns: repeat(3, minmax(0,1fr)); gap: 10px; }}
    .kpi {{ background: var(--panel-2); padding: 12px; border-radius: 8px; border: 1px solid var(--border); }}
    .kpi .label {{ color: var(--muted); font-size: 12px; }}
    .kpi .value {{ font-size: 20px; font-weight: 600; }}
    .summary-rows {{ margin-top: 12px; }}
    .summary-row {{ display: grid; grid-template-columns: 180px 1fr 60px; gap: 10px; align-items: center; margin: 6px 0; }}
    .summary-type {{ color: var(--text); font-size: 13px; }}
    .summary-bar {{ background: #0f172a; border-radius: 999px; height: 10px; overflow: hidden; border: 1px solid var(--border); }}
    .summary-fill {{ height: 100%; background: linear-gradient(90deg, #22c55e, #38bdf8); }}
    .summary-count {{ text-align: right; color: var(--muted); font-size: 13px; }}

    .clusters {{ display: grid; gap: 16px; }}
    .cluster-card {{ background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 16px; }}
    .cluster-header {{ display: flex; justify-content: space-between; align-items: baseline; gap: 10px; }}
    .cluster-header h3 {{ margin: 0; font-size: 18px; }}
    .cluster-type {{ padding: 4px 10px; border-radius: 999px; background: #0ea5e9; color: #051827; font-weight: 600; font-size: 12px; }}
    .cluster-meta {{ margin: 10px 0 12px; display: grid; grid-template-columns: repeat(4, minmax(0,1fr)); gap: 8px; font-size: 13px; color: var(--muted); }}
    .probs {{ display: grid; grid-template-columns: repeat(6, minmax(0,1fr)); gap: 8px; margin-bottom: 10px; }}
    .prob {{ background: var(--panel-2); border: 1px solid var(--border); border-radius: 8px; padding: 6px 8px; font-size: 12px; color: var(--muted); display: flex; justify-content: space-between; }}
    .gene-panel {{ background: #0c1424; border-radius: 10px; padding: 8px; border: 1px solid var(--border); }}
    .gene-track {{ width: 100%; height: 64px; }}
    .gene-track polygon, .gene-track rect {{ cursor: pointer; }}
    .gene-track polygon:hover, .gene-track rect:hover {{ stroke: #ffffff; stroke-width: 1; }}
    .gene-info {{ margin-top: 10px; padding: 10px; border: 1px solid var(--border); border-radius: 10px; background: #0c1424; }}
    .gene-info-title {{ font-size: 12px; color: var(--muted); margin-bottom: 6px; }}
    .gene-info-body {{ display: grid; gap: 6px; }}
    .gene-row {{ display: grid; grid-template-columns: 120px 1fr; gap: 8px; font-size: 12px; }}
    .seq-block {{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace; font-size: 11px; white-space: pre-wrap; word-break: break-all; background: #0b1324; border: 1px solid var(--border); border-radius: 8px; padding: 6px; color: #cbd5f5; }}
    .expr-panel {{ display: grid; gap: 6px; margin-top: 6px; }}
    .expr-item {{ display: grid; grid-template-columns: 140px 1fr; gap: 8px; font-size: 12px; align-items: center; }}
    .badge {{ display: inline-flex; align-items: center; justify-content: center; padding: 2px 6px; border-radius: 999px; font-size: 11px; font-weight: 600; }}
    .badge.up {{ background: #ef4444; color: #1f0a0a; }}
    .badge.down {{ background: #22c55e; color: #052912; }}
    .badge.ns {{ background: #64748b; color: #0b1220; }}
    .badge.na {{ background: #334155; color: #e2e8f0; }}
    .expr-summary {{ margin-top: 10px; padding: 10px; border: 1px solid var(--border); border-radius: 10px; background: #0c1424; }}
    .expr-title {{ font-size: 12px; color: var(--muted); margin-bottom: 6px; }}
    .expr-head, .expr-row {{ display: grid; grid-template-columns: 160px 1fr; gap: 8px; align-items: center; margin-bottom: 4px; }}
    .expr-gene {{ font-size: 12px; color: #e2e8f0; }}
    .expr-badges {{ display: grid; gap: 6px; }}
    .cond {{ font-size: 11px; color: #94a3b8; text-align: center; }}
    .muted {{ color: var(--muted); font-size: 12px; padding: 6px; }}

    @media (max-width: 980px) {{
      .summary-grid {{ grid-template-columns: 1fr; }}
      .summary-kpis {{ grid-template-columns: 1fr 1fr; }}
      .cluster-meta {{ grid-template-columns: 1fr 1fr; }}
      .probs {{ grid-template-columns: repeat(3, minmax(0,1fr)); }}
    }}
    @media (max-width: 640px) {{
      main {{ padding: 18px; }}
      .summary-kpis {{ grid-template-columns: 1fr; }}
      .probs {{ grid-template-columns: 1fr 1fr; }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>{html_escape(args.title)}</h1>
    <div class=\"sub\">输出目录：{base}</div>
    {filter_note_html}
  </header>
  <main>
    <section class=\"summary\">
      <div class=\"summary-grid\">
        <div>
          <div class=\"summary-kpis\">
            <div class=\"kpi\"><div class=\"label\">BGC 总数</div><div class=\"value\">{total_clusters}</div></div>
            <div class=\"kpi\"><div class=\"label\">类型数量</div><div class=\"value\">{unique_types}</div></div>
            <div class=\"kpi\"><div class=\"label\">含 GBK 可视化</div><div class=\"value\">{sum(1 for c in clusters if c.get('genes'))}</div></div>
          </div>
          <div class=\"summary-rows\">
            {"".join(summary_rows)}
          </div>
        </div>
        <div>
          <div class=\"muted\">颜色说明</div>
          <div class=\"muted\">基因颜色表示平均置信度 (avg_p)，从蓝色(低)到琥珀色(高)。</div>
          <div class=\"muted\">表达模式使用红色(上调)与绿色(下调)，与基因颜色区分开。</div>
          <div class=\"muted\">点击基因箭头后会在下方展示基因信息。</div>
        </div>
      </div>
    </section>

    <section class=\"clusters\">
      {"".join(cluster_cards)}
    </section>
  </main>

  <script>
    function formatSeq(seq) {{
      if (!seq) return '';
      var out = [];
      for (var i = 0; i < seq.length; i += 60) {{
        out.push(seq.slice(i, i + 60));
      }}
      return out.join('\\n');
    }}
    document.addEventListener('click', function (e) {{
      var target = e.target.closest('[data-gene]');
      if (!target) return;
      var card = target.closest('.cluster-card');
      var panel = card ? card.querySelector('.gene-info-body') : null;
      if (!panel) return;

      var data = target.dataset;
      var expr = [];
      try {{
        expr = JSON.parse(data.expr || '[]');
      }} catch (e) {{
        expr = [];
      }}
      var rows = [
        ['Locus tag', data.locus || '-'],
        ['Gene', data.geneName || '-'],
        ['Function', data.function || '-'],
        ['Product', data.product || '-'],
        ['GFF3 ID', data.gffid || '-'],
        ['Location', (data.start || '-') + '..' + (data.end || '-') + ' (' + (data.strand || '-') + ')'],
        ['Avg_p', data.avgp || 'NA'],
        ['DB xref', data.dbx || '-'],
        ['Notes', data.notes || '-']
      ];
      var html = rows.map(function (r) {{
        return '<div class="gene-row"><span>' + r[0] + '</span><span>' + r[1] + '</span></div>';
      }}).join('');

      if (expr.length) {{
        var exprRows = expr.map(function (e) {{
          var status = e.status || 'na';
          var badge = '<span class="badge ' + status + '">' + status.toUpperCase() + '</span>';
          var detail = 'log2fc=' + (e.log2fc === null || e.log2fc === undefined ? 'NA' : e.log2fc) +
            ', padj=' + (e.padj === null || e.padj === undefined ? 'NA' : e.padj);
          return '<div class="expr-item"><span>' + e.condition + '</span><span>' + badge + ' ' + detail + '</span></div>';
        }}).join('');
        html += '<div class="expr-panel">' + exprRows + '</div>';
      }} else {{
        html += '<div class="muted">表达谱：未提供或未匹配到该基因</div>';
      }}

      var seq = data.protein || '';
      if (seq) {{
        html += '<div class="gene-row"><span>Protein length</span><span>' + (data.protlen || seq.length) + ' aa</span></div>';
        html += '<div class="seq-block">' + formatSeq(seq) + '</div>';
      }} else {{
        html += '<div class="muted">序列：GBK 中未找到 /translation</div>';
      }}

      panel.innerHTML = html;
    }});
  </script>
</body>
</html>
"""

    output_html.write_text(html)
    print(f"Wrote {output_html}")


if __name__ == "__main__":
    main()
