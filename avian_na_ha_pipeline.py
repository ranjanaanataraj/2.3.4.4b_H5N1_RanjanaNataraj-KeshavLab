#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
avian_na_ha_pipeline.py

Subcommands:
  na     → compute NA stalk lengths from aligned N1 FASTA
  ha     → count HA N-X-[ST] motifs (X != P)
  merge  → join HA/NA by EPI and (optionally) save contingency table
  plot   → make a bubble plot from a contingency table (or compute on the fly)

Notes:
- Requires EPI_ISL_* in FASTA headers. Dates (YYYY-MM-DD) are optional but recommended.
- Gaps ('-') tolerated; removed for length/GLS calculations.
"""

from Bio import SeqIO
import pandas as pd
import re, argparse, sys, json, os, time, logging
from pathlib import Path
from typing import Tuple, List, Optional

LOG = logging.getLogger("na-ha")

# ---------- FASTA & metadata ----------

def fasta_to_dataframe(fasta_file: str) -> pd.DataFrame:
    headers, sequences, lengths = [], [], []
    for record in SeqIO.parse(fasta_file, "fasta"):
        hdr = record.id if record.id else str(record.description).split()[0]
        seq = str(record.seq)
        headers.append(hdr)
        sequences.append(seq)
        lengths.append(len(seq))
    return pd.DataFrame({"Header": headers, "Sequence": sequences, "Length": lengths})

_EPI_RE = re.compile(r'EPI_ISL_\d+', re.IGNORECASE)
_DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}')

def extract_epi_and_date(header: str) -> Tuple[Optional[str], Optional[str]]:
    epi = (_EPI_RE.search(header).group(0)) if _EPI_RE.search(header) else None
    dt  = (_DATE_RE.search(header).group(0)) if _DATE_RE.search(header) else None
    return epi, dt

def annotate_epi_and_date(df: pd.DataFrame, label: str = "") -> pd.DataFrame:
    epis, dates = [], []
    for h in df["Header"]:
        epi, dt = extract_epi_and_date(h)
        if epi is None:
            LOG.warning("[%s] Missing EPI in header; dropped: %s", label, h)
        epis.append(epi); dates.append(dt)
    out = df.copy()
    out["EPI"] = epis; out["Date"] = dates
    out = out.dropna(subset=["EPI"])
    return out

def filter_by_date(df: pd.DataFrame, min_date: Optional[str], max_date: Optional[str]) -> pd.DataFrame:
    if "Date" not in df.columns:
        return df
    out = df.copy()
    if min_date:
        out = out[(out["Date"].isna()) | (out["Date"] >= min_date)]
    if max_date:
        out = out[(out["Date"].isna()) | (out["Date"] <= max_date)]
    return out

# ---------- NA stalk ----------

def find_stalk_bounds(reference_seq: str, begin_regex: str, end_regex: str,
                      start_offset: int = 0, end_offset: int = 0) -> Tuple[int,int]:
    m_begin = re.search(begin_regex, reference_seq, re.IGNORECASE)
    m_end   = re.search(end_regex, reference_seq, re.IGNORECASE)
    if not m_begin or not m_end:
        raise ValueError("Begin and/or end motifs not located in reference. Adjust regex/offsets.")
    start_idx = m_begin.end() + start_offset
    end_idx   = m_end.start() + end_offset
    if start_idx >= end_idx:
        raise ValueError("Computed stalk start >= end. Check regex/offsets.")
    return start_idx, end_idx

def extract_stalk_lengths(df: pd.DataFrame, stalk_start: int, stalk_end: int) -> pd.DataFrame:
    out = df.copy()
    out["Stalk_length"] = [
        seq[stalk_start:stalk_end].replace("-", "").__len__() for seq in out["Sequence"]
    ]
    return out

# ---------- HA GLS ----------

_GLS_RE = re.compile(r'[Nn][^Pp][SsTt]')

def count_gls_motifs(sequences: List[str]) -> List[int]:
    counts = []
    for seq in sequences:
        seq = seq.replace("-", "")
        counts.append(len(_GLS_RE.findall(seq)))
    return counts

# ---------- Utilities ----------

def write_audit_log(path: str, args: argparse.Namespace):
    meta = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "cmd": args.cmd,
        "args": vars(args),
    }
    with open(path, "w") as f:
        json.dump(meta, f, indent=2)

def setup_logging(level: str):
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

# ---------- Subcommands ----------

def run_na(a: argparse.Namespace):
    df = annotate_epi_and_date(fasta_to_dataframe(a.fasta), "NA")
    before = len(df)
    df = filter_by_date(df, a.min_date, a.max_date)
    LOG.info("[NA] rows: %d → %d after date filters", before, len(df))
    if df.empty: sys.exit("[NA] No valid NA sequences.")

    ref_seq = df["Sequence"].iloc[0]
    s_idx, e_idx = find_stalk_bounds(ref_seq, a.begin_regex, a.end_regex, a.start_offset, a.end_offset)
    if a.drop_first:
        df = df.iloc[1:].copy()
    df = extract_stalk_lengths(df, s_idx, e_idx)
    df[["EPI","Date","Stalk_length"]].to_csv(a.out, index=False)
    LOG.info("[NA] wrote %s (%d rows)", a.out, len(df))
    if a.audit: write_audit_log(a.audit, a)

def run_ha(a: argparse.Namespace):
    df = annotate_epi_and_date(fasta_to_dataframe(a.fasta), "HA")
    before = len(df)
    df = filter_by_date(df, a.min_date, a.max_date)
    LOG.info("[HA] rows: %d → %d after date filters", before, len(df))
    if df.empty: sys.exit("[HA] No valid HA sequences.")
    df = df.copy()
    df["GLS_count"] = count_gls_motifs(df["Sequence"].tolist())
    df[["EPI","Date","GLS_count"]].to_csv(a.out, index=False)
    LOG.info("[HA] wrote %s (%d rows)", a.out, len(df))
    if a.audit: write_audit_log(a.audit, a)

def run_merge(a: argparse.Namespace):
    ha = pd.read_csv(a.ha); na = pd.read_csv(a.na)
    miss = []
    for cols, dfname in [({"EPI","GLS_count"}, ha), ({"EPI","Stalk_length"}, na)]:
        need = cols - set(dfname.columns)
        if need: miss.append(str(need))
    if miss: sys.exit(f"[MERGE] Missing columns: {', '.join(miss)}")
    merged = ha.merge(na[["EPI","Stalk_length"]], on="EPI", how="inner")
    merged.to_csv(a.out, index=False)
    LOG.info("[MERGE] wrote %s (%d rows)", a.out, len(merged))
    if a.counts:
        counts = (merged.groupby(["Stalk_length","GLS_count"])
                  .size().reset_index(name="Frequency")
                  .sort_values(["Stalk_length","GLS_count"]))
        counts.to_csv(a.counts, index=False)
        LOG.info("[MERGE] wrote counts %s", a.counts)
    if a.audit: write_audit_log(a.audit, a)

def run_plot(a: argparse.Namespace):
    import matplotlib.pyplot as plt
    if a.counts_csv:
        df = pd.read_csv(a.counts_csv)
    else:
        merged = pd.read_csv(a.merged_csv)
        df = (merged.groupby(["Stalk_length","GLS_count"])
              .size().reset_index(name="Frequency"))
    if df.empty: sys.exit("[PLOT] No data to plot.")
    x = df["GLS_count"]; y = df["Stalk_length"]; s = df["Frequency"] * a.scale
    plt.figure()
    plt.scatter(x, y, s=s)
    plt.xlabel("HA GLS count")
    plt.ylabel("NA stalk length (aa)")
    plt.title("Stalk length × GLS count (bubble ~ frequency)")
    plt.grid(True, linestyle=":")
    plt.savefig(a.out, dpi=300, bbox_inches="tight")
    LOG.info("[PLOT] wrote %s", a.out)
    if a.audit: write_audit_log(a.audit, a)

# ---------- CLI ----------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="NA stalk-length + HA GLS pipeline")
    p.add_argument("--log-level", default="INFO", help="DEBUG|INFO|WARNING|ERROR")
    sub = p.add_subparsers(dest="cmd", required=True)

    pa = sub.add_parser("na", help="Compute NA stalk lengths")
    pa.add_argument("--fasta", required=True)
    pa.add_argument("--begin-regex", required=True)
    pa.add_argument("--end-regex", required=True)
    pa.add_argument("--start-offset", type=int, default=0)
    pa.add_argument("--end-offset", type=int, default=0)
    pa.add_argument("--drop-first", action="store_true")
    pa.add_argument("--min-date"); pa.add_argument("--max-date")
    pa.add_argument("--out", required=True)
    pa.add_argument("--audit", help="Write JSON audit log to this path")
    pa.set_defaults(func=run_na)

    ph = sub.add_parser("ha", help="Count HA glycosylation motifs")
    ph.add_argument("--fasta", required=True)
    ph.add_argument("--min-date"); ph.add_argument("--max-date")
    ph.add_argument("--out", required=True)
    ph.add_argument("--audit")
    ph.set_defaults(func=run_ha)

    pm = sub.add_parser("merge", help="Merge NA/HA CSVs; optional counts table")
    pm.add_argument("--ha", required=True)
    pm.add_argument("--na", required=True)
    pm.add_argument("--out", required=True)
    pm.add_argument("--counts")
    pm.add_argument("--audit")
    pm.set_defaults(func=run_merge)

    pp = sub.add_parser("plot", help="Bubble plot from counts or merged CSV")
    g = pp.add_mutually_exclusive_group(required=True)
    g.add_argument("--counts-csv", help="Precomputed contingency table")
    g.add_argument("--merged-csv", help="Merged HA-NA table")
    pp.add_argument("--scale", type=float, default=15.0, help="Bubble size scale factor")
    pp.add_argument("--out", required=True)
    pp.add_argument("--audit")
    pp.set_defaults(func=run_plot)

    return p

def main():
    args = build_parser().parse_args()
    setup_logging(args.log_level)
    args.func(args)

if __name__ == "__main__":
    main()
