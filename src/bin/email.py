#!/usr/bin/env python3

from __future__ import annotations

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.1"

"""
Send HLorf results by email.

Inputs
------
--bed     : final kept ORF predictions (BED12)
--tsv     : merged prediction table with header
--counts  : per-chunk counts (8 tab-separated columns, no header)
--yml     : merged versions.yml

The script builds a compact HTML summary (and plaintext fallback) and
delivers it through SMTP or mailx, also saving a copy under
<outdir>/reports/HLorf_email_summary.html.
"""

import argparse
import collections
import csv
import html
import os
import re
import sys
from statistics import mean
from typing import Any, Dict, List, Tuple

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR in sys.path:
    sys.path.remove(_SCRIPT_DIR)
    sys.path.append(_SCRIPT_DIR)


def e(s: Any) -> str:
    """HTML-escape a value; safe for None."""
    return html.escape("" if s is None else str(s))


# ---------------------------- Data loading ----------------------------- #

CountRow = collections.namedtuple(
    "CountRow",
    ["meta", "initial", "tai", "blast", "samba", "all", "unique", "kept"],
)


def parse_counts(path: str) -> List[CountRow]:
    rows: List[CountRow] = []
    with open(path, newline="") as fh:
        for raw in fh:
            raw = raw.strip()
            if not raw:
                continue
            cols = raw.split("\t")
            if len(cols) < 8:
                raise ValueError(f"Counts row has <8 columns: {raw}")
            cols = cols[:8]
            meta = cols[0]
            ints = [int(x) if x != "" else 0 for x in cols[1:]]
            rows.append(CountRow(meta, *ints))
    return rows


def parse_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row for row in reader]


def parse_bed_ids(path: str) -> List[str]:
    ids: List[str] = []
    with open(path, newline="") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            ids.append(cols[3])
    return ids


def parse_versions(path: str) -> List[Tuple[str, str, str]]:
    """
    Preserve possible duplicate blocks by parsing line-by-line.
    Returns list of unique (section, tool, version), keeping first occurrence.
    """
    triples: List[Tuple[str, str, str]] = []
    section = None
    section_re = re.compile(r'^"([^"]+)":\s*$')
    entry_re = re.compile(r"^\s+([^:]+):\s*(.+)\s*$")
    with open(path) as fh:
        for line in fh:
            m = section_re.match(line)
            if m:
                section = m.group(1)
                continue
            m = entry_re.match(line)
            if m and section:
                tool, ver = m.group(1).strip(), m.group(2).strip()
                triples.append((section, tool, ver))
    seen = set()
    unique: List[Tuple[str, str, str]] = []
    for t in triples:
        if t in seen:
            continue
        seen.add(t)
        unique.append(t)
    return unique


# ---------------------------- Statistics ------------------------------ #


def iqr_outliers(values: List[float]) -> Tuple[float, float, List[int]]:
    """Return (lower, upper, indices) where indices are positions flagged as outliers."""
    if not values:
        return (0.0, 0.0, [])
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    q1 = sorted_vals[int(0.25 * (n - 1))]
    q3 = sorted_vals[int(0.75 * (n - 1))]
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    idx = [i for i, v in enumerate(values) if v < lower or v > upper]
    return (lower, upper, idx)


def count_stats(rows: List[CountRow]) -> Dict[str, Any]:
    totals = collections.Counter()
    ratios: List[float] = []
    ratio_by_meta: Dict[str, float] = {}
    for r in rows:
        totals.update(
            {
                "initial": r.initial,
                "tai": r.tai,
                "blast": r.blast,
                "samba": r.samba,
                "all": r.all,
                "unique": r.unique,
                "kept": r.kept,
            }
        )
        ratio = (r.kept / r.initial) if r.initial else 0.0
        ratios.append(ratio)
        ratio_by_meta[r.meta] = ratio
    per_stage = dict(totals)

    lower, upper, out_idx = iqr_outliers(ratios)

    return {
        "per_stage": per_stage,
        "mean_retained": mean(ratios) if ratios else 0.0,
        "ratio_outliers": [
            {
                "meta": rows[i].meta,
                "kept": rows[i].kept,
                "initial": rows[i].initial,
                "ratio": ratios[i],
            }
            for i in out_idx
        ],
        "ratio_bounds": (lower, upper),
    }


def tsv_stats(rows: List[Dict[str, str]]) -> Dict[str, Any]:
    def f(name: str) -> List[float]:
        vals = []
        for r in rows:
            try:
                vals.append(float(r.get(name, "") or 0.0))
            except ValueError:
                vals.append(0.0)
        return vals

    def b(name: str) -> List[int]:
        vals = []
        for r in rows:
            v = r.get(name, "")
            if isinstance(v, str):
                v = v.strip()
            vals.append(1 if str(v) in {"1", "True", "true"} else 0)
        return vals

    tai_mean = f("tai_mean_score")
    rna_score = f("rna_score")
    blast_pid = f("blast_pid")
    blast_pct = f("blast_percentage_aligned")
    neglog = f("neg_log10_blast_evalue")
    canonical_start = b("has_canonical_start")
    canonical_stop = b("has_canonical_stop")
    # Derive class labels; fall back to probs when missing/empty
    class_labels: List[str] = []
    for r in rows:
        cls = (r.get("predicted_class") or "").strip()
        if not cls or cls.lower() == "predicted_class":
            try:
                pc = float(r.get("prob_coding", 0) or 0)
                pn = float(r.get("prob_noncoding", 0) or 0)

                if pc == pn == 0:
                    print(f"WARN: skipping row {r} due to missing probs")
            except ValueError:
                if not r["id"] == "id":
                    print(f"WARN: skipping row {r} due to missing probs")
                continue
            if pc == pn == 0:
                cls = "unknown"
            else:
                cls = "coding" if pc >= pn else "noncoding"

        if cls == "1":
            cls = "coding"
        elif cls == "0":
            cls = "noncoding"

        class_labels.append(cls)

    def avg(xs: List[float]) -> float:
        return mean(xs) if xs else 0.0

    return {
        "mean_tai": avg(tai_mean),
        "mean_rna": avg(rna_score),
        "mean_blast_pid": avg(blast_pid),
        "mean_blast_pct": avg(blast_pct),
        "mean_neglog": avg(neglog),
        "canonical_start_rate": avg(canonical_start),
        "canonical_stop_rate": avg(canonical_stop),
        "class_counts": collections.Counter(class_labels),
    }


def bed_stats(ids: List[str]) -> Dict[str, Any]:
    base_counts = collections.Counter()
    for ident in ids:
        base = ident
        if "_ORF" in ident:
            base = ident.split("_ORF", 1)[0]
        m = re.search(r"(.*)\.p\d+$", ident)
        if m:
            base = m.group(1)
        base_counts[base] += 1
    total = len(ids)
    top = base_counts.most_common(5)
    return {"total": total, "top": top}


# ------------------------- Report rendering --------------------------- #


def render_summary(counts: Dict[str, Any]) -> str:
    kept = counts["per_stage"].get("kept", 0)
    initial = counts["per_stage"].get("initial", 0) or 1
    ratio = kept / initial
    return f"""
    <table style="border-collapse:collapse;margin-bottom:12px;">
      <tr>
        <td style="padding:6px 10px;border:1px solid #ddd;background:#f5f5f5;">Initial</td>
        <td style="padding:6px 10px;border:1px solid #ddd;">{kept + (initial - kept)}</td>
        <td style="padding:6px 10px;border:1px solid #ddd;background:#f5f5f5;">Kept</td>
        <td style="padding:6px 10px;border:1px solid #ddd;">{kept}</td>
        <td style="padding:6px 10px;border:1px solid #ddd;background:#f5f5f5;">Kept / Initial</td>
        <td style="padding:6px 10px;border:1px solid #ddd;">{ratio:.3f}</td>
      </tr>
    </table>
    """


def render_counts(counts: Dict[str, Any]) -> str:
    per = counts["per_stage"]
    outliers = counts["ratio_outliers"]
    lower, upper = counts["ratio_bounds"]
    labels = {
        "initial": "Initial amount of records",
        "tai": "TranslationAi predictions",
        "blast": "BLAST hits",
        "samba": "RNASamba predictions",
        "all": "All predictions",
        "unique": "Unique predictions (non-duplicated)",
        "kept": "Predictions kept",
    }
    rows = "".join(
        f"<tr><td style='border:1px solid #ddd;padding:4px 8px'>{e(labels.get(stage, stage))}</td>"
        f"<td style='border:1px solid #ddd;padding:4px 8px;text-align:right'>{val}</td></tr>"
        for stage, val in per.items()
    )
    if outliers:
        details = "; ".join(
            f"{o['meta']} kept={o['kept']} initial={o['initial']} ratio={o['ratio']:.3f}"
            for o in outliers
        )
        out_html = (
            f"<p style='margin:4px 0'>Ratio outliers (kept/initial outside IQR fence "
            f"[{lower:.3f}, {upper:.3f}]): {details}</p>"
        )
    else:
        out_html = "<p style='margin:4px 0'>Ratio outliers: none</p>"
    return f"""
    <h3 style="margin:10px 0 6px;">Counts (per stage totals)</h3>
    <table style="border-collapse:collapse">{rows}</table>
    {out_html}
    """


def render_tsv(stats: Dict[str, Any]) -> str:
    class_rows = "".join(
        f"<tr><td style='border:1px solid #ddd;padding:4px 8px'>{e(k)}</td>"
        f"<td style='border:1px solid #ddd;padding:4px 8px;text-align:right'>{v}</td></tr>"
        for k, v in stats["class_counts"].items()
    )

    return f"""
    <h3 style="margin:10px 0 6px;">Prediction table stats</h3>
    <p style="margin:4px 0">
       Mean tai_mean_score (translation confidence): {stats["mean_tai"]:.3f}<br/>
       Mean rna_score (RNA coding potential): {stats["mean_rna"]:.3f}<br/>
       Mean blast % aligned (alignment coverage): {stats["mean_blast_pct"]:.3f}<br/>
       Mean blast PID (alignment identity): {stats["mean_blast_pid"]:.3f}<br/>
       Mean -log10(evalue) (higher is better): {stats["mean_neglog"]:.3f}
    </p>
    <p style="margin:4px 0">
       Canonical start rate: {stats["canonical_start_rate"]:.3f} &nbsp;|&nbsp;
       Canonical stop rate: {stats["canonical_stop_rate"]:.3f}
    </p>
    """
    # <h3 style="margin:8px 0 4px;">Class distribution</h4>
    # <table style="border-collapse:collapse;margin-bottom:8px">{class_rows}</table>


def render_versions(triples: List[Tuple[str, str, str]]) -> str:
    rows = "".join(
        f"<tr>"
        f"<td style='border:1px solid #ddd;padding:4px 8px'>{e(sec)}</td>"
        f"<td style='border:1px solid #ddd;padding:4px 8px'>{e(tool)}</td>"
        f"<td style='border:1px solid #ddd;padding:4px 8px'>{e(ver)}</td>"
        f"</tr>"
        for sec, tool, ver in triples
    )
    return f"""
    <h3 style="margin:12px 0 6px;">Versions</h3>
    <table style="border-collapse:collapse">{rows}</table>
    """


def build_email_html(
    subject: str,
    counts_stats: Dict[str, Any],
    tsv_stats_obj: Dict[str, Any],
    versions: List[Tuple[str, str, str]],
) -> str:
    return f"""<!doctype html>
<html>
<body style="font-family:Arial,Helvetica,sans-serif;font-size:13px;color:#222;">
  <h2 style="margin:6px 0 12px;">{e(subject)}</h2>
  {render_summary(counts_stats)}
  {render_counts(counts_stats)}
  {render_tsv(tsv_stats_obj)}
  {render_versions(versions)}
  <p style="margin-top:16px;color:#555;">Sent automatically by HLorf pipeline.</p>
</body>
</html>"""


# ----------------------------- Email send ----------------------------- #


def to_plaintext(html_body: str) -> str:
    text = re.sub(r"<br\s*/?>", "\n", html_body)
    text = re.sub(r"</p>", "\n\n", text)
    text = re.sub(r"<[^>]+>", "", text)
    return html.unescape(text).strip()


def send_email(
    to_addr: str,
    subject: str,
    html_body: str,
    plaintext: bool,
    smtp_server: str,
    smtp_port: int,
    smtp_user: str,
    smtp_password: str,
    from_addr: str,
    smtp_security: str,
    use_mailx: bool,
) -> None:
    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText

    if use_mailx:
        cmd = ["mailx"]
        if from_addr:
            cmd.extend(["-r", from_addr])
        cmd.extend(
            [
                "-S",
                "mime=1",
                "-S",
                f"content-type=text/{'plain' if plaintext else 'html'}",
                "-S",
                "charset=UTF-8",
                "-s",
                subject,
                to_addr,
            ]
        )
        import subprocess

        body = to_plaintext(html_body) if plaintext else html_body
        subprocess.run(cmd, input=body, text=True, check=True)
        return

    body = to_plaintext(html_body) if plaintext else html_body
    msg = MIMEMultipart()
    msg["From"] = from_addr or smtp_user
    msg["To"] = to_addr
    msg["Subject"] = subject
    msg.attach(MIMEText(body, "plain" if plaintext else "html", "utf-8"))

    if smtp_security not in {"ssl", "tls"}:
        raise ValueError("smtp_security must be ssl or tls")

    connector = smtplib.SMTP_SSL if smtp_security == "ssl" else smtplib.SMTP
    with connector(smtp_server, smtp_port) as server:
        if smtp_security == "tls":
            server.starttls()
        if smtp_user:
            server.login(smtp_user, smtp_password)
        server.send_message(msg)


# ------------------------------- CLI ---------------------------------- #


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Send HLorf results by email.")
    ap.add_argument("--bed", required=True, help="Final kept predictions (BED12)")
    ap.add_argument("--tsv", required=True, help="Merged ORF table (tsv with header)")
    ap.add_argument("--counts", required=True, help="Per-chunk counts (tsv, 8 cols)")
    ap.add_argument("--yml", required=True, help="Merged versions.yml")

    ap.add_argument("--outdir", required=True, help="Pipeline outdir to store report")
    ap.add_argument("--run-name", default="", help="Run name for subject")
    ap.add_argument(
        "--status",
        choices=["success", "failed"],
        default="success",
        help="Pipeline status label",
    )
    ap.add_argument("--email", default="", help="Recipient email")
    ap.add_argument(
        "--email-on-fail",
        default="",
        help="Recipient when status=failed and --email empty",
    )
    ap.add_argument("--plaintext-email", action="store_true", help="Send plaintext")
    ap.add_argument("--use-mailx", action="store_true", help="Send via mailx")

    ap.add_argument("--smtp-server", default="smtp.gmail.com")
    ap.add_argument("--smtp-port", type=int, default=465)
    ap.add_argument("--smtp-user", default="")
    ap.add_argument("--smtp-password", default="")
    ap.add_argument("--from-addr", default="")
    ap.add_argument("--smtp-security", choices=["ssl", "tls"], default="ssl")
    return ap.parse_args()


def main() -> int:
    args = parse_args()

    to_addr = args.email
    if not to_addr and args.email_on_fail and args.status == "failed":
        to_addr = args.email_on_fail

    subject = f"[HLorf] {args.status.capitalize()}: {args.run_name or 'run'}"

    # Load data
    counts_rows = parse_counts(args.counts)
    tsv_rows = parse_tsv(args.tsv)
    bed_ids = parse_bed_ids(args.bed)
    versions = parse_versions(args.yml)

    counts_summary = count_stats(counts_rows)
    tsv_summary = tsv_stats(tsv_rows)
    email_html = build_email_html(subject, counts_summary, tsv_summary, versions)

    reports_dir = os.path.join(args.outdir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    report_path = os.path.join(reports_dir, "HLorf_email_summary.html")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(email_html)

    if to_addr:
        send_email(
            to_addr=to_addr,
            subject=subject,
            html_body=email_html,
            plaintext=args.plaintext_email,
            smtp_server=args.smtp_server,
            smtp_port=args.smtp_port,
            smtp_user=args.smtp_user,
            smtp_password=args.smtp_password,
            from_addr=args.from_addr,
            smtp_security=args.smtp_security,
            use_mailx=args.use_mailx,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
