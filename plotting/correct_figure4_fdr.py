#!/usr/bin/env python
"""Synchronize the Figure 4 cGAS-STING title with the versioned statistics.

The reviewer-updated composite was exported as a single PDF by the latest
analysis notebook. This deterministic post-processing step prevents a manually
transcribed FDR in that composite from diverging from the result table.
"""
from __future__ import annotations

import csv
from io import BytesIO
from pathlib import Path

import pypdfium2
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "figures" / "source" / "Figure4_reviewed_base.pdf"
STATS = ROOT / "results" / "figure_inputs" / "module_scores_liao_stats_v4.csv"
OUTPUT_PDF = ROOT / "figures" / "final" / "Figure4.pdf"
OUTPUT_PNG = ROOT / "figures" / "final" / "Figure4.png"


def cgas_fdr() -> float:
    with STATS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["module"] == "cGAS_STING":
                return float(row["Liao_KW_fdr"])
    raise RuntimeError("cGAS_STING was not found in the module-statistics table")


def main() -> None:
    reader = PdfReader(BASE)
    if len(reader.pages) != 1:
        raise RuntimeError("Figure4_reviewed_base.pdf must contain exactly one page")
    page = reader.pages[0]
    width = float(page.mediabox.width)
    height = float(page.mediabox.height)

    overlay_bytes = BytesIO()
    overlay = canvas.Canvas(overlay_bytes, pagesize=(width, height))
    overlay.setFillColorRGB(1, 1, 1)
    overlay.rect(820, 695, width - 820, 30, fill=1, stroke=0)
    overlay.setFillColorRGB(0.18, 0.18, 0.18)
    overlay.setFont("Helvetica", 8.2)
    overlay.drawCentredString(985, 706.5, f"cGAS-STING (Liao exact KW FDR = {cgas_fdr():.2f})")
    overlay.save()
    overlay_bytes.seek(0)
    page.merge_page(PdfReader(overlay_bytes).pages[0])

    writer = PdfWriter()
    writer.add_page(page)
    writer.add_metadata(
        {
            "/Title": "Figure 4",
            "/Subject": "Reviewer-corrected analysis figure with table-derived FDR",
        }
    )
    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PDF.open("wb") as handle:
        writer.write(handle)

    pdf = pypdfium2.PdfDocument(OUTPUT_PDF)
    bitmap = pdf[0].render(scale=300 / 72)
    bitmap.to_pil().convert("RGB").save(OUTPUT_PNG, dpi=(300, 300))
    print(f"Figure 4 synchronized to cGAS-STING FDR={cgas_fdr():.2f}.")


if __name__ == "__main__":
    main()
