#!/usr/bin/env python3
"""Build Master's presentation PPTX matching the LVI thesis template style."""

from __future__ import annotations

from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Emu, Inches, Pt

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
PREVIEW = ROOT / "_template_preview"
DISCARD = ROOT.parent / "downstream" / "problematic_discard" / "plots"
OUT = ROOT / "CRC_Mutational_Signatures_Presentation.pptx"

# Template palette (sampled from Fatinikun LVI deck)
NAVY = RGBColor(27, 57, 107)
RED = RGBColor(226, 28, 60)
DARK = RGBColor(45, 45, 55)
GREY = RGBColor(90, 95, 105)
LIGHT_GREY = RGBColor(243, 244, 248)
WHITE = RGBColor(255, 255, 255)
LINE = RGBColor(220, 223, 230)
TABLE_HEADER = RGBColor(27, 57, 107)

SECTIONS = [
    "Introduction",
    "Materials",
    "Methods",
    "Results & Discussion",
    "Conclusion",
]

# Slide size matches template PDF mediabox 959.76 x 540 pt ≈ 13.33" x 7.5"
SLIDE_W = Inches(13.33)
SLIDE_H = Inches(7.5)


def set_run(run, *, size=18, bold=False, italic=False, color=DARK, font="Calibri"):
    run.font.name = font
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color


def add_textbox(slide, left, top, width, height, text, *, size=18, bold=False, italic=False,
                color=DARK, align=PP_ALIGN.LEFT, font="Calibri", anchor=MSO_ANCHOR.TOP):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    try:
        tf.auto_size = None
    except Exception:
        pass
    tf.paragraphs[0].alignment = align
    # vertical anchor via bodyPr
    bodyPr = tf._txBody.bodyPr
    if anchor == MSO_ANCHOR.MIDDLE:
        bodyPr.set("anchor", "ctr")
    elif anchor == MSO_ANCHOR.BOTTOM:
        bodyPr.set("anchor", "b")
    else:
        bodyPr.set("anchor", "t")
    p = tf.paragraphs[0]
    run = p.add_run()
    run.text = text
    set_run(run, size=size, bold=bold, italic=italic, color=color, font=font)
    return box


def add_paragraph(tf, text, *, size=16, bold=False, italic=False, color=DARK,
                  align=PP_ALIGN.LEFT, space_before=6, space_after=4):
    p = tf.add_paragraph()
    p.alignment = align
    p.space_before = Pt(space_before)
    p.space_after = Pt(space_after)
    run = p.add_run()
    run.text = text
    set_run(run, size=size, bold=bold, italic=italic, color=color)
    return p


def rounded_rect(slide, left, top, width, height, fill=LIGHT_GREY, line=None):
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill
    if line is None:
        shape.line.fill.background()
    else:
        shape.line.color.rgb = line
        shape.line.width = Pt(1)
    # softer corners
    try:
        adj = shape.adjustments
        adj[0] = 0.08
    except Exception:
        pass
    return shape


def circle_badge(slide, left, top, size, number):
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, left, top, size, size)
    circ.fill.solid()
    circ.fill.fore_color.rgb = RED
    circ.line.fill.background()
    tf = circ.text_frame
    tf.word_wrap = False
    bodyPr = tf._txBody.bodyPr
    bodyPr.set("anchor", "ctr")
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    run = p.add_run()
    run.text = str(number)
    set_run(run, size=14, bold=True, color=WHITE)
    return circ


def add_nav(slide, active: str):
    """Breadcrumb: Introduction > Materials > ... with active in red."""
    parts = []
    for s in SECTIONS:
        parts.append((s, s == active))
    # build as single textbox with multiple runs
    box = slide.shapes.add_textbox(Inches(0.45), Inches(0.22), Inches(9.5), Inches(0.35))
    tf = box.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    first = True
    for name, is_active in parts:
        if not first:
            sep = p.add_run()
            sep.text = "  >  "
            set_run(sep, size=11, color=GREY)
        run = p.add_run()
        run.text = name
        set_run(run, size=11, bold=is_active, color=RED if is_active else GREY)
        first = False


def add_logos(slide, title_style=False):
    logo = PREVIEW / ("embed_p1_1_xref7.png" if title_style else "embed_p2_0_xref22.png")
    if not logo.exists():
        return
    if title_style:
        slide.shapes.add_picture(str(logo), Inches(7.6), Inches(0.25), height=Inches(0.55))
    else:
        slide.shapes.add_picture(str(logo), Inches(10.55), Inches(0.18), height=Inches(0.42))


def add_page_number(slide, n: int):
    add_textbox(
        slide, Inches(12.5), Inches(7.1), Inches(0.6), Inches(0.3),
        str(n), size=11, color=GREY, align=PP_ALIGN.RIGHT,
    )


def add_title_block(slide, title: str, subtitle: str | None = None):
    add_textbox(
        slide, Inches(0.45), Inches(0.55), Inches(12.2), Inches(0.55),
        title, size=32, bold=True, color=NAVY,
    )
    if subtitle:
        add_textbox(
            slide, Inches(0.45), Inches(1.1), Inches(12.2), Inches(0.4),
            subtitle, size=16, color=GREY,
        )


def new_content_slide(prs, active: str, page: int, title: str, subtitle: str | None = None):
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    add_nav(slide, active)
    add_logos(slide, title_style=False)
    add_title_block(slide, title, subtitle)
    add_page_number(slide, page)
    return slide


def numbered_cards(slide, items: list[str], top=Inches(1.7), left=Inches(0.45),
                   width=Inches(12.4), card_h=Inches(1.15), gap=Inches(0.18)):
    for i, text in enumerate(items):
        y = top + i * (card_h + gap)
        rounded_rect(slide, left, y, width, card_h, fill=LIGHT_GREY)
        circle_badge(slide, left + Inches(0.22), y + Inches(0.35), Inches(0.45), i + 1)
        add_textbox(
            slide, left + Inches(0.9), y + Inches(0.25), width - Inches(1.2), card_h - Inches(0.35),
            text, size=16, color=DARK, anchor=MSO_ANCHOR.MIDDLE,
        )


def stat_cards(slide, stats: list[tuple[str, str, str]], top=Inches(1.65)):
    n = len(stats)
    gap = Inches(0.25)
    total_w = Inches(12.4)
    card_w = Emu((total_w - gap * (n - 1)) // n)
    left0 = Inches(0.45)
    for i, (big, label, note) in enumerate(stats):
        left = left0 + i * (card_w + gap)
        rounded_rect(slide, left, top, card_w, Inches(1.55), fill=WHITE, line=LINE)
        add_textbox(
            slide, left + Inches(0.15), top + Inches(0.2), card_w - Inches(0.3), Inches(0.65),
            big, size=28, bold=True, color=NAVY, align=PP_ALIGN.CENTER,
        )
        add_textbox(
            slide, left + Inches(0.15), top + Inches(0.85), card_w - Inches(0.3), Inches(0.3),
            label, size=13, bold=True, color=DARK, align=PP_ALIGN.CENTER,
        )
        add_textbox(
            slide, left + Inches(0.15), top + Inches(1.15), card_w - Inches(0.3), Inches(0.3),
            note, size=11, color=GREY, align=PP_ALIGN.CENTER,
        )


def add_table(slide, rows: list[list[str]], left, top, width, row_h=Inches(0.38)):
    nrows = len(rows)
    ncols = len(rows[0])
    table_shape = slide.shapes.add_table(nrows, ncols, left, top, width, row_h * nrows)
    table = table_shape.table
    col_w = Emu(width // ncols)
    for c in range(ncols):
        table.columns[c].width = col_w
    for r, row in enumerate(rows):
        for c, val in enumerate(row):
            cell = table.cell(r, c)
            cell.text = ""
            p = cell.text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if c == 0 or c == ncols - 1 else PP_ALIGN.CENTER
            run = p.add_run()
            run.text = val
            is_header = r == 0
            set_run(
                run,
                size=11 if not is_header else 11,
                bold=is_header or c == 0,
                color=WHITE if is_header else DARK,
            )
            cell.fill.solid()
            cell.fill.fore_color.rgb = TABLE_HEADER if is_header else (WHITE if r % 2 == 0 else LIGHT_GREY)
            # vertical center
            cell.text_frame._txBody.bodyPr.set("anchor", "ctr")
    return table_shape


def add_picture_fit(slide, path: Path, left, top, max_w, max_h):
    if not path.exists():
        add_textbox(slide, left, top, max_w, Inches(0.4), f"[Missing: {path.name}]", size=12, color=RED)
        return None
    pic = slide.shapes.add_picture(str(path), left, top)
    # scale to fit
    scale = min(max_w / pic.width, max_h / pic.height)
    pic.width = Emu(int(pic.width * scale))
    pic.height = Emu(int(pic.height * scale))
    return pic


def build():
    prs = Presentation()
    prs.slide_width = SLIDE_W
    prs.slide_height = SLIDE_H
    page = 0

    # ----- 1 Title -----
    page += 1
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    deco = PREVIEW / "embed_p1_0_xref6.png"
    if deco.exists():
        # large faded decoration on left (as in template)
        pic = slide.shapes.add_picture(str(deco), Inches(-0.8), Inches(0.8), height=Inches(6.2))
        # reduce opacity via alpha not easily available; leave as-is (template image already soft)
    add_logos(slide, title_style=True)
    add_textbox(
        slide, Inches(1.8), Inches(2.0), Inches(9.7), Inches(1.6),
        "Algorithm, Reference Catalogue, and NGS Scope\nShape Mutational-Signature Readouts in CRC",
        size=28, bold=True, color=NAVY, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.8), Inches(3.7), Inches(9.7), Inches(0.4),
        "Independent replication and expansion of Battuello et al. (2024)",
        size=16, italic=True, color=GREY, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.8), Inches(4.3), Inches(9.7), Inches(0.35),
        "Heritage Fatinikun",
        size=18, italic=True, color=DARK, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.8), Inches(4.7), Inches(9.7), Inches(0.55),
        "Institute of Tissue Medicine and Pathology, University of Bern\n"
        "Replication of Briefings in Bioinformatics 25:bbae249",
        size=13, color=GREY, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.8), Inches(5.9), Inches(9.7), Inches(0.55),
        "Master's / Internship Presentation\n20 October 2026",
        size=14, bold=True, color=RED, align=PP_ALIGN.CENTER,
    )
    add_page_number(slide, page)

    # ----- 2 Agenda -----
    page += 1
    slide = new_content_slide(prs, "Introduction", page, "Agenda", "What this talk covers")
    numbered_cards(slide, [
        "Motivation — why fitting algorithm, COSMIC catalogue, and NGS assay change CRC signature readouts.",
        "Materials & methods — cohorts, ten tools, six catalogues, dual criteria, downsampling, QC discard.",
        "Results — replication of Battuello et al. benchmarks, then expansions (tools, v3.4–v3.6, thresholds, discard).",
        "Conclusions — what is robust, what changes with workflow choice, and what we do not over-claim.",
    ], top=Inches(1.65), card_h=Inches(1.05))

    # ----- 3 Problem -----
    page += 1
    slide = new_content_slide(
        prs, "Introduction", page, "Scientific Problem",
        "Signature assignment is not a single standardised procedure",
    )
    numbered_cards(slide, [
        "CRC stratification uses MMR- and POLE-linked SBS signatures, but laboratories choose different fitting tools, COSMIC releases, and assays (panel / WES / WGS).",
        "Mathematical fit (cosine ≥ 0.9) and biological subtype contrast (ΔMMR, ΔPOLE) can move in opposite directions.",
        "Battuello et al. (2024) quantified these effects; we independently reimplement, expand, and audit paper-versus-replicate numbers.",
    ])

    # ----- 4 Framework -----
    page += 1
    slide = new_content_slide(
        prs, "Introduction", page, "Battuello Framework",
        "Three axes × two evaluation criteria",
    )
    # three axis cards
    axes = [
        ("Algorithm", "How the catalogue is fitted\n(MP, SPA, MuSiCal, …)"),
        ("Reference", "Which COSMIC catalogue\n(v2 → v3.6, CRC TS)"),
        ("NGS scope", "Panel vs WES vs WGS\n(coding-enriched vs genome-wide)"),
    ]
    for i, (t, d) in enumerate(axes):
        left = Inches(0.45) + i * Inches(4.15)
        rounded_rect(slide, left, Inches(1.7), Inches(3.95), Inches(2.0), fill=LIGHT_GREY)
        add_textbox(slide, left + Inches(0.2), Inches(1.9), Inches(3.55), Inches(0.45),
                    t, size=20, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
        add_textbox(slide, left + Inches(0.2), Inches(2.45), Inches(3.55), Inches(1.0),
                    d, size=14, color=DARK, align=PP_ALIGN.CENTER)
    rounded_rect(slide, Inches(0.45), Inches(4.0), Inches(6.05), Inches(2.4), fill=WHITE, line=LINE)
    add_textbox(slide, Inches(0.65), Inches(4.15), Inches(5.65), Inches(0.4),
                "Mathematical criterion", size=16, bold=True, color=NAVY)
    add_textbox(slide, Inches(0.65), Inches(4.6), Inches(5.65), Inches(1.5),
                "Cosine similarity between observed and reconstructed SBS96.\nAssigned if cosine ≥ 0.9.",
                size=14, color=DARK)
    rounded_rect(slide, Inches(6.75), Inches(4.0), Inches(6.05), Inches(2.4), fill=WHITE, line=LINE)
    add_textbox(slide, Inches(6.95), Inches(4.15), Inches(5.65), Inches(0.4),
                "Biological criterion", size=16, bold=True, color=NAVY)
    add_textbox(slide, Inches(6.95), Inches(4.6), Inches(5.65), Inches(1.5),
                "ΔMMR = median(MSI) − median(MSS)\nΔPOLE = median(POLE-MUT) − median(POLE-WT among MSS)",
                size=14, color=DARK)

    # ----- 5 Aims -----
    page += 1
    slide = new_content_slide(prs, "Introduction", page, "Project Aims", "Replication, expansion, and QC stress-test")
    numbered_cards(slide, [
        "Replicate key Battuello et al. benchmarks on CRC SBS96 matrices within stated numerical tolerance.",
        "Expand across additional fitting tools, COSMIC v3.4–v3.6, a CRC tissue-specific catalogue, and multi-catalogue mutation thresholds.",
        "Test whether discarding “problematic” samples (low cosine, low mutations, high flat load; Maura et al. 2019) systematically improves ΔMMR / ΔPOLE or clinical OS association.",
    ])

    # ----- 6 Cohorts -----
    page += 1
    slide = new_content_slide(
        prs, "Materials", page, "Cohorts",
        "Same CRC SBS96 matrices as Battuello et al. (GRCh38)",
    )
    stat_cards(slide, [
        ("230", "Preclinical WES", "primary algorithm grid"),
        ("230", "Preclinical TSO-500", "NGS-scope contrast"),
        ("63", "Preclinical WGS", "NGS-scope contrast"),
        ("152", "Clinical WES", "TCGA; thresholds + OS"),
    ])
    add_table(
        slide,
        [
            ["Cohort", "Assay", "n", "Role"],
            ["Preclinical cell lines", "WES", "230", "Algorithm × reference grid"],
            ["Preclinical cell lines", "TSO-500", "230", "Coding-enriched panel"],
            ["Preclinical cell lines", "WGS", "63", "Genome-wide contrast"],
            ["TCGA clinical CRC", "WES", "152", "Clinical Δ / survival / thresholds"],
        ],
        Inches(0.45), Inches(3.5), Inches(12.4), row_h=Inches(0.42),
    )
    add_textbox(
        slide, Inches(0.45), Inches(5.8), Inches(12.4), Inches(0.8),
        "Labels: MSI / MSS and POLE-MUT / POLE-WT from project metadata (paper supplementary tables). "
        "Master cosine table: 34,132 rows (10 tools × 6 catalogues × available datasets).",
        size=13, color=GREY,
    )

    # ----- 7 Tools & catalogues -----
    page += 1
    slide = new_content_slide(
        prs, "Materials", page, "Tools and Catalogues",
        "Paper core retained; five expansion tools; SignatureAnalyzer disabled",
    )
    rounded_rect(slide, Inches(0.45), Inches(1.65), Inches(6.05), Inches(4.8), fill=LIGHT_GREY)
    add_textbox(slide, Inches(0.7), Inches(1.85), Inches(5.55), Inches(0.4),
                "Fitting tools (n = 10)", size=18, bold=True, color=NAVY)
    add_textbox(slide, Inches(0.7), Inches(2.4), Inches(5.55), Inches(3.8),
                "Paper core\n• MutationalPatterns\n• deconstructSigs\n• signature.tools.lib\n"
                "• SigProfilerAssignment (SPA)\n• MuSiCal\n\n"
                "Expansion\n• sigLASSO  • sigfit  • sigminer\n• SigsPack  • YAPSA\n\n"
                "SignatureAnalyzer off — de novo NMF path, not catalogue-fitting equivalent.",
                size=14, color=DARK)
    rounded_rect(slide, Inches(6.75), Inches(1.65), Inches(6.05), Inches(4.8), fill=LIGHT_GREY)
    add_textbox(slide, Inches(7.0), Inches(1.85), Inches(5.55), Inches(0.4),
                "COSMIC references", size=18, bold=True, color=NAVY)
    add_textbox(slide, Inches(7.0), Inches(2.4), Inches(5.55), Inches(3.8),
                "• COSMIC v2\n• COSMIC v3.2\n• COSMIC v3.4\n• COSMIC v3.5\n• COSMIC v3.6\n"
                "• TS — CRC tissue-specific (24 SBS from paper GitHub list;\n"
                "  paper text cites 26 — disclosed)\n\n"
                "Shared contribution / reconstruction schema → cosine ≥ 0.9.",
                size=14, color=DARK)

    # ----- 8 Methods pipeline -----
    page += 1
    slide = new_content_slide(
        prs, "Methods", page, "Analysis Pipeline",
        "R + Python fitters → shared cosine schema → dual criteria → QC",
    )
    steps = [
        ("0", "SBS96\nmatrices"),
        ("1", "Fit\n10 tools"),
        ("2", "Cosine\n+ Δ scores"),
        ("3", "Downsample\nthresholds"),
        ("4", "Problematic\ndiscard"),
        ("5", "Paper-style\nfigures"),
    ]
    for i, (num, lab) in enumerate(steps):
        left = Inches(0.45) + i * Inches(2.15)
        rounded_rect(slide, left, Inches(2.0), Inches(1.95), Inches(2.0), fill=WHITE, line=LINE)
        add_textbox(slide, left, Inches(2.15), Inches(1.95), Inches(0.45),
                    num, size=22, bold=True, color=RED, align=PP_ALIGN.CENTER)
        add_textbox(slide, left + Inches(0.1), Inches(2.7), Inches(1.75), Inches(1.1),
                    lab, size=14, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
        if i < len(steps) - 1:
            add_textbox(slide, left + Inches(1.85), Inches(2.7), Inches(0.35), Inches(0.4),
                        "→", size=20, bold=True, color=GREY, align=PP_ALIGN.CENTER)
    add_textbox(
        slide, Inches(0.45), Inches(4.4), Inches(12.4), Inches(2.2),
        "Downsampling: multinomial SBS96 subsampling 5–95% × 5 replicates; MutationalPatterns refit for every COSMIC key "
        "(same scientific question as paper variant shuffling; thresholds approximate).\n\n"
        "Problematic discard (Maura-style): flag if cosine < 0.9 OR mutations below Fig. 6 threshold OR flat contribution ≥ 0.25; "
        "recompute ΔMMR / ΔPOLE and exploratory TCGA OS (n = 152; 36 events).\n\n"
        "Not reproduced (no local data): metanormal PBMC rematch; exon / intron / extragenic regional SBS96 pies.",
        size=14, color=DARK,
    )

    # ----- 9 Algorithm results -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Algorithm Choice",
        "Preclinical WES · COSMIC v2  (paper Fig. 4 analogue)",
    )
    add_picture_fit(slide, FIG / "Fig01_alg_cosine_by_tool_v2.png",
                    Inches(0.4), Inches(1.55), Inches(6.2), Inches(5.2))
    rounded_rect(slide, Inches(6.9), Inches(1.65), Inches(5.9), Inches(5.1), fill=LIGHT_GREY)
    add_textbox(slide, Inches(7.15), Inches(1.85), Inches(5.4), Inches(0.4),
                "Replication check", size=16, bold=True, color=NAVY)
    add_textbox(slide, Inches(7.15), Inches(2.35), Inches(5.4), Inches(4.1),
                "• Most tools: median cosine ≥ 0.9\n\n"
                "• SPA assigned 182/230 = 79.1%\n  (paper ≈ 79%)\n\n"
                "• Among paper tools, SPA highest ΔMMR\n  0.675 (MP 0.344; paper ≈ 0.34)\n\n"
                "• SPA ΔPOLE ≈ 0.70; MP ≈ 0.59\n  (paper ≈ 0.59)\n\n"
                "Expansion: sigLASSO ΔMMR 0.734\nexceeds SPA — qualify tool universe.",
                size=14, color=DARK)

    # ----- 10 Algorithm Δ -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Biological Separation by Tool",
        "ΔMMR and ΔPOLE on preclinical WES · v2",
    )
    add_picture_fit(slide, FIG / "Fig03_alg_delta_MMR_violin_v2.png",
                    Inches(0.35), Inches(1.55), Inches(6.3), Inches(5.2))
    add_picture_fit(slide, FIG / "Fig05_alg_delta_POLE_violin_v2.png",
                    Inches(6.7), Inches(1.55), Inches(6.2), Inches(5.2))

    # ----- 11 Reference -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Reference Catalogue",
        "MutationalPatterns · preclinical WES  (paper Fig. 5 analogue)",
    )
    add_table(
        slide,
        [
            ["Reference", "Paper ΔMMR", "This run", "Match"],
            ["COSMIC v2", "≈ 0.34", "0.344", "yes"],
            ["COSMIC v3.2", "≈ 0.26", "0.257", "yes"],
            ["Tissue-specific (TS)", "≈ 0.27", "0.272", "yes"],
        ],
        Inches(0.45), Inches(1.6), Inches(12.4), row_h=Inches(0.42),
    )
    add_picture_fit(slide, FIG / "Fig11_ref_cosine_v2_v32_TS.png",
                    Inches(0.4), Inches(3.5), Inches(6.0), Inches(3.4))
    add_textbox(
        slide, Inches(6.7), Inches(3.5), Inches(6.1), Inches(3.4),
        "Larger catalogues raise cosine but shrink ΔMMR — signal dilution / overfitting.\n\n"
        "Expansion to v3.4–v3.6 continues the pattern (MP ΔMMR ≈ 0.246 on v3.5/v3.6).\n\n"
        "Takeaway: mathematical and biological criteria trade off as the reference grows.",
        size=15, color=DARK,
    )

    # ----- 12 NGS -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "NGS Scope",
        "MutationalPatterns · COSMIC v3.2  (paper Fig. 3 analogue)",
    )
    add_table(
        slide,
        [
            ["Assay", "ΔMMR", "Flat in MSI (median)", "vs paper"],
            ["TSO-500", "0.323", "≈ 0", "highest (≈ 0.32)"],
            ["WES", "0.257", "≈ 0", "intermediate (≈ 0.26)"],
            ["WGS", "0.115", "0.269", "lowest; flat ↑ (≈ 0.11)"],
        ],
        Inches(0.45), Inches(1.55), Inches(12.4), row_h=Inches(0.4),
    )
    add_picture_fit(slide, FIG / "Fig10_ngs_signature_classes_MP_v3.2.png",
                    Inches(0.4), Inches(3.35), Inches(7.2), Inches(3.6))
    add_textbox(
        slide, Inches(7.8), Inches(3.5), Inches(5.0), Inches(3.3),
        "Ordering TSO > WES > WGS matches the paper.\n\n"
        "WGS cosine is high — not “technically worse” — but subtype contrast falls as flat signatures absorb MSI weight.\n\n"
        "Coding-enriched assays favour MMR/POLE stratification here.",
        size=14, color=DARK,
    )

    # ----- 13 Mutation threshold -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Mutation Burden Thresholds",
        "Multi-COSMIC expansion of paper Fig. 6 (metanormal panels blocked)",
    )
    add_picture_fit(slide, FIG / "Fig6C_preclinical_threshold_by_cosmic.png",
                    Inches(0.3), Inches(1.55), Inches(6.3), Inches(5.2))
    add_picture_fit(slide, FIG / "Fig6D_clinical_threshold_by_cosmic.png",
                    Inches(6.7), Inches(1.55), Inches(6.2), Inches(5.2))

    # ----- 14 Threshold numbers -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Thresholds Depend on Catalogue and Cohort",
        "Inferred mutations for cosine ≥ 0.9 (MutationalPatterns)",
    )
    add_table(
        slide,
        [
            ["Dataset", "COSMIC", "Inferred threshold", "Paper (v2)"],
            ["Preclinical", "v2", "≈ 351", "323"],
            ["Preclinical", "TS", "≈ 284", "—"],
            ["Preclinical", "v3.2–v3.6", "≈ 146", "—"],
            ["Clinical", "v2", "≈ 49", "64"],
            ["Clinical", "v3.2–v3.6", "≈ 17–21", "—"],
        ],
        Inches(0.45), Inches(1.65), Inches(12.4), row_h=Inches(0.45),
    )
    numbered_cards(
        slide,
        [
            "Richer catalogues reach cosine ≥ 0.9 at lower mutation counts — thresholds are not a sample-only property.",
            "Clinical WES needs far fewer mutations than unmatched preclinical cell lines (same direction as the paper).",
            "Paper metanormal rematch (323 → 145) was not reproduced — no local PBMC / metanormal data.",
        ],
        top=Inches(4.5), card_h=Inches(0.7), gap=Inches(0.12),
    )

    # ----- 15 Problematic discard -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Problematic-Sample Discard",
        "Maura-style QC: low cosine · low mutations · high flat (≥ 0.25)",
    )
    left_plot = DISCARD / "delta_MMR_change_preclinical.png"
    alt = DISCARD / "delta_MMR_all_vs_filtered_preclinical.png"
    add_picture_fit(slide, left_plot if left_plot.exists() else alt,
                    Inches(0.35), Inches(1.55), Inches(6.5), Inches(5.2))
    rounded_rect(slide, Inches(7.05), Inches(1.65), Inches(5.75), Inches(5.1), fill=LIGHT_GREY)
    add_textbox(slide, Inches(7.3), Inches(1.85), Inches(5.25), Inches(0.4),
                "Headline findings", size=16, bold=True, color=NAVY)
    add_textbox(slide, Inches(7.3), Inches(2.4), Inches(5.25), Inches(4.1),
                "• ≈ 10% of sample × tool × COSMIC rows flagged\n\n"
                "• Median preclinical ΔMMR change after discard: 0\n\n"
                "• ΔMMR increased in only 18 / 102 blocks (≈ 18%)\n\n"
                "• SPA on large catalogues: heavy flat flags "
                "(≈ 87% preclinical loss) with modest ΔMMR gains\n\n"
                "• TCGA OS (36 events): exploratory; no significant "
                "MMR high-vs-low association after filtering\n\n"
                "QC flags are useful diagnostics — not a default fix.",
                size=13, color=DARK)

    # ----- 16 Paper vs replicate -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Paper versus Replicate",
        "Key Battuello et al. benchmarks recovered within tight tolerance",
    )
    add_table(
        slide,
        [
            ["Metric", "Setting", "Paper", "This run"],
            ["ΔMMR", "MP · WES · v2", "≈ 0.34", "0.344"],
            ["ΔMMR", "MP · WES · v3.2", "≈ 0.26", "0.257"],
            ["ΔMMR", "MP · WES · TS", "≈ 0.27", "0.272"],
            ["Assigned ≥ 0.9", "SPA · WES · v2", "≈ 79%", "79.1%"],
            ["ΔMMR ranking", "SPA among paper tools", "highest", "0.675"],
            ["NGS ΔMMR order", "MP · v3.2", "TSO>WES>WGS", "0.323>0.257>0.115"],
            ["Flat in MSI", "MP · v3.2 · WGS", "elevated", "median 0.269"],
        ],
        Inches(0.45), Inches(1.55), Inches(12.4), row_h=Inches(0.48),
    )
    add_textbox(
        slide, Inches(0.45), Inches(5.6), Inches(12.4), Inches(1.2),
        "Full audit: results/thesis/paper_vs_replicate.tsv  ·  Manuscript: results/thesis/MANUSCRIPT.md",
        size=13, color=GREY,
    )

    # ----- 17 Synthesis -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "What Changes — and What Does Not",
        "Practical synthesis across axes",
    )
    numbered_cards(slide, [
        "Algorithm and reference are first-order drivers of both assignment rate and subtype Δ scores.",
        "Among paper tools, SPA maximises biological contrast at the cost of more unassigned samples; sigLASSO can exceed SPA when the tool set expands.",
        "Coding-enriched assays (TSO-500, WES) preserve MMR/POLE contrast better than WGS here, despite high WGS cosine.",
        "Discarding problematic samples does not systematically strengthen ΔMMR or OS association — report flags, do not auto-exclude by default.",
    ], card_h=Inches(1.05))

    # ----- 18 Limitations -----
    page += 1
    slide = new_content_slide(
        prs, "Results & Discussion", page, "Limitations",
        "What this work does not claim",
    )
    cards = [
        ("No metanormal", "Cannot claim paper Fig. 6C–D artefact rematch (323 → 145)."),
        ("No regional WGS pies", "Exon / intron / extragenic SBS96 unavailable; Fig10 class bars as proxy."),
        ("TS 24 vs “26”", "GitHub list used; ΔMMR still matches ≈ 0.27."),
        ("Survival power", "Only 36 OS events — exploratory Cox / log-rank only."),
        ("Downsampling method", "Multinomial SBS96 vs paper variant shuf — thresholds approximate."),
        ("Flat cut-off 0.25", "Aggressive for SPA on large catalogues; discard rates disclosed."),
    ]
    for i, (t, d) in enumerate(cards):
        col = i % 3
        row = i // 3
        left = Inches(0.45) + col * Inches(4.2)
        top = Inches(1.65) + row * Inches(2.4)
        rounded_rect(slide, left, top, Inches(4.0), Inches(2.15), fill=LIGHT_GREY)
        add_textbox(slide, left + Inches(0.2), top + Inches(0.25), Inches(3.6), Inches(0.45),
                    t, size=16, bold=True, color=NAVY)
        add_textbox(slide, left + Inches(0.2), top + Inches(0.8), Inches(3.6), Inches(1.1),
                    d, size=13, color=DARK)

    # ----- 19 Conclusions -----
    page += 1
    slide = new_content_slide(
        prs, "Conclusion", page, "Conclusions",
        "What we would defend",
    )
    numbered_cards(slide, [
        "Key Battuello benchmarks replicate: MP ΔMMR for v2/v3.2/TS, SPA assignment 79.1%, NGS order TSO > WES > WGS with WGS flat inflation.",
        "Expansions refine the story: newer catalogues further trade ΔMMR for cosine; mutation thresholds fall as catalogues grow; tool rankings shift (sigLASSO > SPA).",
        "Maura-style discard is useful QC reporting but not a blanket robustness fix in this CRC setting.",
        "Report tool + COSMIC + assay together. High cosine ≠ clinical subtype utility. Sequencing breadth does not replace biologically appropriate genomic scope.",
    ], card_h=Inches(1.05))

    # ----- 20 Thank you -----
    page += 1
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    if deco.exists():
        slide.shapes.add_picture(str(deco), Inches(-0.8), Inches(0.8), height=Inches(6.2))
    add_logos(slide, title_style=True)
    add_textbox(
        slide, Inches(1.5), Inches(2.6), Inches(10.3), Inches(0.7),
        "Thank you for your attention",
        size=36, bold=True, color=NAVY, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.5), Inches(3.5), Inches(10.3), Inches(0.45),
        "Feel free to ask any questions",
        size=18, italic=True, color=GREY, align=PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, Inches(1.5), Inches(4.5), Inches(10.3), Inches(1.2),
        "Heritage Fatinikun\n"
        "Institute of Tissue Medicine and Pathology · University of Bern\n"
        "Replication of Battuello et al., Brief Bioinform 2024 (bbae249)",
        size=14, color=DARK, align=PP_ALIGN.CENTER,
    )
    add_page_number(slide, page)

    # ----- 21 Backup -----
    page += 1
    slide = new_content_slide(
        prs, "Conclusion", page, "Backup — Exact Numbers",
        "For questions that ask “how many, exactly?”",
    )
    add_table(
        slide,
        [
            ["Item", "Value"],
            ["Master cosine rows", "34,132"],
            ["Tools enabled / catalogues", "10 / 6 (v2, v3.2–v3.6, TS)"],
            ["SPA assigned (preclinical WES · v2)", "182 / 230 (79.1%)"],
            ["MP ΔMMR v2 / v3.2 / TS", "0.344 / 0.257 / 0.272"],
            ["SPA / sigLASSO ΔMMR (v2)", "0.675 / 0.734"],
            ["NGS ΔMMR TSO / WES / WGS (MP · v3.2)", "0.323 / 0.257 / 0.115"],
            ["WGS flat in MSI (median)", "0.269"],
            ["Preclinical / clinical threshold (v2)", "≈ 351 / ≈ 49"],
            ["Rows flagged (any QC rule)", "≈ 10.1%"],
            ["ΔMMR↑ after discard (blocks)", "18 / 102"],
            ["TCGA OS events", "36 / 152"],
        ],
        Inches(0.45), Inches(1.5), Inches(12.4), row_h=Inches(0.38),
    )

    prs.save(OUT)
    print(f"Wrote {OUT} ({page} slides)")


if __name__ == "__main__":
    build()
