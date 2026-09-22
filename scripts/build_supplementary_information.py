#!/usr/bin/env python3
"""Assemble the Supplementary Information document (markdown -> docx + pdf).

Sources, all of which stay authoritative; nothing is retyped here:
  docs/supplementary_methods.md    extended procedures + supplementary references
  docs/figure_table_captions.md    Figure S1..Sn and Table S1..Sn legends
  results/figures/fig_s<n>.png     the figure images

Output (results/supplementary/):
  Supplementary_Information.md     assembled source
  Supplementary_Information.docx   12 pt Times New Roman
  Supplementary_Information.pdf    12 pt Times New Roman (Nimbus Roman substituted
                                   at render time if Times is not installed; the
                                   two are metric-compatible)

Typography note: the .docx stores the literal font name "Times New Roman", so it
resolves natively in Word on the co-authors' machines regardless of what is
installed on the build host.

    python3 scripts/build_supplementary_information.py            # md + docx + pdf
    python3 scripts/build_supplementary_information.py --md-only  # just the markdown
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent       # repo root (scripts/..)
DOCS = ROOT / "docs"                                 # SI source text
FIGDIR = ROOT / "results" / "figures"
SI_DIR = ROOT / "results" / "supplementary"         # SI outputs
BUILD = SI_DIR / "_si_build"

# pandoc 3.1.13 was used for the published SI (e.g. `module load Pandoc/3.1.13`).
# Override with the PANDOC environment variable if it is not on PATH.
PANDOC = os.environ.get("PANDOC") or shutil.which("pandoc") or "pandoc"
WEASYPRINT = "weasyprint"

# A4 with 2.5 cm margins leaves 16 cm of text width.
TEXT_WIDTH_IN = 6.3
MAX_FIG_HEIGHT_IN = 8.2

# Kept as a constant so the SI builds without the manuscript text. Update it if
# the article title changes.
TITLE = "Benchmarking virus identification tools on vaginal metagenomes"
AUTHORS = ("Faruk Dube, Anna-Ursula Happel, Heather B. Jaspan, "
           "Luisa Warchavchik Hugerth")


# ------------------------------------------------------------------ parsing
def split_methods() -> tuple[str, str]:
    """Return (body, references) from supplementary_methods.md.

    The references block is lifted out so it can be placed at the very end of the
    SI, after the tables, rather than stranded in the middle of the document.
    """
    text = (DOCS / "supplementary_methods.md").read_text(encoding="utf-8")
    text = re.sub(r"\A# Supplementary Methods\n+", "", text)
    # Drop the build-time pointer line; it is meaningless to a reader.
    text = re.sub(r"\AExtended procedures for `docs/manuscript\.md`\.[^\n]*\n+", "", text)
    parts = re.split(r"\n## References\s*\n", text, maxsplit=1)
    return (parts[0].strip(), parts[1].strip() if len(parts) > 1 else "")


def legends(kind: str, section: str) -> dict[int, str]:
    """Extract {n: legend_text} for '**Figure S<n>. ...**' style blocks.

    `section` is the '## ' heading whose block is scanned, so supplementary
    legends are never confused with the main-figure legends in the same file.
    """
    cap = (DOCS / "figure_table_captions.md").read_text(encoding="utf-8")
    m = re.search(rf"^## {re.escape(section)}\s*$(.*?)(?=^## |\Z)", cap, re.M | re.S)
    if not m:
        sys.exit(f"section '## {section}' not found in figure_table_captions.md")
    out = {}
    for blk in re.finditer(rf"^\*\*{kind} S(\d+)\..*?(?=^\*\*{kind} S\d+\.|\Z)",
                           m.group(1), re.M | re.S):
        out[int(blk.group(1))] = blk.group(0).strip()
    return out


# A page break that survives BOTH writers. `\newpage` is LaTeX: with no LaTeX in
# the toolchain the docx and html writers each drop it silently, so the breaks
# vanish without any error and figures run together mid-page. Emitting one raw
# block per target format lets each writer take its own and ignore the other.
PAGEBREAK = (
    '```{=openxml}\n<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n```\n\n'
    '```{=html}\n<div style="page-break-after: always;"></div>\n```'
)

SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def normalise_typography(text: str) -> str:
    """Turn Unicode subscript digits into pandoc subscript markup.

    The source files write log₁₀ because it is readable and greppable, but no
    Times-family serif carries U+2081/U+2080, so the renderer silently falls back
    to DejaVu Sans for those two glyphs alone: a sans-serif "10" sitting inside a
    serif word. Converting to ~10~ renders a true subscript in the body font in
    both the .docx and the .pdf.
    """
    text = re.sub(r"[₀-₉]+",
                  lambda m: "~" + m.group(0).translate(SUBSCRIPTS) + "~", text)
    # The sources spell it both ways; render both as a true subscript.
    return re.sub(r"\blog10\b", "log~10~", text)


def fig_width(path: Path) -> float:
    """Fit to text width, but shrink further if that would overrun the page."""
    with Image.open(path) as im:
        w, h = im.size
    width = TEXT_WIDTH_IN
    if h / w * width > MAX_FIG_HEIGHT_IN:
        width = MAX_FIG_HEIGHT_IN * w / h
    return round(width, 2)


# ------------------------------------------------------------------ assembly
def build_markdown() -> Path:
    methods, refs = split_methods()
    figs = legends("Figure", "Supplementary figures")
    tabs = legends("Table", "Supplementary tables")

    out = [
        "# Supplementary Information", "",
        f"**{TITLE}**", "",
        AUTHORS, "",
        f"This document contains Supplementary Methods, Supplementary Figures S1 to "
        f"S{max(figs)}, and legends for Supplementary Tables S1 to S{max(tabs)}. "
        f"The table data are in the accompanying workbook "
        f"`supplementary_tables.xlsx`, one sheet per table.", "",
        PAGEBREAK, "",
        "# Supplementary Methods", "", methods, "",
        PAGEBREAK, "",
        "# Supplementary Figures", "",
    ]

    missing = []
    for n in sorted(figs):
        img = FIGDIR / f"fig_s{n}.png"
        if not img.exists():
            missing.append(n)
            continue
        # Relative to results/supplementary/, where the assembled markdown lives: the
        # public .md must not carry build-host paths. pandoc resolves it via --resource-path.
        rel = Path(os.path.relpath(img, SI_DIR)).as_posix()
        # Alt text names the figure. implicit_figures is disabled in every pandoc
        # call, so the alt text is not turned into an extra visible caption.
        out += [f"![Figure S{n}]({rel}){{width={fig_width(img)}in}}", "", figs[n], "", PAGEBREAK, ""]
    if missing:
        sys.exit(f"missing figure image(s) for: {missing}")

    out += ["# Supplementary Tables", "",
            "Table data are in the accompanying workbook "
            "`supplementary_tables.xlsx`; each table below is one sheet, named "
            "by its table number.", ""]
    for n in sorted(tabs):
        out += [tabs[n], ""]

    if refs:
        out += [PAGEBREAK, "", "# Supplementary References", "", refs, ""]

    SI_DIR.mkdir(parents=True, exist_ok=True)
    md = SI_DIR / "Supplementary_Information.md"
    md.write_text(normalise_typography("\n".join(out)), encoding="utf-8")
    print(f"  markdown : {md.relative_to(ROOT)} "
          f"({len(figs)} figures, {len(tabs)} table legends)")
    return md


# ------------------------------------------------------------------ styling
def times_reference_docx() -> Path:
    """Pandoc's default reference.docx, restyled to 12 pt Times New Roman.

    Patched by editing styles.xml directly: python-docx is not installed on this
    host, and the edit is a handful of attribute rewrites.
    """
    BUILD.mkdir(exist_ok=True)
    base, patched = BUILD / "reference-default.docx", BUILD / "reference-times.docx"
    with base.open("wb") as fh:
        subprocess.run([PANDOC, "--print-default-data-file", "reference.docx"],
                       stdout=fh, check=True)

    with zipfile.ZipFile(base) as zin:
        items = {n: zin.read(n) for n in zin.namelist()}

    xml = items["word/styles.xml"].decode("utf-8")
    # Replace each <w:rFonts/> wholesale rather than rewriting font ATTRIBUTES.
    # Attribute-level rewriting is wrong twice over: the theme attribute is spelled
    # w:cstheme (not w:csTheme) so a cased alternation misses it and leaves a live
    # theme reference that outranks any explicit name; and a blanket
    # w:eastAsia="..." rewrite also hits <w:lang w:eastAsia="en-US"/>, silently
    # setting the document language to "Times New Roman".
    xml = re.sub(r"<w:rFonts\b[^>]*/>",
                 '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" '
                 'w:cs="Times New Roman" w:eastAsia="Times New Roman"/>', xml)
    # Size is deliberately NOT patched: pandoc's default reference.docx already
    # sets <w:sz w:val="24"/> (24 half-points == 12 pt) in docDefaults, and
    # re-inserting it by hand appended it after <w:lang/>, which violates the
    # OOXML child-element order for w:rPr.
    items["word/styles.xml"] = xml.encode("utf-8")

    with zipfile.ZipFile(patched, "w", zipfile.ZIP_DEFLATED) as zout:
        for name, data in items.items():
            zout.writestr(name, data)
    return patched


CSS = """
@page { size: A4; margin: 2.5cm; }
body { font-family: "Times New Roman", "Nimbus Roman", "Liberation Serif", serif;
       font-size: 12pt; line-height: 1.5; text-align: left; }
/* No page-break rule on h1: breaks are emitted explicitly in the markdown so
   the .docx and .pdf paginate identically. A CSS rule here would add a second
   break next to each explicit one and leave blank pages in the PDF only. */
h1 { font-size: 14pt; }
h2 { font-size: 12pt; }
h1, h2, h3 { font-weight: bold; }
img { max-width: 100%; height: auto; display: block; margin: 0 auto 0.6em auto; }
p { margin: 0 0 0.6em 0; orphans: 2; widows: 2; }
code { font-family: "Nimbus Mono PS", monospace; font-size: 11pt; }
"""


PANDOC_FROM = ["-f", "markdown-implicit_figures"]

# Messages that mean an image did not make it into the output. Both tools report
# these as warnings and still exit 0, so they are promoted to errors here.
UNRESOLVED = re.compile(r"Could not fetch resource|Failed to load image|"
                        r"Failed to load|not found|No such file", re.I)


def run_checked(cmd: list[str], what: str) -> None:
    """Run a build step; fail on a non-zero exit or on any unresolved-image message."""
    res = subprocess.run(cmd, capture_output=True, text=True)
    msgs = (res.stdout + res.stderr).strip()
    if msgs:
        print(msgs)
    if res.returncode != 0:
        sys.exit(f"{what} failed (exit {res.returncode})")
    bad = [l for l in msgs.splitlines() if UNRESOLVED.search(l)]
    if bad:
        sys.exit(f"{what}: unresolved resource(s):\n  " + "\n  ".join(bad))


def count_pdf_images(path: Path) -> tuple[int, int]:
    """Return (primary images, transparency masks) embedded in the PDF."""
    from pypdf import PdfReader
    from pypdf.generic import IndirectObject

    primary, masks, seen = set(), set(), set()

    def key(ref):
        return ref.idnum if isinstance(ref, IndirectObject) else id(ref)

    def walk(res):
        if not res or "/XObject" not in res:
            return
        for ref in res["/XObject"].values():
            if key(ref) in seen:
                continue
            seen.add(key(ref))
            obj = ref.get_object()
            if obj.get("/Subtype") == "/Image":
                primary.add(key(ref))
                if obj.get("/SMask") is not None:
                    masks.add(key(obj["/SMask"]))
            elif obj.get("/Subtype") == "/Form":
                walk(obj.get("/Resources"))

    for page in PdfReader(str(path)).pages:
        walk(page.get("/Resources"))
    return len(primary - masks), len(masks)


def build_docx(md: Path, n_figs: int) -> None:
    out = SI_DIR / "Supplementary_Information.docx"
    run_checked([PANDOC, *PANDOC_FROM, str(md), "-o", str(out),
                 f"--reference-doc={times_reference_docx()}",
                 "--resource-path", str(SI_DIR)], "docx build")
    with zipfile.ZipFile(out) as z:
        media = [n for n in z.namelist() if n.startswith("word/media/")]
    if len(media) != n_figs:
        sys.exit(f"docx embeds {len(media)} images, expected {n_figs}")
    print(f"  docx     : {out.relative_to(ROOT)} ({out.stat().st_size/1e6:.1f} MB, "
          f"{len(media)} images)")


def build_pdf(md: Path, n_figs: int) -> None:
    BUILD.mkdir(exist_ok=True)
    css, html = BUILD / "si.css", BUILD / "si.html"
    css.write_text(CSS, encoding="utf-8")
    run_checked([PANDOC, *PANDOC_FROM, str(md), "-t", "html5", "-s",
                 "-o", str(html), "--css", str(css),
                 "--resource-path", str(SI_DIR)], "html build")
    out = SI_DIR / "Supplementary_Information.pdf"
    # The HTML sits in _si_build/ but its image paths are relative to results/supplementary/.
    run_checked([WEASYPRINT, "--base-url", f"{SI_DIR}/", str(html), str(out)], "pdf build")
    n_img, n_mask = count_pdf_images(out)
    if n_img != n_figs:
        sys.exit(f"pdf embeds {n_img} primary images, expected {n_figs}")
    print(f"  pdf      : {out.relative_to(ROOT)} ({out.stat().st_size/1e6:.1f} MB, "
          f"{n_img} images + {n_mask} transparency masks)")


def main() -> int:
    md = build_markdown()
    if "--md-only" in sys.argv:
        return 0
    if not shutil.which(PANDOC):
        sys.exit(f"pandoc not found ({PANDOC}); try `module load Pandoc/3.1.13` or set PANDOC")
    n_figs = len(legends("Figure", "Supplementary figures"))
    build_docx(md, n_figs)
    if shutil.which(WEASYPRINT):
        build_pdf(md, n_figs)
    else:
        sys.exit("weasyprint not on PATH; the PDF cannot be built")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
