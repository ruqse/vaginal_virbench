# Archived: Graphviz version of Figure 1 (study design)

Archived 2026-09-22. Moved out of `results/figures/`.

| File | What it is |
|---|---|
| `fig1_study_design.dot` | Graphviz source. This is the only code that generated Figure 1; no R/Python script draws it, and it does not produce the current `results/figures/fig1.png`, which is hand-made. |
| `fig1.svg` | SVG rendered from the .dot (Graphviz 2.44.0). |
| `fig1.png.bak-20260922-194925` | Earlier PNG render, backed up before `results/figures/fig1.png` was replaced at 19:49. |

Re-render (the paths in the .dot header point to `docs/figures/`, which is out of date; use these):

    dot -Tsvg fig1_study_design.dot -o fig1.svg
    dot -Tpng -Gdpi=180 fig1_study_design.dot -o fig1.png

`results/figures/fig1.png` (current version) was not archived.
