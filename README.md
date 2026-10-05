# Kevin H. Wong's Open Lab Notebook

Live site: <https://kevinhwong1.github.io/KevinHWong_Notebook/>

Built with [Quarto](https://quarto.org) and published by GitHub Actions on every push to `master`.

## Layout

| Folder / file | What goes there |
|---|---|
| `notebook/` | Dated entries (`YYYY-MM-DD-title.qmd`): lab work, analyses, protocols, troubleshooting |
| `pipelines/` | Complete, maintained workflows; their scripts live in `pipelines/scripts/` |
| `protocols/` | Wet lab protocols (one `.qmd` per protocol, from `_templates/protocol.qmd`), reference protocols and manuals (PDF) |
| `projects.qmd` | Entries grouped by the `project:` field |
| `images/` | Figures used in posts (link as `/images/file.png`) |
| `_labels.yml` | The allowed types, organisms, methods and projects |
| `_templates/` | Templates used by `new_post.sh` |
| `_migration/` | One-time Jekyll → Quarto conversion script and label table |

## Everyday use

```bash
./new_post.sh                    # new notebook entry (menus for type, project, labels)
./new_post.sh troubleshooting    # error → cause → fix template
./new_post.sh pipeline           # new pipeline page + script stub
quarto preview                   # live preview in the browser
python3 check_labels.py          # catch label typos before pushing
git add -A && git commit -m "..." && git push
```

## Front matter

```yaml
---
title: "Pdam 31298: Cell Ranger → GeneExt → CellBender run"
date: 2026-05-07
type: "Analysis"                       # one of _labels.yml types
project: "Cnidarian Stem Cells"    # optional; must match projects.qmd
categories: ["Analysis", "Pocillopora damicornis", "scRNA-seq"]  # type + organisms + methods
tools: ["Cell Ranger", "GeneExt"]      # free text, shown under the title
protocols: ["dna-rna-coextraction-zymo"] # optional: wet-lab protocol page(s) used (file name in protocols/)
---
```

Old Jekyll URLs (`/Post-Title/`) redirect to the new pages through each post's `aliases:` field.

---

*Previously built with Jekyll Now (Barry Clark) and material from many open lab notebooks; thanks to all.*
