#!/usr/bin/env bash
# One-time Jekyll -> Quarto switch for this repo (already run on the `quarto` branch, kept for the record).
#   bash _migration/setup_quarto.sh
set -euo pipefail
cd "$(dirname "$0")/.."

# 1. Park the Jekyll machinery (Quarto ignores folders starting with "_")
mkdir -p _jekyll_archive
for f in _config.yml _layouts _includes _sass style.scss index.html about.md 404.md \
         archivebycategory.md archivebydate.md archivebytag.md CNAME; do
  if [[ -e "$f" ]]; then mv -n "$f" _jekyll_archive/; fi
done

# 2. Give the two Markdown protocols a title and keep their old URLs
add_fm() {  # file title  (edits in place)
  python3 - "$1" "$2" <<'PY'
import os, sys
f, title = sys.argv[1], sys.argv[2]
txt = open(f, encoding="utf-8").read()
if not txt.startswith("---"):
    lines = txt.split("\n")
    if lines and lines[0].startswith("# "):
        lines = lines[1:]
    slug = os.path.basename(f)[:-3]
    fm = f'---\ntitle: "{title}"\naliases:\n  - "/protocols/{slug}/"\n---\n'
    open(f, "w", encoding="utf-8").write(fm + "\n".join(lines))
PY
}
add_fm protocols/BS_DNA_Methylation_Analysis_Overview.md "Steps for QC and Analysis of Bisulfite Sequencing Data"
add_fm protocols/MeDIP_Protocol.md "Preparation of MeDIP-enriched, Bisulfite Converted Illumina Libraries"

# 3. Convert posts using _migration/labels.csv (generated from the rules if it doesn't exist yet)
[[ -e _migration/labels.csv ]] || python3 _migration/migrate_posts.py labels
python3 _migration/migrate_posts.py convert

# 4. Ignore Quarto build output
grep -qxF '/.quarto/' .gitignore || printf '\n/.quarto/\n/_site/\n' >> .gitignore
echo "done - run: quarto preview"
