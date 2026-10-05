#!/usr/bin/env bash
# Create a new notebook entry (or pipeline page) with labels picked from _labels.yml.
#
#   ./new_post.sh                    # notebook entry
#   ./new_post.sh troubleshooting    # notebook entry using the troubleshooting template
#   ./new_post.sh pipeline           # new page in pipelines/
#
# Then preview with `quarto preview`, and commit + push when happy (the site rebuilds itself).
# Works with the bash 3.2 that ships with macOS.
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-post}"
TODAY=$(date '+%Y-%m-%d')

# print the items of one section of _labels.yml
labels() {
  awk -v sec="$1" '
    /^[a-z_]+:/ { insec = ($1 == sec ":") ; next }
    insec && /^  - / { sub(/^  - /, ""); print }
  ' _labels.yml
}

# pick_many SECTION PROMPT -> prints chosen labels as "a", "b"
pick_many() {
  local i=1 line choice out="" items=()
  while IFS= read -r line; do items+=("$line"); done < <(labels "$1")
  echo "" >&2; echo "$2 (numbers separated by spaces, Enter for none):" >&2
  for line in "${items[@]}"; do printf "  %2d) %s\n" "$i" "$line" >&2; i=$((i+1)); done
  read -r -p "> " choice
  for i in $choice; do
    if [[ "$i" =~ ^[0-9]+$ ]] && (( i >= 1 && i <= ${#items[@]} )); then
      out+="${out:+, }\"${items[$((i-1))]}\""
    fi
  done
  echo "$out"
}

# pick_one SECTION PROMPT -> prints one label (or a new free-text value)
pick_one() {
  local i=1 line choice items=()
  while IFS= read -r line; do items+=("$line"); done < <(labels "$1")
  echo "" >&2; echo "$2:" >&2
  for line in "${items[@]}"; do printf "  %2d) %s\n" "$i" "$line" >&2; i=$((i+1)); done
  echo "   0) none / type a new one" >&2
  read -r -p "> " choice
  if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#items[@]} )); then
    echo "${items[$((choice-1))]}"
  elif [[ "$choice" == "0" ]]; then
    read -r -p "New value (Enter for none): " choice
    echo "$choice"
  else
    echo "$choice"
  fi
}

read -r -p "Title: " TITLE
SLUG=$(echo "$TITLE" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/-/g; s/^-+|-+$//g' | cut -c1-70)

if [[ "$MODE" == "pipeline" ]]; then
  TEMPLATE=_templates/pipeline.qmd
  FILE="pipelines/${SLUG}.qmd"
  TYPE="Pipeline"; PROJECT=""
  CATS=$(pick_many methods "Methods")
else
  if [[ "$MODE" == "troubleshooting" ]]; then
    TEMPLATE=_templates/troubleshooting-post.qmd; TYPE="Troubleshooting"
  else
    TEMPLATE=_templates/notebook-post.qmd
    TYPE=$(pick_one types "Type of entry")
  fi
  PROJECT=$(pick_one projects "Project")
  ORGS=$(pick_many organisms "Organisms")
  METHODS=$(pick_many methods "Methods")
  CATS="\"${TYPE}\"${ORGS:+, $ORGS}${METHODS:+, $METHODS}"
  FILE="notebook/${TODAY}-${SLUG}.qmd"
fi

read -r -p $'\nTools (comma separated, e.g. Cell Ranger, Seurat; Enter for none): ' TOOLS_IN
TOOLS=$(echo "$TOOLS_IN" | awk -F',' '{for(i=1;i<=NF;i++){gsub(/^ +| +$/,"",$i); if($i!="") printf "%s\"%s\"", (n++?", ":""), $i}}')

if [[ -e "$FILE" ]]; then echo "$FILE already exists, not overwriting." >&2; exit 1; fi

# escape characters that are special to sed's replacement text
esc() { printf '%s' "$1" | sed -e 's/[\/&|]/\\&/g'; }
sed -e "s|__TITLE__|$(esc "$TITLE")|" \
    -e "s|__DATE__|${TODAY}|g" \
    -e "s|__TYPE__|$(esc "$TYPE")|" \
    -e "s|__PROJECT__|$(esc "$PROJECT")|" \
    -e "s|__CATEGORIES__|$(esc "$CATS")|" \
    -e "s|__TOOLS__|$(esc "$TOOLS")|" \
    -e "s|__SLUG__|${SLUG}|g" "$TEMPLATE" > "$FILE"

if [[ "$MODE" == "pipeline" ]]; then
  mkdir -p pipelines/scripts
  [[ -e "pipelines/scripts/${SLUG}.sh" ]] || printf '#!/bin/bash\n# %s\n' "$TITLE" > "pipelines/scripts/${SLUG}.sh"
  echo "Also created pipelines/scripts/${SLUG}.sh (the page includes it)."
fi

if [[ -n "$PROJECT" ]] && ! grep -qxF "  - $PROJECT" _labels.yml; then
  echo "Note: '$PROJECT' is a new project. Add it to _labels.yml and give it a section on projects.qmd."
fi

echo "Created $FILE"
"${EDITOR:-nano}" "$FILE"
