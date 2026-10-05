#!/usr/bin/env python3
"""Warn about notebook labels that aren't in _labels.yml (typos, near-duplicates).
Usage: python3 check_labels.py"""
import glob, re, sys
import yaml  # pip install pyyaml

vocab = yaml.safe_load(open("_labels.yml"))
allowed_cats = set(vocab["types"]) | set(vocab["organisms"]) | set(vocab["methods"])
problems = 0
for f in sorted(glob.glob("notebook/*.md") + glob.glob("notebook/*.qmd")):
    m = re.match(r"^---\s*\n(.*?)\n---", open(f, encoding="utf-8").read(), re.S)
    if not m:
        print(f"{f}: no front matter"); problems += 1; continue
    fm = yaml.safe_load(m.group(1)) or {}
    for c in fm.get("categories") or []:
        if c not in allowed_cats:
            print(f"{f}: unknown category '{c}'"); problems += 1
    if fm.get("type") not in vocab["types"]:
        print(f"{f}: unknown type '{fm.get('type')}'"); problems += 1
    if fm.get("project") and fm["project"] not in vocab["projects"]:
        print(f"{f}: unknown project '{fm['project']}'"); problems += 1
print("all labels OK" if not problems else f"{problems} problem(s)")
sys.exit(1 if problems else 0)
