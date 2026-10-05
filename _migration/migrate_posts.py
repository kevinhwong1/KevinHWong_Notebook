#!/usr/bin/env python3
"""
One-time migration of the old Jekyll notebook (_posts/) to Quarto (notebook/).

Two steps, so the labels can be reviewed before anything is converted:

  python3 _migration/migrate_posts.py labels    # writes _migration/labels.csv (rule-based guesses)
  # ...open labels.csv in Excel/R, fix any type/project/categories/tools...
  python3 _migration/migrate_posts.py convert   # writes notebook/*.md using labels.csv

Re-running `convert` overwrites notebook/ files generated from _posts, so edit
labels.csv (not the converted posts) until you are happy, then delete _posts.

What `convert` does for every post:
  * new front matter: title, date, type, project, categories, tools, aliases
  * alias = the old Jekyll URL (/Old-Title/) so every existing link redirects
  * rewrites {{ site.baseurl }}/images/... and github.com/.../blob/master/images/... to /images/...
  * rewrites links to other posts (old site URLs and github _posts links) to the new pages
"""
import csv, glob, os, re, sys, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "_posts")
OUT = os.path.join(ROOT, "notebook")
CSV = os.path.join(ROOT, "_migration", "labels.csv")

FNAME_RE = re.compile(r"^(\d{4})-(\d{2})-(\d{1,2})[-_](.*)\.md$")

# ---------------------------------------------------------------------------
# Controlled vocabulary (keep in sync with _labels.yml)
# ---------------------------------------------------------------------------
TYPES = ["Lab work", "Protocol", "Analysis", "Troubleshooting", "Central doc", "Planning"]

ORGANISM_RULES = [
    ("Porites astreoides", r"porites|astreoides|p\.?\s?astreoides|\bpast\b|pjb|thermal.transplant|mansour"),
    ("Galaxea fascicularis", r"galaxea|gfas"),
    ("Mnemiopsis leidyi", r"mnemi"),
    ("Nematostella vectensis", r"nvec|nematostella"),
    ("Acropora cervicornis", r"\bacer\b|cervicornis|ehrens"),
    ("Pocillopora damicornis", r"\bpdam\b|damicornis|ehrens"),
    ("Montipora capitata", r"montipora|mcap"),
    ("Pocillopora acuta", r"\bacuta\b"),
    ("Astrangia poculata", r"astrangia|poculata"),
    ("Geoduck", r"geoduck"),
    ("Symbiodiniaceae", r"symbiodin|dtrenchii|symbiont.mapping|\bits2\b|symbiont.density"),
]

METHOD_RULES = [
    ("DNA/RNA extraction", r"extract(ion|s)?\b.*(dna|rna)|(dna|rna).*extract|zymo|extraxtion|nanodrop"),
    ("Physiology assays", r"protein|lipid|carbohydrate|citrate|sam.elisa|\bsam\b|antioxidant|\btac\b|symbiodinium.density"),
    ("Tissue & cell prep", r"airbrush|homogen|nuclei|cell.count|welcome.party|separation"),
    ("Metabolomics", r"metabolom"),
    ("ITS2", r"\bits2\b"),
    ("DNA methylation", r"wgbs|wbgs|methyl|bissnp|bs.snper|epidiverse|picomethyl"),
    ("Transcriptomics", r"tag.?seq|transcriptom|hisat"),
    ("scRNA-seq", r"scrna|single.cell|cellranger|cell.ranger|velocity|velocyto|geneext|symbiont.mapping|hcr"),
    ("Genome annotation", r"genome|maker|augustus|busco|annotation"),
    ("Functional annotation", r"kofam"),
    ("Comparative genomics", r"orthofinder"),
    ("Imaging", r"\btem\b|\bhcrs?\b|microscop"),
    ("Nanopore", r"minion|nanopore"),
    ("16S", r"\b16s\b"),
    ("HPC & computing", r"pegasus|conda|\bhpc\b|lab.notebook|computing"),
]

TOOL_RULES = [
    ("Cell Ranger", r"cellranger|cell.ranger"), ("GeneExt", r"geneext"), ("velocyto", r"velocyto"),
    ("OrthoFinder", r"orthofinder"), ("BUSCO", r"busco"), ("MAKER", r"maker"), ("AUGUSTUS", r"augustus"),
    ("KofamScan", r"kofam"), ("nf-core/methylseq", r"nf.core|methylseq"), ("BS-SNPer", r"bs.snper"),
    ("EpiDiverse", r"epidiverse"), ("BisSNP", r"bissnp"), ("HISAT2", r"hisat"), ("conda", r"conda"),
]

PROJECT_RULES = [  # first match wins
    ("Porites Thermal Transplant", r"thermal.transplant|thermal_transp|methylseq|wgbs|wbgs|bissnp|bs.snper|epidiverse|m.bias"),
    ("Porites July Bleaching", r"july.bleaching|pjb|porites.bleaching|kw.ah.es|tag.seq.samples|mansour"),
    ("Porites Nutrition", r"porites.nutrition|june.\(?patch"),
    ("Astrangia Nutrition", r"astrangia|poculata.nutrition"),
    ("Symbiont Integration", r"symbiont.integration|ah.tag.seq|mcap2020"),
    ("Porites astreoides Genome", r"p\.?\s?astreoides.genome|past.genome|maker|augustus|busco.on.p|busco.*astreoides|kofam.*astreoides|porites.astreoides.genome"),
    ("Cnidarian Stem Cells", r"ehrens|\bpdam\b|\bacer\b|nvec|nematostella"),
    ("Mnemiopsis Phagocytes", r"mnemi"),
    ("Dark Genes", r"gfas|galaxea.thermal|galaxea.airbrush|galaxea.fasicularis|khw.nar|orthofinder|dark.genes|dtrenchii"),
]


# Hand corrections where the rules above guess wrong (filename prefix -> fields to replace)
OVERRIDES = {
    "2019-03-13-Zymo-DNA-RNA-Extract": dict(methods="DNA/RNA extraction"),
    "2020-10-29-20201029-DNA-RNA": dict(project="Porites July Bleaching"),
    "2020-11-09-20201106-DNA-RNA": dict(project="Porites July Bleaching"),
    "2021-03-22-20210321-Symbiont": dict(organisms="Montipora capitata"),
    "2021-03-25-20210325-Symbiont": dict(organisms="Montipora capitata"),
    "2021-08-26-20210826-Lipid": dict(project=""),
    "2022-11-21-Testing-BS-SNPer": dict(organisms="Porites astreoides"),
    "2023-03-15-EpiDiverse": dict(organisms="Porites astreoides"),
    "2023-09-05-Gfas-Dark-Genes": dict(tools="Cell Ranger"),
    "2024-03-27-Testing-CellRanger-mito": dict(methods="scRNA-seq"),
    "2024-07-24-CellRanger-on-KHW-NAR": dict(organisms="Galaxea fascicularis"),
    "2024-07-24-scRNAseq-Velocity": dict(tools="velocyto; scVelo"),
    "2024-10-30-Pegasus-conda": dict(type="Protocol"),
    "2024-11-27-Orthofinder-between": dict(methods="scRNA-seq; Comparative genomics"),
    "2025-03-04_Troubleshooting_velocyto": dict(methods="scRNA-seq", tools="velocyto"),
}


def jekyll_slug(s):
    """Reproduce Jekyll's :title slug (mode: pretty, cased) for old-URL aliases."""
    s = re.sub(r"[^A-Za-z0-9._~!$&'()+,;=@]+", "-", s)
    return s.strip("-")


def new_slug(s):
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", "-", s)
    return s.strip("-")[:80].strip("-")


def read_post(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.match(r"^---\s*\n(.*?)\n---\s*\n?", txt, re.S)
    fm, body = (m.group(1), txt[m.end():]) if m else ("", txt)
    meta = {}
    for line in fm.splitlines():
        mm = re.match(r"^(\w+):\s*(.*)$", line)
        if mm:
            meta[mm.group(1)] = mm.group(2).strip().strip("'\"")
    return meta, body


def posts():
    for path in sorted(glob.glob(os.path.join(SRC, "*.md"))):
        fn = os.path.basename(path)
        m = FNAME_RE.match(fn)
        if not m:
            print("  skipping (not a dated post):", fn)
            continue
        y, mo, d, rest = m.groups()
        date = f"{y}-{mo}-{int(d):02d}"
        meta, body = read_post(path)
        title = meta.get("title") or rest.replace("-", " ")
        # Jekyll only published files named YYYY-MM-DD-title.md; others never had a URL
        old = jekyll_slug(rest) if fn[10] == "-" and len(fn.split("-")[2]) == 2 else ""
        yield dict(file=fn, path=path, date=date, title=title, meta=meta, body=body,
                   old_slug=old, new_name=f"{date}-{new_slug(rest)}")


def find(rules, text, many=True):
    hits = [lab for lab, rx in rules if re.search(rx, text, re.I)]
    return hits if many else (hits[0] if hits else "")


def guess(p):
    meta = p["meta"]
    text = " ".join([p["file"], p["title"], meta.get("tags", ""), meta.get("categories", "")])
    cats = meta.get("categories", "").lower()
    t = text.lower()
    if "central" in t:
        typ = "Central doc"
    elif "goals" in t:
        typ = "Planning"
    elif "troubleshoot" in t:
        typ = "Troubleshooting"
    elif "analysis" in cats:
        typ = "Analysis"
    elif "process" in cats or "proessing" in cats:
        typ = "Lab work"
    else:
        typ = "Protocol"
    project = "" if typ in ("Central doc", "Planning") else find(PROJECT_RULES, text, many=False)
    if "welcome-party" in t:
        project = ""
    g = dict(type=typ, project=project,
             organisms="; ".join(find(ORGANISM_RULES, text)),
             methods="; ".join(find(METHOD_RULES, text)),
             tools="; ".join(find(TOOL_RULES, text)))
    for prefix, fix in OVERRIDES.items():
        if p["file"].startswith(prefix):
            g.update(fix)
    return g


def cmd_labels():
    rows = []
    for p in posts():
        g = guess(p)
        rows.append(dict(file=p["file"], title=p["title"], **g,
                         old_tags=p["meta"].get("tags", ""), old_categories=p["meta"].get("categories", "")))
    with open(CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {CSV}")


def yaml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def split(s):
    return [x.strip() for x in s.split(";") if x.strip()]


def cmd_convert():
    labels = {r["file"]: r for r in csv.DictReader(open(CSV, encoding="utf-8"))}
    plist = list(posts())
    by_old = {p["old_slug"]: p["new_name"] for p in plist if p["old_slug"]}
    by_file = {p["file"]: p["new_name"] for p in plist}
    os.makedirs(OUT, exist_ok=True)

    def to_post(target):
        return f"/notebook/{target}.md"

    def fix_links(body):
        # {{ site.baseurl }}/images -> /images
        body = re.sub(r"\{\{\s*site\.baseurl\s*\}\}/", "/", body)
        # github blob/raw links into this repo's images/ and protocols/ -> local
        body = re.sub(r"https?://(?:github\.com/kevinhwong1/KevinHWong_Notebook/(?:blob|raw)/master|"
                      r"raw\.githubusercontent\.com/kevinhwong1/KevinHWong_Notebook/master)/"
                      r"((?:images|protocols)/[^)\s\"'?]+)(?:\?raw=true)?", r"/\1", body)

        # github links to _posts/<file>.md -> new page
        def gh_post(m):
            fn = urllib.parse.unquote(m.group(1))
            return to_post(by_file[fn]) if fn in by_file else m.group(0)
        body = re.sub(r"https?://github\.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/([^)\s\"'#]+)", gh_post, body)

        # old site URLs -> new page
        def site_post(m):
            slug = urllib.parse.unquote(m.group(1))
            return to_post(by_old[slug]) + (m.group(2) or "") if slug in by_old else m.group(0)
        body = re.sub(r"https?://kevinhwong1\.github\.io/KevinHWong_Notebook/([^/)\s\"'#]+)/?(#[^)\s\"']*)?", site_post, body)
        return body

    n = 0
    for p in plist:
        lab = labels.get(p["file"])
        if lab is None:
            print("  no labels.csv row for", p["file"], "- run `labels` again")
            continue
        cats = [lab["type"]] + split(lab["organisms"]) + split(lab["methods"])
        fm = ["---", f"title: {yaml_str(p['title'])}", f"date: {p['date']}", f"type: {yaml_str(lab['type'])}"]
        fm.append(f"project: {yaml_str(lab['project'])}")  # keep even if empty (listing shows blank)
        fm.append("categories: [" + ", ".join(yaml_str(c) for c in cats) + "]")
        if split(lab.get("protocols") or ""):
            fm.append("protocols: [" + ", ".join(yaml_str(c) for c in split(lab["protocols"])) + "]")
        if split(lab["tools"]):
            fm.append("tools: [" + ", ".join(yaml_str(c) for c in split(lab["tools"])) + "]")
        if p["old_slug"]:
            fm.append(f"aliases:\n  - {yaml_str('/' + p['old_slug'] + '/')}")
        fm.append("---\n")
        with open(os.path.join(OUT, p["new_name"] + ".md"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(fm) + fix_links(p["body"]).lstrip("\n"))
        n += 1
    print(f"converted {n} posts into {OUT}")


if __name__ == "__main__":
    {"labels": cmd_labels, "convert": cmd_convert}.get(
        sys.argv[1] if len(sys.argv) > 1 else "", lambda: print(__doc__))()
