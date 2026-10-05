#!/usr/bin/env python3
"""
Record which old notebook posts used which wet-lab protocol page, by filling the
`protocols` column of _migration/labels.csv. Then re-run:

    python3 _migration/protocol_links.py
    python3 _migration/migrate_posts.py convert

For NEW posts you don't need this: just add  protocols: ["<protocol-file-name>"]
to the post's front matter.
"""
import csv, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "_migration", "labels.csv")

# protocol page (file name in protocols/, without .qmd) -> old post filename prefixes
PROTOCOL_POSTS = {
    "dna-rna-coextraction-zymo": [
        "2019-02-13-Zymo-DNA-RNA-Extraction-Protocol",
        "2019-03-13-Zymo-DNA-RNA-Extract-P.astreoides-Genome",
        "2019-12-04-DNA-RNA-Extractions-Thermal-Transplant",
        "2019-12-05-DNA-RNA-Extractions-Thermal-Transplant",
        "2020-01-13-DNA-RNA-Extractions-Thermal-Transplant",
        "2020-07-30-DNA-RNA-Extractions-Thermal-Transplant",
        "2020-07-31-DNA-RNA-Extraxtion-Thermal-Transplant",
        "2020-08-05-DNA-RNA-Extractions-Porites-astreoides",
        "2020-08-06-DNA-RNA-extractions-Porites-Thermal-Transplant-Larvae",
        "2020-08-12-DNA-RNA-Extractions-on-P.-astreoides-larvae",
        "2020-08-13-DNA-RNA-Extractions-on-P.-astreoides-larvae",
        "2020-08-16-DNA-RNA-Extractions-on-P.-astreoides-larvae",
        "2020-08-19-DNA-RNA-Extractions-on-P.-astreoides-larvae",
        "2020-08-20-DNA-RNA-Extractions-on-P.-astreoides-larvae",
        "2020-10-23-20201022-DNA-RNA-Extractions",
        "2020-10-26-20201025-DNA-RNA-Extractions",
        "2020-10-27-20201027-DNA-RNA-Extractions",
        "2020-10-29-20201029-DNA-RNA-Extractions",
        "2020-11-09-20201106-DNA-RNA-Extractions",
        "2020-11-12-20201110-DNA-RNA-Extractions",
        "2020-11-12-20201111-DNA-RNA-Extractions",
        "2020-11-18-20201117-DNA-RNA-Extractions",
        "2020-11-23-20201119-DNA-RNA-Extractions",
        "2020-11-23-20201120-DNA-RNA-Extractions",
        "2020-11-25-20201125-DNA-RNA-Extractions",
        "2020-11-30-20201126-DNA-RNA-Extractions",
        "2020-12-03-20201202-DNA-RNA-Extractions",
    ],
    "coral-tissue-removal-homogenization": [
        "2019-05-14-Adult-coral-Homogenate",
        "2021-06-11-Airbrushing-Protocol",
        "2019-10-04-Airbrushing-and-Homogenizing",
        "2022-11-29-Galaxea-Airbrush",
        "2019-05-15-Citrate-Synthase-for",
        "2019-11-26-Total-Antioxidant",
        "2021-06-23-20210622-Lipid",
        "2021-06-25-20210624-Lipid",
        "2021-06-30-20210629-Total-Protein",
        "2021-07-14-20210713-Carbohydrate",
    ],
    "symbiont-density-hemocytometer": [
        "2018-08-14-Symbiodinium-Density",
        "2022-11-29-Galaxea-Airbrush",
        "2023-01-23-Coral-processing-for-TEM",
    ],
    "total-protein-bca": [
        "2018-10-05-Total-Protein-Extraction",
        "2021-06-30-20210629-Total-Protein",
        "2019-03-14-Citrate-Synthase-Troubleshooting",
        "2019-05-15-Citrate-Synthase-for",
        "2019-11-26-Total-Antioxidant",
    ],
    "citrate-synthase-activity": [
        "2019-03-14-Citrate-Synthase-Troubleshooting",
        "2019-05-15-Citrate-Synthase-for",
    ],
    "sam-elisa": ["2019-05-15-SAM-ELISA"],
    "total-antioxidant-capacity": ["2019-11-26-Total-Antioxidant"],
    "metabolite-extraction": [
        "2020-02-03-Metabolomics",
        "2020-10-07-20201006", "2020-10-09-20201008", "2020-10-12-20201009",
        "2020-10-14-20201013", "2020-10-15-20201015", "2020-10-20-20201020",
        "2021-03-22-20210321-Symbiont", "2021-03-25-20210325-Symbiont",
        "2021-04-16-Symbiont-Integration",
    ],
    "wgbs-library-prep-picomethyl": [
        "2021-04-16-Thermal-Transplant-WGBS",
        "2021-04-01-20210401-WGBS", "2021-04-15-20210415-WGBS", "2021-04-23-20210422-WGBS",
        "2021-04-27-20210426-WGBS", "2021-04-29-20210428-WGBS", "2021-05-04-20210504-WGBS",
        "2021-05-06-20210506-WGBS", "2021-05-10-20210510-WGBS", "2021-05-12-20210512-WGBS",
        "2021-05-26-20210526-WGBS",
    ],
    "its2-amplicon-pcr": [
        "2022-03-16-Touchdown-PCR",
        "2021-02-18-Thermal-Transplant-ITS2", "2021-03-04-20210302-Thermal-Transplant-ITS2",
        "2021-03-10-20210309-Thermal-Transplant-ITS2", "2021-03-25-20210325-Thermal-Transplant-ITS2",
        "2021-11-04-20211104-ITS2", "2021-11-09-20211109-ITS2", "2022-03-21-20220319-PJB-ITS2",
    ],
    "coral-fixation-tem": ["2023-01-23-Coral-processing-for-TEM"],
    "nanopore-direct-rna-minion": ["2023-10-24-Direct-RNA"],
}

# Old "Protocol" posts that have been rewritten as a protocol page: they stay in the
# notebook as the record of what was done, so they are no longer typed "Protocol".
REWRITTEN = {
    "2019-02-13-Zymo-DNA-RNA-Extraction-Protocol": "Lab work",
    "2020-08-12-DNA-RNA-Extractions-on-P.-astreoides-larvae": "Lab work",
    "2018-08-14-Symbiodinium-Density": "Lab work",
    "2018-10-05-Total-Protein-Extraction": "Lab work",
    "2019-03-14-Citrate-Synthase-Troubleshooting": "Troubleshooting",
    "2019-05-14-Adult-coral-Homogenate": "Lab work",
    "2019-05-15-SAM-ELISA": "Lab work",
    "2019-11-26-Total-Antioxidant": "Lab work",
    "2020-02-03-Metabolomics": "Lab work",
    "2021-04-16-Thermal-Transplant-WGBS": "Lab work",
    "2021-06-11-Airbrushing-Protocol": "Lab work",
    "2022-03-16-Touchdown-PCR": "Lab work",
    "2023-01-23-Coral-processing-for-TEM": "Lab work",
    "2023-10-24-Direct-RNA": "Lab work",
}

rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
fields = list(rows[0])
if "protocols" not in fields:
    fields.append("protocols")
for r in rows:
    used = [slug for slug, prefixes in PROTOCOL_POSTS.items()
            if any(r["file"].startswith(p) for p in prefixes)]
    r["protocols"] = "; ".join(used)
    for prefix, new_type in REWRITTEN.items():
        if r["file"].startswith(prefix):
            r["type"] = new_type
with open(CSV, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
for slug in PROTOCOL_POSTS:
    n = sum(slug in r["protocols"] for r in rows)
    print(f"{slug}: {n} notebook posts")
