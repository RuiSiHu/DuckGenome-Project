#!/usr/bin/env python3
import sys
import gzip
from collections import Counter

if len(sys.argv) != 3:
    sys.exit("Usage: python stat_panduck_vcf.py <input.vcf.gz> <output.tsv>")

vcf_file = sys.argv[1]
out_file = sys.argv[2]

SV_THRESHOLD = 50  # length difference >= 50 bp as SV
type_categories = ["SNP", "MNP", "InDels", "SV", "BND", "Other"]
direction_categories = ["Insertion", "Deletion"]
site_direction_categories = ["Insertion", "Deletion", "Mixed_or_NA"]
multiallelic_site_categories = [
    "All-SNP multiallelic sites",
    "All-MNP multiallelic sites",
    "All-InDel multiallelic sites",
    "All-SV multiallelic sites",
    "All-BND multiallelic sites",
    "Mixed-type multiallelic sites",
    "Other multiallelic sites",
]
at_branch_categories = ["no_AT", "single_path", "bifurcation", "high_branching"]


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_info(info_str):
    d = {}
    if not info_str or info_str == ".":
        return d
    for x in info_str.split(";"):
        if not x:
            continue
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
        else:
            d[x] = True
    return d


def is_symbolic(alt):
    return alt.startswith("<") and alt.endswith(">")


def is_bnd(alt):
    return ("[" in alt) or ("]" in alt)


def symbolic_label(alt):
    if is_symbolic(alt):
        return alt[1:-1].split(":")[0].upper()
    return ""


def get_info_value_by_alt(info_dict, key, idx):
    val = str(info_dict.get(key, "")).strip()
    if not val or val == ".":
        return ""
    parts = [x.strip() for x in val.split(",")]
    if idx < len(parts) and parts[idx] and parts[idx] != ".":
        return parts[idx]
    return parts[0] if parts and parts[0] != "." else ""


def classify_alt(ref, alt, svtype):
    if alt in {".", "*"}:
        return "Other"
    if is_bnd(alt):
        return "BND"
 
    svtype = str(svtype).upper()
    sym_label = symbolic_label(alt)
    if is_symbolic(alt):
        if svtype in {"BND", "TRA"} or sym_label in {"BND", "TRA"}:
            return "BND"
        return "SV"

    rlen = len(ref)
    alen = len(alt)

    if rlen == 1 and alen == 1:
        return "SNP"
    if rlen == alen and rlen > 1:
        return "MNP"
    if rlen != alen:
        if abs(alen - rlen) >= SV_THRESHOLD:
            return "SV"
        return "InDels"
    # Equal length complex substitutions are MNP; otherwise Other.
    return "MNP" if rlen > 1 else "Other"


def parse_svlen_for_alt(info_dict, idx):
    raw = get_info_value_by_alt(info_dict, "SVLEN", idx)
    if not raw:
        return None
    try:
        return int(float(raw))
    except ValueError:
        return None


def indel_direction(ref, alt, svtype="", svlen=None):
    if is_bnd(alt):
        return "NA"
    if is_symbolic(alt):
        svtype = str(svtype).upper()
        sym_label = symbolic_label(alt)

        if sym_label in {"INS", "DUP"}:
            return "Insertion"
        if sym_label == "DEL":
            return "Deletion"

        if svtype == "INS":
            return "Insertion"
        if svtype == "DEL":
            return "Deletion"
        if svtype == "DUP":
            return "Insertion"
        if isinstance(svlen, int):
            if svlen > 0:
                return "Insertion"
            if svlen < 0:
                return "Deletion"
        return "NA"
    rlen = len(ref)
    alen = len(alt)
    if alen > rlen:
        return "Insertion"
    elif alen < rlen:
        return "Deletion"
    else:
        return "NA"


def is_transition(ref, alt):
    pair = (ref.upper(), alt.upper())
    return pair in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def at_branch_type(info_dict):
    at = str(info_dict.get("AT", "")).strip()
    if not at or at == ".":
        return "no_AT"
    n = len([x.strip() for x in at.split(",") if x.strip() and x.strip() != "."])
    if n == 0:
        return "no_AT"
    if n <= 1:
        return "single_path"
    elif n == 2:
        return "bifurcation"
    else:
        return "high_branching"


def get_multiallelic_site_label(class_set):
    if len(class_set) != 1:
        return "Mixed-type multiallelic sites"

    only = next(iter(class_set))
    mapping = {
        "SNP": "All-SNP multiallelic sites",
        "MNP": "All-MNP multiallelic sites",
        "InDels": "All-InDel multiallelic sites",
        "SV": "All-SV multiallelic sites",
        "BND": "All-BND multiallelic sites",
    }
    return mapping.get(only, "Other multiallelic sites")


def update_direction_subcounts(
    variant_class,
    ref,
    alt,
    info_dict,
    idx,
    sv_sub_counter,
    indel_sub_counter,
):
    if variant_class == "SV":
        svtype = get_info_value_by_alt(info_dict, "SVTYPE", idx)
        svlen = parse_svlen_for_alt(info_dict, idx)
        direction = indel_direction(ref, alt, svtype=svtype, svlen=svlen)
        if direction in {"Insertion", "Deletion"}:
            sv_sub_counter[direction] += 1
    elif variant_class == "InDels":
        direction = indel_direction(ref, alt)
        if direction in {"Insertion", "Deletion"}:
            indel_sub_counter[direction] += 1


def get_site_direction(directions):
    valid = [d for d in directions if d in {"Insertion", "Deletion"}]
    if len(valid) != len(directions) or not valid:
        return "Mixed_or_NA"
    return valid[0] if len(set(valid)) == 1 else "Mixed_or_NA"


site_cnt = Counter()
biallelic_cnt = Counter()
multiallelic_site_cnt = Counter()
multiallelic_alt_cnt = Counter()

biallelic_sv_sub = Counter()
biallelic_indel_sub = Counter()
multiallelic_alt_sv_sub = Counter()
multiallelic_alt_indel_sub = Counter()
multiallelic_site_sv_sub = Counter()
multiallelic_site_indel_sub = Counter()

alt_count_dist = Counter()     
at_branch_cnt = Counter()       # no_AT / single / bi / high

sample_count = 0
ts = 0
tv = 0

with open_vcf(vcf_file) as f:
    for line in f:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            sample_count = max(0, len(cols) - 9)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        ref = fields[3]
        alt_field = fields[4]
        info_dict = parse_info(fields[7])

        if alt_field == ".":
            continue

        alts = alt_field.split(",")
        site_cnt["Total variant sites"] += 1
        alt_count_dist[str(len(alts))] += 1
        at_branch_cnt[at_branch_type(info_dict)] += 1

        alt_classes = []
        for idx, alt in enumerate(alts):
            svtype = get_info_value_by_alt(info_dict, "SVTYPE", idx)
            c = classify_alt(ref, alt, svtype)
            alt_classes.append(c)

            if c == "SNP":
                if is_transition(ref, alt):
                    ts += 1
                else:
                    tv += 1

        # biallelic
        if len(alts) == 1:
            site_cnt["Biallelic sites"] += 1
            c = alt_classes[0]
            biallelic_cnt[c] += 1

            update_direction_subcounts(
                c, ref, alts[0], info_dict, 0, biallelic_sv_sub, biallelic_indel_sub
            )

        # multiallelic
        else:
            site_cnt["Multiallelic sites"] += 1

            for idx, (alt, c) in enumerate(zip(alts, alt_classes)):
                multiallelic_alt_cnt[c] += 1

                update_direction_subcounts(
                    c,
                    ref,
                    alt,
                    info_dict,
                    idx,
                    multiallelic_alt_sv_sub,
                    multiallelic_alt_indel_sub,
                )

            class_set = set(alt_classes)
            multiallelic_site_cnt[get_multiallelic_site_label(class_set)] += 1
            if len(class_set) == 1:
                only_class = next(iter(class_set))
                if only_class in {"InDels", "SV"}:
                    directions = []
                    for idx, alt in enumerate(alts):
                        if only_class == "SV":
                            svtype = get_info_value_by_alt(info_dict, "SVTYPE", idx)
                            svlen = parse_svlen_for_alt(info_dict, idx)
                            direction = indel_direction(ref, alt, svtype=svtype, svlen=svlen)
                        else:
                            direction = indel_direction(ref, alt)
                        directions.append(direction)

                    site_direction = get_site_direction(directions)
                    if only_class == "SV":
                        multiallelic_site_sv_sub[site_direction] += 1
                    else:
                        multiallelic_site_indel_sub[site_direction] += 1


total_sites = site_cnt["Total variant sites"]
biallelic_sites = site_cnt["Biallelic sites"]
multiallelic_sites = site_cnt["Multiallelic sites"]

tstv_ratio = "NA"
if tv > 0:
    tstv_ratio = f"{ts / tv:.4f}"


with open(out_file, "w") as out:
    out.write("Section\tCategory\tCount\n")

    # Overall
    overall_rows = [
        ("Overall", "Total variant sites", total_sites),
        ("Overall", "Biallelic sites", biallelic_sites),
        ("Overall", "Multiallelic sites", multiallelic_sites),
    ]
    for sec, cat, cnt in overall_rows:
        out.write(f"{sec}\t{cat}\t{cnt}\n")

    # Biallelic (final table format)
    biallelic_rows = [
        ("Biallelic", "SNP", biallelic_cnt["SNP"]),
        ("Biallelic", "MNP", biallelic_cnt["MNP"]),
        ("Biallelic", "Small insertion", biallelic_indel_sub["Insertion"]),
        ("Biallelic", "Small deletion", biallelic_indel_sub["Deletion"]),
        ("Biallelic", "Insertion (SV)", biallelic_sv_sub["Insertion"]),
        ("Biallelic", "Deletion (SV)", biallelic_sv_sub["Deletion"]),
    ]
    for sec, cat, cnt in biallelic_rows:
        out.write(f"{sec}\t{cat}\t{cnt}\n")

    # Multiallelic_sites (site-level final table format)
    multiallelic_site_rows = [
        ("Multiallelic_sites", "SNP", multiallelic_site_cnt["All-SNP multiallelic sites"]),
        ("Multiallelic_sites", "MNP", multiallelic_site_cnt["All-MNP multiallelic sites"]),
        ("Multiallelic_sites", "Small insertion", multiallelic_site_indel_sub["Insertion"]),
        ("Multiallelic_sites", "Small deletion", multiallelic_site_indel_sub["Deletion"]),
        ("Multiallelic_sites", "Insertion (SV)", multiallelic_site_sv_sub["Insertion"]),
        ("Multiallelic_sites", "Deletion (SV)", multiallelic_site_sv_sub["Deletion"]),
        (
            "Multiallelic_sites",
            "Mixed-type",
            multiallelic_site_cnt["Mixed-type multiallelic sites"]
            + multiallelic_site_indel_sub["Mixed_or_NA"]
            + multiallelic_site_sv_sub["Mixed_or_NA"],
        ),
    ]
    for sec, cat, cnt in multiallelic_site_rows:
        out.write(f"{sec}\t{cat}\t{cnt}\n")

    # Multiallelic ALT (ALT-level final table format)
    multiallelic_alt_rows = [
        ("Multiallelic ALT", "SNP", multiallelic_alt_cnt["SNP"]),
        ("Multiallelic ALT", "MNP", multiallelic_alt_cnt["MNP"]),
        ("Multiallelic ALT", "Small insertion", multiallelic_alt_indel_sub["Insertion"]),
        ("Multiallelic ALT", "Small deletion", multiallelic_alt_indel_sub["Deletion"]),
        ("Multiallelic ALT", "Insertion (SV)", multiallelic_alt_sv_sub["Insertion"]),
        ("Multiallelic ALT", "Deletion (SV)", multiallelic_alt_sv_sub["Deletion"]),
    ]
    for sec, cat, cnt in multiallelic_alt_rows:
        out.write(f"{sec}\t{cat}\t{cnt}\n")

    for k in sorted(alt_count_dist.keys(), key=lambda x: int(x)):
        cnt = alt_count_dist[k]
        out.write(f"Graph_ALT_count\tALT={k}\t{cnt}\n")

    out.write(f"Graph_AT\tbifurcation\t{at_branch_cnt['bifurcation']}\n")
    out.write(f"Graph_AT\thigh_branching\t{at_branch_cnt['high_branching']}\n")

    # Summary
    out.write(f"Summary\tTs/Tv ratio\t{tstv_ratio}\n")
    out.write(f"Summary\tNumber of genomes\t{sample_count}\n")
    out.write(f"Summary\tSV threshold (>bp)\t{SV_THRESHOLD}\n")

print(f"Output written to: {out_file}")
