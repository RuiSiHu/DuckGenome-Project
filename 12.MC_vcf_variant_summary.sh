#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: bash $0 -v <input.vcf.gz> -o <output_dir>"
    exit 1
}

VCF=""
OUTDIR=""

while getopts "v:o:h" opt; do
    case $opt in
        v) VCF="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -z "$VCF" || -z "$OUTDIR" ]] && usage
[[ ! -f "$VCF" ]] && echo "[ERROR] VCF not found: $VCF" && exit 1

mkdir -p "$OUTDIR"

# 1. basic info
bcftools query -l "$VCF" > "$OUTDIR/sample_names.txt"
bcftools view -h "$VCF" | grep '^##INFO' > "$OUTDIR/vcf_info_header.txt" || true

# 2. standard bcftools stats
bcftools stats "$VCF" > "$OUTDIR/vcf.stats.txt"
grep '^SN' "$OUTDIR/vcf.stats.txt" > "$OUTDIR/summary_SN.txt"
grep '^TSTV' "$OUTDIR/vcf.stats.txt" > "$OUTDIR/summary_TSTV.txt"

{
    echo -e "Category\tCount"
    echo -e "SNP\t$(bcftools view -v snps -H "$VCF" | wc -l)"
    echo -e "MNP\t$(bcftools view -v mnps -H "$VCF" | wc -l)"
    echo -e "INDEL\t$(bcftools view -v indels -H "$VCF" | wc -l)"
    echo -e "OTHER\t$(bcftools view -v other -H "$VCF" | wc -l)"
} > "$OUTDIR/bcftools_variant_counts.tsv"

bcftools view -v other -H "$VCF" | head -20 > "$OUTDIR/other_variants_head20.txt" || true

# 3. length-based site classification
bcftools query -f '%REF\t%ALT\n' "$VCF" | awk '
BEGIN{FS="\t"}
{
    ref=$1; split($2,alts,","); maxdiff=0; has_symbolic=0
    for(i in alts){
        alt=alts[i]
        if(alt ~ /^<.*>$/){has_symbolic=1; continue}
        diff=length(alt)-length(ref)
        if(diff<0) diff=-diff
        if(diff>maxdiff) maxdiff=diff
    }
    if(has_symbolic && maxdiff==0) symbolic++
    else if(maxdiff>=50) large_var++
    else if(maxdiff>=1) small_indel++
    else if(length(ref)==1) snp_like++
    else mnp_like++
}
END{
    print "Category\tCount"
    print "SNP_like\t"snp_like
    print "MNP_like\t"mnp_like
    print "Small_indels_1_49bp\t"small_indel
    print "Large_variants_ge50bp\t"large_var
    print "Symbolic_or_unresolved\t"symbolic
}' > "$OUTDIR/length_based_site_classification.tsv"

# 4. length-based allele classification
bcftools query -f '%REF\t%ALT\n' "$VCF" | awk '
BEGIN{FS="\t"}
{
    ref=$1; split($2,alts,",")
    for(i in alts){
        alt=alts[i]
        if(alt ~ /^<.*>$/) continue
        diff=length(alt)-length(ref)
        if(diff>=50) large_ins++
        else if(diff<=-50) large_del++
        else if(diff>0) small_ins++
        else if(diff<0) small_del++
        else same_len++
    }
}
END{
    print "Category\tCount"
    print "large_insertions_ge50bp\t"large_ins
    print "large_deletions_ge50bp\t"large_del
    print "small_insertions_1_49bp\t"small_ins
    print "small_deletions_1_49bp\t"small_del
    print "same_length_substitutions\t"same_len
}' > "$OUTDIR/length_based_allele_classification.tsv"

# 5. combined report
{
    echo "### sample_names"
    cat "$OUTDIR/sample_names.txt"
    echo

    echo "### bcftools_summary_SN"
    cat "$OUTDIR/summary_SN.txt"
    echo

    echo "### bcftools_summary_TSTV"
    cat "$OUTDIR/summary_TSTV.txt"
    echo

    echo "### bcftools_variant_counts"
    cat "$OUTDIR/bcftools_variant_counts.tsv"
    echo

    echo "### length_based_site_classification"
    cat "$OUTDIR/length_based_site_classification.tsv"
    echo

    echo "### length_based_allele_classification"
    cat "$OUTDIR/length_based_allele_classification.tsv"
    echo
} > "$OUTDIR/variant_summary_report.txt"

echo "[INFO] Done: $OUTDIR"