#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  04_prepare_fuma_pair_uploads.sh <separated_clump_pair_dir> <output_root> [pair_name]

Prepare positive and negative lead-SNP and candidate-SNP files for FUMA from
the directionality-separated PleioFDR clump outputs.

Arguments:
  separated_clump_pair_dir  Directory containing *loci_positive*, *loci_negative*,
                            *snps_positive*, and *snps_negative* files for one pair.
  output_root               Directory where pair-specific FUMA inputs are written.
  pair_name                 Optional pair name. Defaults to the first two
                            hyphen-delimited fields of the input directory name,
                            matching the historical cluster script.

This script is path-neutral. To run it as a Slurm array, use the generic helper:

  bash scripts/slurm/submit_array.sh pair_dirs.txt \
    bash scripts/05_fuma/04_prepare_fuma_pair_uploads.sh '$CVP_ARRAY_ITEM' /path/to/fuma_inputs
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

indir=$1
out_root=$2

if [[ $# -eq 3 ]]; then
  pair_name=$3
else
  base_name=$(basename "$indir")
  pair_name=$(printf '%s\n' "$base_name" | cut -d "-" -f1,2)
fi
out_dir="${out_root}/${pair_name}"

mkdir -p "$out_dir"

first_match() {
  local pattern=$1
  find "$indir" -maxdepth 1 -type f -name "$pattern" | sort | head -n 1
}

write_lead_file() {
  local infile=$1
  local outfile=$2

  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == 1 { print "SNP", "CHR", "BP" }
    NR > 1 { print $3, $2, $4 }
  ' "$infile" > "$outfile"
}

write_candidate_file() {
  local infile=$1
  local outfile=$2

  awk -F'\t' 'BEGIN { OFS="\t" }
    NR == 1 { print "CHR", "BP", "SNP", "A1", "A2", "P" }
    NR > 1 { print $2, $5, $6, $9, $10, $11 }
  ' "$infile" > "$outfile"
  gzip -f "$outfile"
}

lead_pos=$(first_match "*${pair_name}*loci_positive*.csv")
lead_neg=$(first_match "*${pair_name}*loci_negative*.csv")
snp_pos=$(first_match "*${pair_name}*snps_positive*.csv")
snp_neg=$(first_match "*${pair_name}*snps_negative*.csv")

if [[ -n "$lead_pos" ]]; then
  write_lead_file "$lead_pos" "${out_dir}/${pair_name}-IndepLeadPos.txt"
fi

if [[ -n "$lead_neg" ]]; then
  write_lead_file "$lead_neg" "${out_dir}/${pair_name}-IndepLeadNeg.txt"
fi

if [[ -n "$snp_pos" ]]; then
  write_candidate_file "$snp_pos" "${out_dir}/${pair_name}-CandidatePos.csv"
fi

if [[ -n "$snp_neg" ]]; then
  write_candidate_file "$snp_neg" "${out_dir}/${pair_name}-CandidateNeg.csv"
fi
