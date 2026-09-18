#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def bh(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    m = len(pvalues)

    order = np.argsort(pvalues)
    ranked = pvalues[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.minimum(ranked, 1.0)

    out = np.empty(m)
    out[order] = ranked
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True,
                    help="Text file with one Ensembl gene ID per line")
    ap.add_argument("--background", required=True,
                    help="FUMA ENSG.genes.txt")
    ap.add_argument("--genesets", required=True,
                    help="Directory containing historical FUMA *.gmt files")
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--adj-p", type=float, default=0.05)
    args = ap.parse_args()

    gene_file = Path(args.genes)
    geneset_dir = Path(args.genesets)
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # FUMA background / Ensembl -> Entrez mapping
    # --------------------------------------------------
    bg = pd.read_csv(args.background, sep="\t", dtype=str)

    bg["entrez_clean"] = (
        pd.to_numeric(bg["entrezID"], errors="coerce")
        .astype("Int64")
        .astype(str)
    )

    valid_bg = bg.loc[
        bg["entrez_clean"].notna()
        & (bg["entrez_clean"] != "<NA>")
    ].copy()

    background = set(valid_bg["entrez_clean"])

    # --------------------------------------------------
    # Input Ensembl IDs
    # --------------------------------------------------
    input_ensg = {
        x.strip()
        for x in gene_file.read_text().splitlines()
        if x.strip()
    }

    input_map = valid_bg[valid_bg["ensembl_gene_id"].isin(input_ensg)].copy()

    input_entrez = set(input_map["entrez_clean"]) & background

    # Entrez -> all input symbols, preserving FUMA-like multiplicity
    entrez_symbols = {}

    for _, r in input_map.iterrows():
        entrez = r["entrez_clean"]

        symbol = r.get("external_gene_name")
        if pd.isna(symbol) or not str(symbol).strip():
            symbol = r.get("hgnc_symbol")

        if pd.notna(symbol) and str(symbol).strip():
            entrez_symbols.setdefault(entrez, [])
            s = str(symbol).strip()
            if s not in entrez_symbols[entrez]:
                entrez_symbols[entrez].append(s)

    M = len(background)
    K = len(input_entrez)

    print("Input Ensembl IDs:", len(input_ensg))
    print("Recognised Ensembl IDs:", input_map["ensembl_gene_id"].nunique())
    print("Unique input Entrez IDs:", K)
    print("Background Entrez IDs:", M)

    all_sig = []

    # --------------------------------------------------
    # Reconstruct each FUMA GENE2FUNC category
    # --------------------------------------------------
    for gmt_path in sorted(geneset_dir.glob("*.gmt")):

        category = gmt_path.stem
        rows = []

        with gmt_path.open() as f:
            for line_no, line in enumerate(f):
                fields = line.rstrip("\n").split("\t")

                if len(fields) < 3:
                    continue

                gene_set = fields[0]
                link = fields[1]

                raw_genes = {x for x in fields[2:] if x}
                genes = raw_genes & background

                n = len(genes)
                overlap_ids = input_entrez & genes
                x = len(overlap_ids)

                p = (
                    hypergeom.sf(x - 1, M, n, K)
                    if x > 0 else 1.0
                )

                symbols = []

                for entrez in overlap_ids:
                    for symbol in entrez_symbols.get(entrez, []):
                        if symbol not in symbols:
                            symbols.append(symbol)

                rows.append({
                    "_line": line_no,
                    "Category": category,
                    "GeneSet": gene_set,
                    "N_genes": n,
                    "N_overlap": x,
                    "p": p,
                    "genes": ":".join(symbols),
                    "link": link,
                })

        df = pd.DataFrame(rows)

        if df.empty:
            continue

        # FUMA behaviour confirmed against historical output:
        # BH correction is across ALL sets in the category.
        df["adjP"] = bh(df["p"].values)

        sig = df[
            (df["N_overlap"] >= args.min_overlap)
            & (df["adjP"] < args.adj_p)
        ].copy()

        if not sig.empty:
            sig = sig[
                [
                    "Category",
                    "GeneSet",
                    "N_genes",
                    "N_overlap",
                    "p",
                    "adjP",
                    "genes",
                    "link",
                ]
            ]

            all_sig.append(sig)

        print(
            f"{category}: "
            f"{len(df)} tested, "
            f"{len(sig)} significant"
        )

    if all_sig:
        result = pd.concat(all_sig, ignore_index=True)
    else:
        result = pd.DataFrame(
            columns=[
                "Category", "GeneSet", "N_genes", "N_overlap",
                "p", "adjP", "genes", "link"
            ]
        )

    outfile = outdir / "GS.txt"
    result.to_csv(outfile, sep="\t", index=False)

    # Also retain mapping information for provenance
    input_map[
        [
            "ensembl_gene_id",
            "entrezID",
            "external_gene_name",
            "hgnc_symbol",
        ]
    ].drop_duplicates().to_csv(
        outdir / "gene_mapping.tsv",
        sep="\t",
        index=False,
    )

    print("\nTotal significant gene sets:", len(result))
    print("Output:", outfile)


if __name__ == "__main__":
    main()
