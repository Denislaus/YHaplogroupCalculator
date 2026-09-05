#!/usr/bin/env python3
"""
predict_y_haplogroup.py
An open-source, client-side Y-DNA haplogroup predictor for raw consumer DNA data.
Supports 23andMe, MyHeritage, FamilyTreeDNA, and AncestryDNA formats (GRCh37/hg19).
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

DEFAULT_DB_NAME = "isogg_ysnp_db.json"


def find_database(custom_path=None):
    """Locates the reference JSON database relative to script execution."""
    if custom_path and os.path.exists(custom_path):
        return custom_path

    # Check script's parent directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_path = os.path.join(script_dir, DEFAULT_DB_NAME)
    if os.path.exists(local_path):
        return local_path

    # Check current working directory
    if os.path.exists(DEFAULT_DB_NAME):
        return DEFAULT_DB_NAME

    return None


def load_snp_reference(db_path):
    """Loads SNP reference mutations keyed by physical GRCh37 position."""
    if not db_path or not os.path.exists(db_path):
        sys.exit(
            f"[-] Error: Could not locate '{DEFAULT_DB_NAME}'.\n"
            f"    Please place '{DEFAULT_DB_NAME}' in the script directory or pass --db <path>."
        )

    print(f"[+] Loading reference database: {db_path}...")
    with open(db_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    db = defaultdict(list)
    for row in data:
        db[int(row["pos"])].append(row)

    print(f"[+] Database indexed: {len(data):,} reference markers loaded.")
    return db


def detect_delimiter(sample_lines):
    """Detects whether raw data is tab-delimited or comma-delimited."""
    non_comment = [line for line in sample_lines if line.strip() and not line.startswith("#")]
    sample_text = "".join(non_comment)
    return "\t" if sample_text.count("\t") > sample_text.count(",") else ","


def parse_raw_y(file_path):
    """
    Parses consumer microarray raw DNA files.
    Extracts valid non-heterozygous, non-blank Y chromosome calls.
    """
    if not os.path.exists(file_path):
        sys.exit(f"[-] Error: Input file not found: {file_path}")

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        sample_lines = [f.readline() for _ in range(40)]

    delimiter = detect_delimiter(sample_lines)
    fmt_label = "23andMe / TSV" if delimiter == "\t" else "MyHeritage / CSV"
    print(f"[+] Reading sample file: {file_path}")
    print(f"[+] Detected file layout: {fmt_label}")

    y_calls = {}
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 4:
                continue

            # Standard layout: [RSID, CHROMOSOME, POSITION, GENOTYPE/RESULT]
            rsid = row[0].strip(' "')
            chrom = row[1].strip(' "').upper()
            pos_str = row[2].strip(' "')
            call = row[3].strip(' "').upper().replace("-", "").replace("?", "")

            if chrom in ("Y", "CHRY", "24"):
                try:
                    pos = int(pos_str)
                    # Hemizygous filter: ignore no-calls and indels
                    if call and call not in ("N", "DD", "II", "DI"):
                        y_calls[pos] = (rsid, call[0])
                except ValueError:
                    continue

    print(f"[+] Extracted {len(y_calls):,} valid Y-chromosome calls.")
    return y_calls


def evaluate_calls(y_calls, snp_db):
    """Matches sample alleles against reference ancestral and derived states."""
    derived_by_hg = defaultdict(list)
    root_derived_counts = defaultdict(int)

    for pos, (rsid, allele) in y_calls.items():
        if pos not in snp_db:
            continue

        for ref in snp_db[pos]:
            hg = ref["haplogroup"]
            root = hg[0] if hg else ""

            if allele == ref["derived"]:
                derived_by_hg[hg].append(ref)
                if root.isalpha():
                    root_derived_counts[root] += 1

    return derived_by_hg, root_derived_counts


def resolve_lineage(derived_by_hg, root_counts):
    """
    Identifies the true root clade and resolves the deepest valid branch,
    filtering out probe cross-hybridization noise.
    """
    if not root_counts:
        return None, None, []

    # Identify primary root family with the highest density of derived markers
    dominant_root = max(root_counts.items(), key=lambda x: x[1])[0]

    # Filter candidate branches belonging to the confirmed paternal tree
    candidate_hgs = [hg for hg in derived_by_hg.keys() if hg.startswith(dominant_root)]

    if not candidate_hgs:
        return dominant_root, dominant_root, []

    # Sort primarily by phylogenetic depth (ISOGG length), secondarily by marker count
    candidate_hgs.sort(key=lambda hg: (len(hg), len(derived_by_hg[hg])), reverse=True)

    terminal_isogg = candidate_hgs[0]

    # Walk path: sort all confirmed branches from root down to terminal
    path_branches = sorted(candidate_hgs, key=lambda hg: len(hg))

    return dominant_root, terminal_isogg, path_branches


def main():
    parser = argparse.ArgumentParser(
        description="Predict Y-DNA haplogroup from raw consumer DNA files (23andMe, MyHeritage, etc.)."
    )
    parser.add_argument("input_file", nargs="?", help="Path to raw DNA file (.txt or .csv)")
    parser.add_argument("--db", default=None, help=f"Path to reference database (default: {DEFAULT_DB_NAME})")
    args = parser.parse_args()

    # Prompt interactively if file is not provided via CLI
    input_path = args.input_file
    if not input_path:
        input_path = input("Enter path to raw DNA file: ").strip(' "')

    db_path = find_database(args.db)
    snp_db = load_snp_reference(db_path)
    y_calls = parse_raw_y(input_path)

    if not y_calls:
        print("[-] No Y-DNA markers found. Check that the file belongs to a biological male and contains Y data.")
        return

    derived_by_hg, root_counts = evaluate_calls(y_calls, snp_db)
    dominant_root, terminal_isogg, path_branches = resolve_lineage(derived_by_hg, root_counts)

    if not dominant_root:
        print("[-] Could not reliably resolve a paternal lineage from the provided markers.")
        return

    # Find the terminal defining SNP name
    terminal_snps = derived_by_hg[terminal_isogg]
    terminal_snp_names = list({s["name"] for s in terminal_snps})
    commercial_label = f"{dominant_root}-{terminal_snp_names[0]}" if terminal_snp_names else terminal_isogg

    # Display clear, user-friendly results
    print("\n" + "=" * 65)
    print("                Y-HAPLOGROUP PREDICTION RESULT")
    print("=" * 65)
    print(f"  Primary Haplogroup (Common Name) :  \033[1m{commercial_label}\033[0m")
    if len(terminal_snp_names) > 1:
        print(f"  Equivalent Terminal Mutations   :  {', '.join(terminal_snp_names)}")
    print(f"  ISOGG Branch Designation        :  {terminal_isogg}")
    print(f"  Major Paternal Clade            :  Haplogroup {dominant_root}")
    print("=" * 65)

    # Lineage step-by-step path
    print("\n[+] Direct Paternal Lineage Path:")
    seen_steps = set()
    step_num = 1
    for hg in path_branches:
        if hg not in seen_steps:
            snps_on_branch = list({s["name"] for s in derived_by_hg[hg]})
            snp_summary = f"({', '.join(snps_on_branch[:3])})" if snps_on_branch else ""
            indent = "  " * step_num
            print(f"{indent}└── {hg:<22} {snp_summary}")
            seen_steps.add(hg)
            step_num += 1

    # Key supporting markers table
    print("\n[+] Supporting Derived SNPs Confirmed in Sample:")
    print(f"  {'Haplogroup':<22} {'SNP':<12} {'Pos (GRCh37)':<15} {'Derived Call'}")
    print("  " + "-" * 60)

    dominant_snps = []
    for hg in path_branches:
        dominant_snps.extend(derived_by_hg[hg])

    # Show deepest markers first
    seen_markers = set()
    count = 0
    for snp in sorted(dominant_snps, key=lambda x: len(x["haplogroup"]), reverse=True):
        marker_id = (snp["haplogroup"], snp["name"])
        if marker_id not in seen_markers:
            print(f"  {snp['haplogroup']:<22} {snp['name']:<12} {snp['pos']:<15} {snp['derived']}")
            seen_markers.add(marker_id)
            count += 1
            if count >= 15:
                break
    print("=" * 65 + "\n")


if __name__ == "__main__":
    main()