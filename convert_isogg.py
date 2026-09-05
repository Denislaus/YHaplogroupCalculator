#!/usr/bin/env python3
"""
convert_isogg.py
Universal converter for ISOGG Y-SNP master tables (TSV/CSV/TXT) into a JSON database.
Supports multiple ISOGG release formats and column arrangements.
"""

import argparse
import csv
import glob
import json
import os
import re
import sys

DEFAULT_OUTPUT_NAME = "isogg_ysnp_db.json"


def find_input_file(custom_path=None):
    """Locates the raw ISOGG variant file via arguments or local auto-discovery."""
    if custom_path and os.path.exists(custom_path):
        return custom_path

    # Check for common patterns in current directory or script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    search_dirs = [os.getcwd(), script_dir]

    patterns = [
        "isogg*.txt",
        "isogg*.tsv",
        "isogg*.csv",
        "*variants*.txt",
        "*ysnp*.txt"
    ]

    for d in search_dirs:
        for pat in patterns:
            matches = glob.glob(os.path.join(d, pat))
            if matches:
                # Prioritize non-json files
                valid_matches = [m for m in matches if not m.endswith(".json")]
                if valid_matches:
                    return valid_matches[0]

    return None


def detect_delimiter(sample_lines):
    """Detects whether raw data is tab-delimited or comma-delimited."""
    sample_text = "".join(sample_lines)
    return "\t" if sample_text.count("\t") >= sample_text.count(",") else ","


def parse_header_indices(header_row):
    """
    Dynamically resolves column indexes based on header names,
    accommodating changes in ISOGG formats over time.
    """
    col_map = {"name": -1, "hg": -1, "pos": -1, "mut": -1}

    for idx, col in enumerate(header_row):
        col_clean = col.strip().lower()
        if "snp" in col_clean or "marker" in col_clean or "name" in col_clean:
            if col_map["name"] == -1 and "other" not in col_clean:
                col_map["name"] = idx
        elif "haplogroup" in col_clean or "clade" in col_clean:
            col_map["hg"] = idx
        elif "position" in col_clean or "grch37" in col_clean or "b37" in col_clean or "build 37" in col_clean:
            col_map["pos"] = idx
        elif "mutation" in col_clean or "change" in col_clean or "allele" in col_clean:
            col_map["mut"] = idx

    # Fallback default positions for standard 6-column ISOGG table
    if col_map["name"] == -1:
        col_map["name"] = 0
    if col_map["hg"] == -1:
        col_map["hg"] = 1
    if col_map["pos"] == -1:
        col_map["pos"] = 4
    if col_map["mut"] == -1:
        col_map["mut"] = 5

    return col_map


def convert(input_path, output_path):
    if not os.path.exists(input_path):
        sys.exit(f"[-] Error: Could not find input file: {input_path}")

    print(f"[+] Reading input file: {input_path}")

    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        sample_lines = [f.readline() for _ in range(25)]

    delimiter = detect_delimiter(sample_lines)
    sep_label = "Tab (\\t)" if delimiter == "\t" else "Comma (,)"
    print(f"[+] Detected delimiter: {sep_label}")

    db = []
    skipped = 0
    seen_snps = set()

    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter=delimiter)

        header = None
        col_map = None

        for row_idx, row in enumerate(reader):
            if not row or not any(row):
                continue

            row_str = " ".join(row).lower()

            # Identify and parse the header line
            if header is None and ("haplogroup" in row_str or "snp" in row_str):
                header = row
                col_map = parse_header_indices(header)
                continue

            # Fallback if no explicit header row exists
            if col_map is None:
                col_map = parse_header_indices(["SNP", "Haplogroup", "AltNames", "RefSNP", "Position", "Mutation"])

            # Ensure row has enough columns for indexed fields
            max_idx = max(col_map.values())
            if len(row) <= max_idx:
                skipped += 1
                continue

            name = row[col_map["name"]].strip()
            hg = row[col_map["hg"]].strip()
            pos_str = row[col_map["pos"]].strip()
            mut_str = row[col_map["mut"]].strip()

            # Clean haplogroup (strip footnotes like " (Notes)")
            hg = re.sub(r"\s*\(.*?\)", "", hg).strip()

            # Must have a valid numeric coordinate
            pos_digits = re.sub(r"[^\d]", "", pos_str)
            if not pos_digits:
                skipped += 1
                continue
            pos = int(pos_digits)

            # Match single base mutations: "G->C", "G -> C", "G to C", "G/C"
            match = re.search(r"([ACGT])\s*(?:->|>|to|/)\s*([ACGT])", mut_str, re.IGNORECASE)
            if not match:
                skipped += 1
                continue

            anc = match.group(1).upper()
            der = match.group(2).upper()

            # Deduplicate by unique (name, pos, derived)
            unique_key = (name, pos, der)
            if unique_key in seen_snps:
                continue
            seen_snps.add(unique_key)

            db.append({
                "name": name,
                "pos": pos,
                "ancestral": anc,
                "derived": der,
                "haplogroup": hg
            })

    print(f"[+] Converted {len(db):,} valid SNPs (skipped {skipped:,} non-SNP/indel rows).")

    if not db:
        print("[-] Error: 0 records converted. Please verify file format and columns.")
        return

    # Write output
    out_dir = os.path.dirname(os.path.abspath(output_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as out:
        json.dump(db, out, indent=2)

    print(f"[+] Saved database to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert ISOGG Y-SNP master tables into a structured JSON reference database."
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        help="Path to raw ISOGG variant file (.txt, .tsv, .csv). Auto-discovers if omitted."
    )
    parser.add_argument(
        "-o", "--output",
        default=DEFAULT_OUTPUT_NAME,
        help=f"Path to output JSON file (default: {DEFAULT_OUTPUT_NAME})"
    )
    args = parser.parse_args()

    input_file = args.input_file
    if not input_file:
        input_file = find_input_file()
        if not input_file:
            input_file = input("Enter path to raw ISOGG table file: ").strip(' "')

    convert(input_file, args.output)


if __name__ == "__main__":
    main()