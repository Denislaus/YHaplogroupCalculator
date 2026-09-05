Y-Predictor
A lightweight Python tool to determine a paternal Y-DNA haplogroup from consumer raw DNA files (23andMe, MyHeritage, FamilyTreeDNA, AncestryDNA).

Runs completely offline on your machine—no telemetry, no external API calls, and zero external dependencies beyond standard Python.

Features
No proprietary packages: Built entirely using Python's standard library (csv, json, argparse, etc.).

Auto-format detection: Handles both tab-delimited files (23andMe) and comma-separated files (MyHeritage, FTDNA) automatically.

Noise filtering: Microarrays often throw false positives on random probes. The script clusters calls across the phylogenetic tree to suppress noise from contradictory branches.

Readable output: Resolves both classic ISOGG tree codes (I2a1b2a1a) and modern SNP notation (I-S17250, N-Z1936).

Prerequisites
Python 3.8 or higher.

A raw DNA data export from a biological male.

Quick Start
1. Clone the repository
Bash
git clone https://github.com/yourusername/y-predictor.git
cd y-predictor
2. Generate the reference database
The script needs a reference index of Y-chromosome mutations mapped to GRCh37/hg19.

Download the curated ISOGG variant file from the yhaplo repository and save it in the project root as isogg.2016.01.04.txt.

Then run the converter:

Bash
python convert_isogg.py
This parses the raw table and outputs isogg_ysnp_db.json. You only need to do this once.

3. Run the predictor
Pass your raw DNA file directly:

Bash
python predict_y_haplogroup.py path/to/raw_dna.txt
Or run it without arguments to get an interactive prompt:

Bash
python predict_y_haplogroup.py
# Enter path to raw DNA file: MyHeritage_raw_dna_data.csv
Example Output
Plaintext
[+] Loading reference database: isogg_ysnp_db.json...
[+] Database indexed: 20,904 reference markers loaded.
[+] Reading sample file: MyHeritage_raw_dna_data.csv
[+] Detected file layout: MyHeritage / CSV
[+] Extracted 7,181 valid Y-chromosome calls.

=================================================================
                Y-HAPLOGROUP PREDICTION RESULT
=================================================================
  Primary Haplogroup (Common Name) :  I-S17250
  Equivalent Terminal Mutations   :  S17250, YP204
  ISOGG Branch Designation        :  I2a1b2a1a
  Major Paternal Clade            :  Haplogroup I
=================================================================

[+] Direct Paternal Lineage Path:
  └── I                      (CTS48, PF3569)
    └── I2a1b                  (M423, CTS176, L178)
      └── I2a1b2                 (L621, S392)
        └── I2a1b2a1               (CTS10228)
          └── I2a1b2a1a              (S17250, YP204)
Optional CLI Flags
Bash
usage: predict_y_haplogroup.py [-h] [--db DB] [input_file]

positional arguments:
  input_file   Path to raw DNA file (.txt, .csv)

options:
  -h, --help   show this help message and exit
  --db DB      Custom path to reference database JSON (defaults to ./isogg_ysnp_db.json)
Supported File Formats
Provider	Build Used	Default File Extension	Delimiter
23andMe (v3, v4, v5)	GRCh37	.txt	Tab (\t)
MyHeritage	GRCh37	.csv	Comma (,)
FTDNA	GRCh37	.csv	Comma (,)
AncestryDNA	GRCh37	.txt	Tab (\t)
Note: Consumer microarrays test a predefined set of SNPs (usually 2,000–8,000 markers on the Y chromosome). This script will find the deepest branch represented on that chip, but it cannot resolve ultra-recent subclades only detectable via next-generation sequencing (WGS or Big Y).

License
MIT