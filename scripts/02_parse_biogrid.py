from pathlib import Path
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "biogrid"
    / "BIOGRID-ORGANISM-Homo_sapiens-5.0.255.tab3.txt"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "interim"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Target proteins
# ============================================================

TARGET_UNIPROT = {
    "Q13363": "CTBP1",
    "P56545": "CTBP2",
}

TARGET_SYMBOLS = {
    "CTBP1",
    "CTBP2",
}


# ============================================================
# Read BioGRID column names directly from the file
# ============================================================

with open(INPUT_FILE, "r", encoding="utf-8") as f:
    header = f.readline().rstrip("\n").lstrip("#").split("\t")


print(f"BioGRID columns detected: {len(header)}")


# ============================================================
# Parse in chunks
# ============================================================

chunks = []

for chunk in pd.read_csv(
    INPUT_FILE,
    sep="\t",
    skiprows=1,
    names=header,
    dtype=str,
    chunksize=100_000,
    low_memory=False,
):

    # Make sure we're only working with human-human interactions
    human_mask = (
        (chunk["Organism ID Interactor A"] == "9606")
        & (chunk["Organism ID Interactor B"] == "9606")
    )

    chunk = chunk[human_mask].copy()

    # Target based on official gene symbol
    symbol_mask = (
        chunk["Official Symbol Interactor A"].isin(TARGET_SYMBOLS)
        | chunk["Official Symbol Interactor B"].isin(TARGET_SYMBOLS)
    )

    # Target based on UniProt accession
    uniprot_mask = (
        chunk["SWISS-PROT Accessions Interactor A"]
        .fillna("")
        .str.contains(
            r"(^|\|)(Q13363|P56545)(\||$)",
            regex=True,
        )
        |
        chunk["SWISS-PROT Accessions Interactor B"]
        .fillna("")
        .str.contains(
            r"(^|\|)(Q13363|P56545)(\||$)",
            regex=True,
        )
    )

    target_mask = symbol_mask | uniprot_mask

    selected = chunk[target_mask].copy()

    if len(selected) > 0:
        chunks.append(selected)

    print(
        f"Processed {len(chunk):,} human records | "
        f"Found {len(selected):,} target records"
    )


# ============================================================
# Combine target interactions
# ============================================================

if chunks:
    df = pd.concat(chunks, ignore_index=True)
else:
    raise RuntimeError("No CtBP1/CtBP2 interactions were found.")


# ============================================================
# Select useful columns
# ============================================================

df = df[
    [
        "BioGRID Interaction ID",
        "Entrez Gene Interactor A",
        "Entrez Gene Interactor B",
        "Official Symbol Interactor A",
        "Official Symbol Interactor B",
        "Synonyms Interactor A",
        "Synonyms Interactor B",
        "Experimental System",
        "Experimental System Type",
        "Author",
        "Publication Source",
        "Organism ID Interactor A",
        "Organism ID Interactor B",
        "Throughput",
        "Score",
        "Modification",
        "Qualifications",
        "Source Database",
        "SWISS-PROT Accessions Interactor A",
        "SWISS-PROT Accessions Interactor B",
    ]
]


# ============================================================
# Save
# ============================================================

output_file = OUTPUT_DIR / "biogrid_raw.csv"

df.to_csv(
    output_file,
    index=False,
)


# ============================================================
# Summary
# ============================================================

print()
print("=" * 60)
print("BIOGRID PARSING COMPLETE")
print("=" * 60)

print(f"Total CtBP interaction records: {len(df):,}")

print()
print("Records by target:")

for target in TARGET_SYMBOLS:
    mask = (
        (df["Official Symbol Interactor A"] == target)
        | (df["Official Symbol Interactor B"] == target)
    )

    print(f"  {target}: {mask.sum():,}")

print()
print(f"Saved to:")
print(output_file)