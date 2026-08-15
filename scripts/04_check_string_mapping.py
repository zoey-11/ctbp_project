from pathlib import Path
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

STRING_DIR = PROJECT_ROOT / "data" / "raw" / "string"
OUTPUT_DIR = PROJECT_ROOT / "data" / "interim"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PHYSICAL_FILE = (
    STRING_DIR / "9606.protein.physical.links.detailed.v12.0.txt"
)

ALIASES_FILE = (
    STRING_DIR / "9606.protein.aliases.v12.0.txt"
)

INFO_FILE = (
    STRING_DIR / "9606.protein.info.v12.0.txt"
)


# ============================================================
# Target proteins
# ============================================================

TARGETS = {
    "9606.ENSP00000290921": {
        "gene": "CTBP1",
        "uniprot": "Q13363",
    },
    "9606.ENSP00000311825": {
        "gene": "CTBP2",
        "uniprot": "P56545",
    },
}

target_ids = set(TARGETS.keys())


# ============================================================
# Read STRING protein information
# ============================================================

print("Reading STRING protein information...")

info = pd.read_csv(
    INFO_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "string_id",
        "preferred_name",
        "protein_size",
        "annotation",
    ],
    dtype=str,
)

string_to_gene = dict(
    zip(
        info["string_id"],
        info["preferred_name"],
    )
)

print(
    f"STRING proteins in info file: "
    f"{len(string_to_gene):,}"
)


# ============================================================
# Read STRING aliases
# ============================================================

print("Reading STRING aliases...")

aliases = pd.read_csv(
    ALIASES_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "string_id",
        "alias",
        "source",
    ],
    dtype=str,
)

# Keep explicit UniProt accession mappings.
uniprot_aliases = aliases[
    aliases["source"].isin(
        [
            "UniProt_AC",
            "Ensembl_UniProt",
            "Ensembl_HGNC_uniprot_ids",
        ]
    )
].copy()

# Keep accession-like identifiers.
uniprot_aliases = uniprot_aliases[
    uniprot_aliases["alias"].str.match(
        r"^[A-Z][A-Z0-9]{5,9}$",
        na=False,
    )
]

uniprot_aliases = uniprot_aliases.drop_duplicates(
    ["string_id", "alias"]
)

string_to_uniprot = (
    uniprot_aliases
    .groupby("string_id")["alias"]
    .apply(lambda x: ";".join(sorted(set(x))))
    .to_dict()
)

print(
    f"STRING proteins with UniProt mappings: "
    f"{len(string_to_uniprot):,}"
)


# ============================================================
# Read physical interaction network
# ============================================================

print("Reading STRING physical interaction network...")

results = []

for chunk in pd.read_csv(
    PHYSICAL_FILE,
    sep=r"\s+",
    dtype=str,
    chunksize=500_000,
):

    mask = (
        chunk["protein1"].isin(target_ids)
        | chunk["protein2"].isin(target_ids)
    )

    selected = chunk[mask].copy()

    if len(selected) > 0:
        results.append(selected)

    print(
        f"Processed {len(chunk):,} records | "
        f"Found {len(selected):,} target records"
    )


if not results:
    raise RuntimeError(
        "No CtBP1/CtBP2 interactions found."
    )


df = pd.concat(
    results,
    ignore_index=True,
)


# ============================================================
# Identify target and partner
# ============================================================

def identify_target(row):

    if row["protein1"] in target_ids:
        return row["protein1"]

    return row["protein2"]


df["target_string_id"] = df.apply(
    identify_target,
    axis=1,
)

df["partner_string_id"] = df.apply(
    lambda row:
        row["protein2"]
        if row["protein1"] == row["target_string_id"]
        else row["protein1"],
    axis=1,
)


# ============================================================
# Add annotations
# ============================================================

df["target_gene"] = df["target_string_id"].map(
    string_to_gene
)

df["partner_gene"] = df["partner_string_id"].map(
    string_to_gene
)

df["target_uniprot"] = df["target_string_id"].map(
    lambda x: TARGETS[x]["uniprot"]
)

df["partner_uniprot_accessions"] = df[
    "partner_string_id"
].map(
    lambda x: string_to_uniprot.get(x)
)


# ============================================================
# Remove self-interactions
# ============================================================

df = df[
    df["partner_string_id"]
    != df["target_string_id"]
].copy()


# ============================================================
# Reorder columns
# ============================================================

df = df[
    [
        "target_gene",
        "target_uniprot",
        "target_string_id",
        "partner_gene",
        "partner_string_id",
        "partner_uniprot_accessions",
        "experimental",
        "database",
        "textmining",
        "combined_score",
    ]
]


# ============================================================
# Save interaction dataset
# ============================================================

output_file = OUTPUT_DIR / "string_raw.csv"

df.to_csv(
    output_file,
    index=False,
)


# ============================================================
# Save STRING identifier mapping
# ============================================================

mapping_ids = set(df["partner_string_id"])

mapping = info[
    info["string_id"].isin(mapping_ids)
].copy()

mapping["uniprot_accessions"] = mapping[
    "string_id"
].map(
    string_to_uniprot
)

mapping = mapping[
    [
        "string_id",
        "preferred_name",
        "protein_size",
        "uniprot_accessions",
    ]
]

mapping_file = (
    OUTPUT_DIR
    / "string_protein_mapping.csv"
)

mapping.to_csv(
    mapping_file,
    index=False,
)


# ============================================================
# Summary
# ============================================================

print()
print("=" * 60)
print("STRING PARSING COMPLETE")
print("=" * 60)

print(
    f"Total target interaction records: "
    f"{len(df):,}"
)

print()

for target in ["CTBP1", "CTBP2"]:

    subset = df[
        df["target_gene"] == target
    ]

    print(f"{target}:")
    print(
        f"  Interaction records: "
        f"{len(subset):,}"
    )
    print(
        f"  Unique STRING partners: "
        f"{subset['partner_string_id'].nunique():,}"
    )

    print()

print(
    f"Interaction file:\n{output_file}"
)

print(
    f"\nMapping file:\n{mapping_file}"
)