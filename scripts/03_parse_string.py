from pathlib import Path
import re

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
# UniProt accession format
# ============================================================

# Official UniProt accession syntax. The looser pattern
# ^[A-Z][A-Z0-9]{5,9}$ also matches gene symbols such as
# CYP26B1 and HEL161, which the Ensembl_UniProt alias source
# contains alongside genuine accessions.
UNIPROT_AC = re.compile(
    r"^("
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")$"
)

# Alias sources that carry UniProt accessions, in order of
# preference for picking a single canonical accession.
# Ensembl_HGNC_uniprot_ids is HGNC-curated and effectively
# one reviewed accession per gene; UniProt_AC and
# Ensembl_UniProt also list secondary accessions and
# TrEMBL isoform fragments.
CANONICAL_SOURCES = [
    "Ensembl_HGNC_uniprot_ids",
    "UniProt_AC",
    "Ensembl_UniProt",
]


def pick_canonical(accessions):
    """Choose one accession from a set, deterministically.

    Prefers 6-character accessions, which are the legacy
    Swiss-Prot range, over 10-character ones, which are
    predominantly TrEMBL. Ties break alphabetically.
    """

    return sorted(accessions, key=lambda a: (len(a), a))[0]


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

print(f"STRING proteins in info file: {len(string_to_gene):,}")


# ============================================================
# Build STRING -> UniProt mapping
# ============================================================

print("Reading STRING aliases...")

alias_df = pd.read_csv(
    ALIASES_FILE,
    sep="\t",
    comment="#",
    header=None,
    names=["string_id", "alias", "source"],
    dtype=str,
)

# Keep UniProt-derived alias sources.
uniprot_aliases = alias_df[
    alias_df["source"].isin(CANONICAL_SOURCES)
].copy()

# Keep only values that are actually accessions.
rejected = (
    ~uniprot_aliases["alias"].str.match(UNIPROT_AC, na=False)
).sum()

uniprot_aliases = uniprot_aliases[
    uniprot_aliases["alias"].str.match(UNIPROT_AC, na=False)
]

uniprot_aliases = uniprot_aliases.drop_duplicates(
    ["string_id", "alias", "source"]
)

print(f"Non-accession aliases rejected: {rejected:,}")

# Every accession a STRING protein maps to, for provenance.
string_to_all_uniprot = (
    uniprot_aliases
    .groupby("string_id")["alias"]
    .apply(lambda x: sorted(set(x)))
    .to_dict()
)

# One canonical accession per STRING protein, taking the
# first source in CANONICAL_SOURCES that yields a hit.
by_source = {
    source: (
        uniprot_aliases[uniprot_aliases["source"] == source]
        .groupby("string_id")["alias"]
        .apply(set)
        .to_dict()
    )
    for source in CANONICAL_SOURCES
}

string_to_canonical = {}
canonical_source_used = {}

for string_id in string_to_all_uniprot:

    for source in CANONICAL_SOURCES:

        candidates = by_source[source].get(string_id)

        if candidates:
            string_to_canonical[string_id] = pick_canonical(candidates)
            canonical_source_used[string_id] = source
            break

print(
    f"STRING proteins with UniProt mappings: "
    f"{len(string_to_all_uniprot):,}"
)


# ============================================================
# Read STRING physical interaction network
# ============================================================

print("Reading STRING physical interaction network...")

results = []

for chunk in pd.read_csv(
    PHYSICAL_FILE,
    sep=r"\s+",
    dtype=str,
    chunksize=500_000,
):

    # Find interactions where CtBP1 or CtBP2 is involved.
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
        "No CtBP1/CtBP2 interactions found in STRING."
    )


df = pd.concat(results, ignore_index=True)

print(f"Raw target records (both edge directions): {len(df):,}")


# ============================================================
# Orient each edge as target -> partner
# ============================================================

# An edge may name a target in either column. When both
# columns are targets (the CtBP1-CtBP2 edge) protein1 is
# taken as the target, so the pair is retained once in each
# orientation.
p1_is_target = df["protein1"].isin(target_ids)

df["target_string_id"] = df["protein1"].where(
    p1_is_target,
    df["protein2"],
)

df["partner_string_id"] = df["protein2"].where(
    p1_is_target,
    df["protein1"],
)


# ============================================================
# Remove self-interactions
# ============================================================

df = df[
    df["partner_string_id"] != df["target_string_id"]
].copy()


# ============================================================
# Collapse reciprocal edges
# ============================================================

# STRING lists every edge twice, as "A B" and as "B A".
# Both copies match the target filter, so each interaction
# arrives here duplicated.

SCORE_COLUMNS = [
    "experimental",
    "database",
    "textmining",
    "combined_score",
]

for column in SCORE_COLUMNS:
    df[column] = pd.to_numeric(df[column])

pair_key = ["target_string_id", "partner_string_id"]

# The two directions should carry identical scores. Verify
# rather than assume, so a silent STRING format change does
# not quietly discard real differences.
inconsistent = (
    df.groupby(pair_key)[SCORE_COLUMNS]
    .nunique()
    .gt(1)
    .any(axis=1)
    .sum()
)

if inconsistent:
    raise RuntimeError(
        f"{inconsistent} target-partner pairs have conflicting "
        f"scores between edge directions; deduplication would "
        f"lose information."
    )

before = len(df)

df = df.drop_duplicates(pair_key).copy()

print(
    f"Collapsed reciprocal edges: {before:,} -> {len(df):,} "
    f"unique target-partner pairs"
)


# ============================================================
# Annotate
# ============================================================

df["target_gene"] = df["target_string_id"].map(
    lambda x: TARGETS[x]["gene"]
)

df["target_uniprot"] = df["target_string_id"].map(
    lambda x: TARGETS[x]["uniprot"]
)

df["partner_gene"] = df["partner_string_id"].map(
    string_to_gene
)

df["partner_uniprot"] = df["partner_string_id"].map(
    string_to_canonical
)

df["partner_uniprot_all"] = df["partner_string_id"].map(
    lambda x: ";".join(string_to_all_uniprot.get(x, [])) or None
)


# ============================================================
# Select final columns
# ============================================================

df = df[
    [
        "target_gene",
        "target_uniprot",
        "target_string_id",
        "partner_gene",
        "partner_string_id",
        "partner_uniprot",
        "partner_uniprot_all",
    ]
    + SCORE_COLUMNS
]

df = df.sort_values(
    ["target_gene", "combined_score", "partner_gene"],
    ascending=[True, False, True],
).reset_index(drop=True)


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

mapping = info[
    info["string_id"].isin(set(df["partner_string_id"]))
].copy()

mapping["uniprot_canonical"] = mapping["string_id"].map(
    string_to_canonical
)

mapping["uniprot_accessions"] = mapping["string_id"].map(
    lambda x: ";".join(string_to_all_uniprot.get(x, [])) or None
)

mapping["canonical_source"] = mapping["string_id"].map(
    canonical_source_used
)

mapping = mapping[
    [
        "string_id",
        "preferred_name",
        "protein_size",
        "uniprot_canonical",
        "uniprot_accessions",
        "canonical_source",
    ]
].sort_values("preferred_name")

mapping_file = OUTPUT_DIR / "string_protein_mapping.csv"

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

print(f"Unique target-partner pairs: {len(df):,}")

print()

for gene in ["CTBP1", "CTBP2"]:

    subset = df[df["target_gene"] == gene]

    mapped = subset[subset["partner_uniprot"].notna()]

    print(f"{gene}:")
    print(f"  Interactions: {len(subset):,}")
    print(
        f"  Unique STRING partners: "
        f"{subset['partner_string_id'].nunique():,}"
    )
    print(
        f"  Partners with canonical UniProt: "
        f"{mapped['partner_string_id'].nunique():,}"
    )
    print()

ctbp1_partners = set(
    df.loc[df["target_gene"] == "CTBP1", "partner_string_id"]
)

ctbp2_partners = set(
    df.loc[df["target_gene"] == "CTBP2", "partner_string_id"]
)

print(
    f"Partners shared by CTBP1 and CTBP2: "
    f"{len(ctbp1_partners & ctbp2_partners):,}"
)

print(
    f"Distinct partners overall: "
    f"{len(ctbp1_partners | ctbp2_partners):,}"
)

print()
print("Evidence channels:")

has_experimental = df["experimental"] > 0
has_database = df["database"] > 0

print(
    f"  experimental > 0: {has_experimental.sum():,}"
)
print(
    f"  database > 0:     {has_database.sum():,}"
)
print(
    f"  textmining only:  "
    f"{(~has_experimental & ~has_database).sum():,}"
)

print()

fallback = mapping[
    mapping["canonical_source"] != "Ensembl_HGNC_uniprot_ids"
]

print(
    f"Partners whose canonical accession came from a "
    f"fallback source: {len(fallback):,}"
)

if len(fallback) > 0:
    for _, row in fallback.iterrows():
        print(
            f"  {row['preferred_name']} -> "
            f"{row['uniprot_canonical']} "
            f"({row['canonical_source']})"
        )

print()
print(f"Interaction file:\n{output_file}")
print(f"\nMapping file:\n{mapping_file}")
