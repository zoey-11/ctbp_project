from pathlib import Path
import re

import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INTERIM_DIR = PROJECT_ROOT / "data" / "interim"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

INTACT_FILE = INTERIM_DIR / "intact_raw.csv"
BIOGRID_FILE = INTERIM_DIR / "biogrid_raw.csv"
STRING_FILE = INTERIM_DIR / "string_raw.csv"
STRING_MAPPING_FILE = INTERIM_DIR / "string_protein_mapping.csv"


# ============================================================
# Targets
# ============================================================

TARGETS = {
    "Q13363": {
        "gene": "CTBP1",
        "string_id": "9606.ENSP00000290921",
    },
    "P56545": {
        "gene": "CTBP2",
        "string_id": "9606.ENSP00000311825",
    },
}

TARGET_UNIPROTS = set(TARGETS.keys())


# ============================================================
# Helper functions
# ============================================================

def clean_pmid(value):
    """
    Convert a PubMed field to a clean PMID.

    IntAct contains values such as:
        31041561
        unassigned1312

    Non-numeric/unassigned values are treated as missing.
    """
    if pd.isna(value):
        return pd.NA

    value = str(value).strip()

    if value.isdigit():
        return value

    return pd.NA


def split_biogrid_accessions(value):
    """
    Split BioGRID's SWISS-PROT accession field.

    Examples:
        Q13363
        Q13465|Q03112
        Q8NFS9|Q06430|Q8N0V5
    """
    if pd.isna(value):
        return []

    value = str(value).strip()

    if not value or value == "-":
        return []

    return [
        accession.strip()
        for accession in value.split("|")
        if accession.strip()
    ]


def canonical_from_biogrid(value):
    """
    Resolve a BioGRID accession field to a deterministic
    accession.

    BioGRID normally supplies one Swiss-Prot accession.
    When multiple are supplied, prefer the shortest valid
    accession and break ties alphabetically.

    IMPORTANT:
    This is a deterministic harmonization rule, not a claim
    that the first BioGRID accession is biologically superior.
    """
    accessions = split_biogrid_accessions(value)

    if not accessions:
        return pd.NA

    return sorted(set(accessions), key=lambda x: (len(x), x))[0]


def canonical_pair(target, partner):
    """
    Return a target-partner pair only if it is valid.

    CTBP1 and CTBP2 remain separate targets.
    """
    if pd.isna(target) or pd.isna(partner):
        return None

    target = str(target).strip()
    partner = str(partner).strip()

    if target not in TARGET_UNIPROTS:
        return None

    if partner == target:
        return None

    return target, partner


# ============================================================
# 1. Load STRING mapping
# ============================================================

print("Reading STRING protein mapping...")

string_mapping = pd.read_csv(
    STRING_MAPPING_FILE,
    dtype=str,
)

string_to_canonical = dict(
    zip(
        string_mapping["string_id"],
        string_mapping["uniprot_canonical"],
    )
)

string_to_gene = dict(
    zip(
        string_mapping["string_id"],
        string_mapping["preferred_name"],
    )
)

print(
    f"STRING proteins with canonical UniProt: "
    f"{len(string_to_canonical):,}"
)


# ============================================================
# 2. Build evidence records — IntAct
# ============================================================

print("\nReading IntAct...")

intact = pd.read_csv(
    INTACT_FILE,
    dtype=str,
)

intact["intact_miscore"] = pd.to_numeric(
    intact["intact_miscore"],
    errors="coerce",
)

intact["pubmed_clean"] = intact["pubmed"].apply(clean_pmid)

intact_records = []

for _, row in intact.iterrows():

    protein_a = str(row["proteinA"]).strip()
    protein_b = str(row["proteinB"]).strip()

    # Determine which side is CTBP1/CTBP2.
    if protein_a in TARGET_UNIPROTS:
        target = protein_a
        partner = protein_b
    elif protein_b in TARGET_UNIPROTS:
        target = protein_b
        partner = protein_a
    else:
        continue

    pair = canonical_pair(target, partner)

    if pair is None:
        continue

    target, partner = pair

    intact_records.append(
        {
            "source": "IntAct",
            "target_uniprot": target,
            "partner_uniprot": partner,
            "partner_gene": pd.NA,
            "pmid": row["pubmed_clean"],
            "method": row["detection_method"],
            "method_class": pd.NA,
            "intact_miscore": row["intact_miscore"],
            "biogrid_score": pd.NA,
            "biogrid_interaction_id": pd.NA,
            "string_experimental": pd.NA,
            "string_database": pd.NA,
            "string_textmining": pd.NA,
            "string_combined_score": pd.NA,
            "original_partner_id": partner,
        }
    )

intact_evidence = pd.DataFrame(intact_records)

print(
    f"IntAct target-partner evidence records: "
    f"{len(intact_evidence):,}"
)


# ============================================================
# 3. Build evidence records — BioGRID
# ============================================================

print("\nReading BioGRID...")

biogrid = pd.read_csv(
    BIOGRID_FILE,
    dtype=str,
)

biogrid_records = []

for _, row in biogrid.iterrows():

    symbol_a = str(
        row["Official Symbol Interactor A"]
    ).strip()

    symbol_b = str(
        row["Official Symbol Interactor B"]
    ).strip()

    accession_a = canonical_from_biogrid(
        row["SWISS-PROT Accessions Interactor A"]
    )

    accession_b = canonical_from_biogrid(
        row["SWISS-PROT Accessions Interactor B"]
    )

    # Determine target using the canonical accession.
    if accession_a in TARGET_UNIPROTS:
        target = accession_a
        partner = accession_b
        partner_gene = symbol_b

    elif accession_b in TARGET_UNIPROTS:
        target = accession_b
        partner = accession_a
        partner_gene = symbol_a

    else:
        continue

    pair = canonical_pair(target, partner)

    if pair is None:
        continue

    target, partner = pair

    biogrid_records.append(
        {
            "source": "BioGRID",
            "target_uniprot": target,
            "partner_uniprot": partner,
            "partner_gene": partner_gene,
            "pmid": clean_pmid(
                row["Publication Source"]
                .replace("PUBMED:", "")
                if pd.notna(row["Publication Source"])
                else pd.NA
            ),
            "method": row["Experimental System"],
            "method_class": pd.NA,
            "intact_miscore": pd.NA,
            "biogrid_score": row["Score"],
            "biogrid_interaction_id": row[
                "BioGRID Interaction ID"
            ],
            "string_experimental": pd.NA,
            "string_database": pd.NA,
            "string_textmining": pd.NA,
            "string_combined_score": pd.NA,
            "original_partner_id": (
                row["Official Symbol Interactor B"]
                if accession_a in TARGET_UNIPROTS
                else row["Official Symbol Interactor A"]
            ),
        }
    )

biogrid_evidence = pd.DataFrame(biogrid_records)

print(
    f"BioGRID target-partner evidence records: "
    f"{len(biogrid_evidence):,}"
)


# ============================================================
# 4. Build evidence records — STRING
# ============================================================

print("\nReading STRING...")

string = pd.read_csv(
    STRING_FILE,
    dtype=str,
)

for column in [
    "experimental",
    "database",
    "textmining",
    "combined_score",
]:
    string[column] = pd.to_numeric(
        string[column],
        errors="coerce",
    )

string_records = []

for _, row in string.iterrows():

    target = str(row["target_uniprot"]).strip()
    partner = str(row["partner_uniprot"]).strip()

    pair = canonical_pair(target, partner)

    if pair is None:
        continue

    target, partner = pair

    string_records.append(
        {
            "source": "STRING",
            "target_uniprot": target,
            "partner_uniprot": partner,
            "partner_gene": row["partner_gene"],
            "pmid": pd.NA,
            "method": pd.NA,
            "method_class": pd.NA,
            "intact_miscore": pd.NA,
            "biogrid_score": pd.NA,
            "biogrid_interaction_id": pd.NA,
            "string_experimental": row["experimental"],
            "string_database": row["database"],
            "string_textmining": row["textmining"],
            "string_combined_score": row["combined_score"],
            "original_partner_id": row[
                "partner_string_id"
            ],
        }
    )

string_evidence = pd.DataFrame(string_records)

print(
    f"STRING target-partner evidence records: "
    f"{len(string_evidence):,}"
)


# ============================================================
# 5. Combine evidence records
# ============================================================

print("\nCombining evidence records...")

evidence = pd.concat(
    [
        intact_evidence,
        biogrid_evidence,
        string_evidence,
    ],
    ignore_index=True,
)


# ============================================================
# 6. Deduplicate identical evidence records
# ============================================================

before = len(evidence)

evidence = evidence.drop_duplicates(
    subset=[
        "source",
        "target_uniprot",
        "partner_uniprot",
        "pmid",
        "method",
        "biogrid_interaction_id",
        "original_partner_id",
    ]
).copy()

print(
    f"Evidence records after exact deduplication: "
    f"{before:,} -> {len(evidence):,}"
)


# ============================================================
# 7. Add target genes
# ============================================================

evidence["target_gene"] = evidence[
    "target_uniprot"
].map(
    lambda x: TARGETS[x]["gene"]
)


# ============================================================
# 8. Save evidence records
# ============================================================

evidence_columns = [
    "source",
    "target_gene",
    "target_uniprot",
    "partner_gene",
    "partner_uniprot",
    "original_partner_id",
    "pmid",
    "method",
    "method_class",
    "intact_miscore",
    "biogrid_score",
    "biogrid_interaction_id",
    "string_experimental",
    "string_database",
    "string_textmining",
    "string_combined_score",
]

evidence = evidence[evidence_columns]

evidence_file = (
    PROCESSED_DIR / "evidence_records.csv"
)

evidence.to_csv(
    evidence_file,
    index=False,
)

print(f"\nSaved evidence records to:")
print(evidence_file)


# ============================================================
# 9. Aggregate partner-level evidence
# ============================================================

print("\nAggregating partner evidence...")

group_columns = [
    "target_uniprot",
    "partner_uniprot",
]

rows = []

for (target, partner), group in evidence.groupby(
    group_columns,
    dropna=False,
):

    target_gene = TARGETS[target]["gene"]

    # --------------------------------------------------------
    # Publication evidence
    # --------------------------------------------------------

    pmids = set(
        group["pmid"]
        .dropna()
        .astype(str)
    )

    # --------------------------------------------------------
    # Method evidence
    # --------------------------------------------------------

    methods = set(
        group.loc[
            group["method"].notna(),
            "method"
        ].astype(str)
    )

    # --------------------------------------------------------
    # IntAct MI scores
    # --------------------------------------------------------

    intact_scores = pd.to_numeric(
        group["intact_miscore"],
        errors="coerce",
    ).dropna()

    # --------------------------------------------------------
    # STRING evidence
    # --------------------------------------------------------

    string_rows = group[
        group["source"] == "STRING"
    ]

    # --------------------------------------------------------
    # Source support
    # --------------------------------------------------------

    sources = set(group["source"])

    # --------------------------------------------------------
    # Partner gene
    # --------------------------------------------------------

    partner_genes = (
        group["partner_gene"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    partner_gene = (
        partner_genes.iloc[0]
        if len(partner_genes) > 0
        else pd.NA
    )

    rows.append(
        {
            "target_gene": target_gene,
            "target_uniprot": target,
            "partner_uniprot": partner,
            "partner_gene": partner_gene,

            "supported_by_intact":
                int("IntAct" in sources),

            "supported_by_biogrid":
                int("BioGRID" in sources),

            "supported_by_string":
                int("STRING" in sources),

            "n_sources":
                len(sources),

            "n_unique_pmids":
                len(pmids),

            "n_unique_methods":
                len(methods),

            "n_intact_records":
                int(
                    (group["source"] == "IntAct").sum()
                ),

            "n_biogrid_records":
                int(
                    (group["source"] == "BioGRID").sum()
                ),

            "max_intact_miscore":
                (
                    intact_scores.max()
                    if len(intact_scores)
                    else pd.NA
                ),

            "median_intact_miscore":
                (
                    intact_scores.median()
                    if len(intact_scores)
                    else pd.NA
                ),

            "string_experimental":
                (
                    string_rows["string_experimental"]
                    .max()
                    if len(string_rows)
                    else pd.NA
                ),

            "string_database":
                (
                    string_rows["string_database"]
                    .max()
                    if len(string_rows)
                    else pd.NA
                ),

            "string_textmining":
                (
                    string_rows["string_textmining"]
                    .max()
                    if len(string_rows)
                    else pd.NA
                ),

            "string_combined_score":
                (
                    string_rows["string_combined_score"]
                    .max()
                    if len(string_rows)
                    else pd.NA
                ),
        }
    )


partner_evidence = pd.DataFrame(rows)


# ============================================================
# 10. Save partner-level evidence
# ============================================================

partner_evidence = partner_evidence.sort_values(
    [
        "target_gene",
        "n_sources",
        "n_unique_pmids",
        "partner_gene",
    ],
    ascending=[
        True,
        False,
        False,
        True,
    ],
).reset_index(drop=True)

partner_file = (
    PROCESSED_DIR / "partner_evidence.csv"
)

partner_evidence.to_csv(
    partner_file,
    index=False,
)

print(f"Saved partner evidence to:")
print(partner_file)


# ============================================================
# 11. Build interaction-level network
# ============================================================

network = partner_evidence.copy()

network_file = (
    PROCESSED_DIR / "integrated_network.csv"
)

network.to_csv(
    network_file,
    index=False,
)

print(f"Saved integrated network to:")
print(network_file)


# ============================================================
# 12. Summary
# ============================================================

print()
print("=" * 70)
print("INTEGRATION COMPLETE")
print("=" * 70)

print(
    f"Evidence records: "
    f"{len(evidence):,}"
)

print(
    f"Unique target-partner interactions: "
    f"{len(network):,}"
)

for gene in ["CTBP1", "CTBP2"]:

    subset = network[
        network["target_gene"] == gene
    ]

    print()
    print(f"{gene}:")
    print(
        f"  Unique partners: "
        f"{len(subset):,}"
    )

    print(
        f"  IntAct-supported: "
        f"{subset['supported_by_intact'].sum():,}"
    )

    print(
        f"  BioGRID-supported: "
        f"{subset['supported_by_biogrid'].sum():,}"
    )

    print(
        f"  STRING-supported: "
        f"{subset['supported_by_string'].sum():,}"
    )
    all_three = (
        (subset["supported_by_intact"] == 1)
        & (subset["supported_by_biogrid"] == 1)
        & (subset["supported_by_string"] == 1)
    ).sum()
    print(f"  Supported by all 3: {all_three:,}")
    print(
        f"  >=2 unique PMIDs: "
        f"{(subset['n_unique_pmids'] >= 2).sum():,}"
    )

print()
print("Output files:")
print(evidence_file)
print(partner_file)
print(network_file)