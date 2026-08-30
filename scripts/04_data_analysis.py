"""Exploratory analysis of the parsed CtBP1/CtBP2 interaction sources.

Runs before integration. Characterises each source on its own terms,
measures how far the three sources agree, and tests whether the
partners they agree on are enriched for the PxDLS motif through which
CtBP is known to recruit its binding partners.

Nothing here writes to data/processed/. The purpose is to establish
what the integration script has to reconcile, and to provide a
biological sanity check on the network before it is built.
"""

from pathlib import Path
import re

import pandas as pd
from scipy.stats import fisher_exact


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INTERIM_DIR = PROJECT_ROOT / "data" / "interim"
UNIPROT_DIR = PROJECT_ROOT / "data" / "raw" / "uniprot"
RESULTS_DIR = PROJECT_ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

INTACT_FILE = INTERIM_DIR / "intact_raw.csv"
BIOGRID_FILE = INTERIM_DIR / "biogrid_raw.csv"
STRING_FILE = INTERIM_DIR / "string_raw.csv"

UNIPROT_REFERENCE = UNIPROT_DIR / "uniprot_human_reviewed.tsv.gz"


# ============================================================
# Target proteins
# ============================================================

TARGETS = {
    "Q13363": "CTBP1",
    "P56545": "CTBP2",
}

TARGET_ACCESSIONS = set(TARGETS.keys())
TARGET_SYMBOLS = set(TARGETS.values())


# ============================================================
# CtBP recruitment motif
# ============================================================

# CtBP recruits transcriptional regulators through a short linear
# motif of consensus PxDLS. PLDLS is the canonical form; the
# extended PxDLSxK variant adds a downstream lysine.
MOTIFS = {
    "PxDLS": r"P.DL[SV]",
    "PLDLS": r"PLDLS",
    "PxDLSxK": r"P.DL[SV].K",
}


# ============================================================
# UniProt reference proteome
# ============================================================

print("Reading UniProt reviewed human proteome...")

reference = pd.read_csv(
    UNIPROT_REFERENCE,
    sep="\t",
    usecols=[
        "Entry",
        "Gene Names",
        "Sequence",
        "Subcellular location [CC]",
    ],
    dtype=str,
)

reference = reference.dropna(subset=["Sequence"])

# The first token of "Gene Names" is the primary HGNC symbol.
reference["symbol"] = (
    reference["Gene Names"].fillna("").str.split().str[0]
)

reference = reference[reference["symbol"] != ""]

reference["nuclear"] = (
    reference["Subcellular location [CC]"]
    .fillna("")
    .str.contains("Nucleus")
)

# Built from every reviewed accession, before symbols are
# collapsed, so that accession lookup stays complete.
accession_to_symbol = dict(
    zip(reference["Entry"], reference["symbol"])
)

known_symbols = set(reference["symbol"])

# The motif background needs one entry per gene rather than per
# accession, so a symbol shared by several reviewed entries is
# collapsed to a single deterministic representative.
symbol_collisions = int(reference["symbol"].duplicated().sum())

proteome = (
    reference
    .sort_values("Entry")
    .drop_duplicates("symbol")
    .reset_index(drop=True)
)

print(
    f"Reviewed human accessions: {len(reference):,} | "
    f"distinct genes: {len(proteome):,} "
    f"({symbol_collisions:,} duplicate symbols collapsed)"
)


# ============================================================
# Partner extraction
# ============================================================

def summarise(name, records, partners, unmapped, notes):
    """Package one source's inventory into a single row."""

    return {
        "source": name,
        "records": records,
        "unique_partners": len(partners),
        "unmapped_identifiers": unmapped,
        "notes": notes,
    }


# ------------------------------------------------------------
# IntAct
# ------------------------------------------------------------

print("\nReading IntAct...")

intact = pd.read_csv(INTACT_FILE)

intact_records = len(intact)

# Interaction identifiers repeat when a curated complex is
# expanded into pairwise records, so distinct publications are a
# better measure of independent evidence than row count.
intact_publications = intact["pubmed"].nunique()

intact_largest_study = (
    intact["pubmed"].value_counts().iloc[0]
    if intact_publications
    else 0
)

intact = intact[intact["proteinA"] != intact["proteinB"]]

intact_partner_accessions = (
    intact["proteinB"]
    .where(
        intact["proteinA"].isin(TARGET_ACCESSIONS),
        intact["proteinA"],
    )
)

intact_partner_accessions = set(
    intact_partner_accessions.dropna()
) - TARGET_ACCESSIONS

# Accessions absent from the reviewed proteome are unreviewed or
# obsolete. They are counted rather than silently dropped, since
# they represent real integration loss.
intact_unmapped = {
    accession
    for accession in intact_partner_accessions
    if accession not in accession_to_symbol
}

intact_partners = {
    accession_to_symbol[accession]
    for accession in intact_partner_accessions
    if accession in accession_to_symbol
} - TARGET_SYMBOLS

print(
    f"  {intact_records:,} records | "
    f"{intact_publications:,} publications | "
    f"{len(intact_partners):,} partners"
)


# ------------------------------------------------------------
# BioGRID
# ------------------------------------------------------------

print("\nReading BioGRID...")

biogrid = pd.read_csv(BIOGRID_FILE, low_memory=False)

biogrid_records = len(biogrid)

biogrid_genetic = (
    biogrid["Experimental System Type"] == "genetic"
).sum()

# Genetic interactions are not physical binding evidence.
biogrid = biogrid[
    biogrid["Experimental System Type"] == "physical"
]

biogrid_high_throughput = (
    biogrid["Throughput"]
    .fillna("")
    .str.contains("High Throughput")
    .sum()
)

symbol_a = biogrid["Official Symbol Interactor A"]
symbol_b = biogrid["Official Symbol Interactor B"]

biogrid_partner_symbols = symbol_b.where(
    symbol_a.isin(TARGET_SYMBOLS),
    symbol_a,
)

biogrid_partner_symbols = set(
    biogrid_partner_symbols.dropna()
) - TARGET_SYMBOLS

biogrid_unmapped = {
    symbol
    for symbol in biogrid_partner_symbols
    if symbol not in known_symbols
}

biogrid_partners = biogrid_partner_symbols

print(
    f"  {biogrid_records:,} records "
    f"({biogrid_genetic:,} genetic removed) | "
    f"{len(biogrid_partners):,} partners"
)


# ------------------------------------------------------------
# STRING
# ------------------------------------------------------------

print("\nReading STRING...")

string = pd.read_csv(STRING_FILE)

if "partner_gene" not in string.columns:
    raise RuntimeError(
        "string_raw.csv lacks a partner_gene column. "
        "Re-run scripts/03_parse_string.py."
    )

string_records = len(string)

string_direct = string[
    (string["experimental"] > 0) | (string["database"] > 0)
]

string_textmining_only = string_records - len(string_direct)

string_partners = set(
    string["partner_gene"].dropna()
) - TARGET_SYMBOLS

string_direct_partners = set(
    string_direct["partner_gene"].dropna()
) - TARGET_SYMBOLS

string_unmapped = {
    symbol
    for symbol in string_partners
    if symbol not in known_symbols
}

print(
    f"  {string_records:,} pairs "
    f"({string_textmining_only:,} textmining-only) | "
    f"{len(string_partners):,} partners"
)


# ============================================================
# Source inventory
# ============================================================

inventory = pd.DataFrame(
    [
        summarise(
            "IntAct",
            intact_records,
            intact_partners,
            len(intact_unmapped),
            f"{intact_publications} publications; "
            f"largest single study contributes "
            f"{intact_largest_study} records",
        ),
        summarise(
            "BioGRID",
            biogrid_records,
            biogrid_partners,
            len(biogrid_unmapped),
            f"{biogrid_genetic} genetic records excluded; "
            f"{biogrid_high_throughput} high-throughput physical "
            f"records retained",
        ),
        summarise(
            "STRING",
            string_records,
            string_partners,
            len(string_unmapped),
            f"{string_textmining_only} textmining-only pairs; "
            f"{len(string_direct_partners)} partners have "
            f"experimental or database support",
        ),
    ]
)

inventory.to_csv(
    RESULTS_DIR / "source_inventory.csv",
    index=False,
)


# ============================================================
# Cross-source agreement
# ============================================================

all_partners = sorted(
    intact_partners | biogrid_partners | string_partners
)

support = pd.DataFrame(
    {
        "partner_gene": all_partners,
    }
)

support["in_intact"] = support["partner_gene"].isin(intact_partners)
support["in_biogrid"] = support["partner_gene"].isin(biogrid_partners)
support["in_string"] = support["partner_gene"].isin(string_partners)

support["in_string_direct"] = support["partner_gene"].isin(
    string_direct_partners
)

support["n_sources"] = (
    support[["in_intact", "in_biogrid", "in_string"]].sum(axis=1)
)

support = support.sort_values(
    ["n_sources", "partner_gene"],
    ascending=[False, True],
).reset_index(drop=True)

support.to_csv(
    RESULTS_DIR / "partner_support_matrix.csv",
    index=False,
)

high_confidence = support[support["n_sources"] == 3]

high_confidence.to_csv(
    RESULTS_DIR / "high_confidence_partners.csv",
    index=False,
)

high_confidence_symbols = set(high_confidence["partner_gene"])


# ============================================================
# PxDLS motif enrichment
# ============================================================

# Enrichment is tested in the partners all three sources agree on,
# against every other reviewed human protein. A real CtBP network
# should be enriched; a network dominated by screen artefacts
# should not.
in_high_confidence = proteome["symbol"].isin(
    high_confidence_symbols
)

n_high_confidence = int(in_high_confidence.sum())
n_background = int((~in_high_confidence).sum())

enrichment_rows = []

for motif_name, pattern in MOTIFS.items():

    has_motif = proteome["Sequence"].str.contains(pattern)

    hits = int((has_motif & in_high_confidence).sum())
    background_hits = int((has_motif & ~in_high_confidence).sum())

    odds_ratio, p_value = fisher_exact(
        [
            [hits, n_high_confidence - hits],
            [background_hits, n_background - background_hits],
        ],
        alternative="greater",
    )

    enrichment_rows.append(
        {
            "motif": motif_name,
            "pattern": pattern,
            "partner_hits": hits,
            "partner_total": n_high_confidence,
            "partner_fraction": hits / n_high_confidence,
            "background_hits": background_hits,
            "background_total": n_background,
            "background_fraction": background_hits / n_background,
            "odds_ratio": odds_ratio,
            "p_value": p_value,
        }
    )

enrichment = pd.DataFrame(enrichment_rows)

enrichment.to_csv(
    RESULTS_DIR / "motif_enrichment.csv",
    index=False,
)

# Which high-confidence partners carry the canonical motif.
carries_pxdls = proteome["Sequence"].str.contains(MOTIFS["PxDLS"])

motif_partners = sorted(
    proteome.loc[carries_pxdls & in_high_confidence, "symbol"]
)


# ============================================================
# Subcellular localisation
# ============================================================

nuclear_high_confidence = (
    proteome.loc[in_high_confidence, "nuclear"].mean()
)

nuclear_background = (
    proteome.loc[~in_high_confidence, "nuclear"].mean()
)


# ============================================================
# Summary
# ============================================================

print()
print("=" * 60)
print("SOURCE INVENTORY")
print("=" * 60)

for row in inventory.itertuples():
    print(f"\n{row.source}")
    print(f"  Records: {row.records:,}")
    print(f"  Unique partners: {row.unique_partners:,}")
    print(f"  Unmapped identifiers: {row.unmapped_identifiers:,}")
    print(f"  {row.notes}")


print()
print("=" * 60)
print("CROSS-SOURCE AGREEMENT")
print("=" * 60)

print(f"\nDistinct partners across all sources: {len(support):,}")

print()

for n in [3, 2, 1]:
    count = (support["n_sources"] == n).sum()
    print(f"  Supported by {n} source(s): {count:,}")

print("\nPairwise overlap:")
print(
    f"  IntAct  & BioGRID: "
    f"{len(intact_partners & biogrid_partners):,}"
)
print(
    f"  IntAct  & STRING : "
    f"{len(intact_partners & string_partners):,}"
)
print(
    f"  BioGRID & STRING : "
    f"{len(biogrid_partners & string_partners):,}"
)

print("\nSource-exclusive partners:")
print(
    f"  IntAct  only: "
    f"{len(intact_partners - biogrid_partners - string_partners):,}"
)
print(
    f"  BioGRID only: "
    f"{len(biogrid_partners - intact_partners - string_partners):,}"
)
print(
    f"  STRING  only: "
    f"{len(string_partners - intact_partners - biogrid_partners):,}"
)


print()
print("=" * 60)
print("PxDLS MOTIF ENRICHMENT")
print("=" * 60)

print(
    f"\nHigh-confidence partners resolved in UniProt: "
    f"{n_high_confidence:,} of {len(high_confidence_symbols):,}"
)

print()
print(
    f"{'motif':10s} {'partners':>16s} {'background':>18s} "
    f"{'OR':>8s} {'p':>11s}"
)

for row in enrichment.itertuples():
    print(
        f"{row.motif:10s} "
        f"{row.partner_hits:5d}/{row.partner_total:<5d} "
        f"({100 * row.partner_fraction:4.1f}%) "
        f"{row.background_hits:6d}/{row.background_total:<6d} "
        f"({100 * row.background_fraction:4.1f}%) "
        f"{row.odds_ratio:8.2f} {row.p_value:11.2e}"
    )

print(f"\nHigh-confidence partners carrying PxDLS ({len(motif_partners)}):")
print("  " + ", ".join(motif_partners))

print()
print("Nuclear localisation:")
print(f"  High-confidence partners: {100 * nuclear_high_confidence:.1f}%")
print(f"  Background proteome:      {100 * nuclear_background:.1f}%")


print()
print("=" * 60)
print("OUTPUTS")
print("=" * 60)

for filename in [
    "source_inventory.csv",
    "partner_support_matrix.csv",
    "high_confidence_partners.csv",
    "motif_enrichment.csv",
]:
    print(f"  {RESULTS_DIR / filename}")
