from pathlib import Path
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_DIR = PROJECT_ROOT / "data" / "raw" / "intact"
OUTPUT_DIR = PROJECT_ROOT / "data" / "interim"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# MITAB helper functions
# ============================================================

def extract_uniprot(value):
    """Extract the first UniProt accession from a MITAB field."""

    if pd.isna(value):
        return None

    for item in str(value).split("|"):
        item = item.strip()

        if item.startswith("uniprotkb:"):
            accession = item.split(":", 1)[1]

            # Remove isoform suffix if present
            accession = accession.split("-")[0]

            return accession

    return None


def extract_taxid(value):
    """Extract taxonomy ID from a MITAB field."""

    if pd.isna(value):
        return None

    value = str(value)

    if value.startswith("taxid:"):
        return value.split(":", 1)[1].split("(")[0]

    return None


def extract_pubmed(value):
    """Extract PubMed IDs from a MITAB publication field."""

    if pd.isna(value):
        return None

    pmids = []

    for item in str(value).split("|"):
        item = item.strip()

        if item.startswith("pubmed:"):
            pmids.append(item.split(":", 1)[1])

    return ";".join(pmids) if pmids else None


def extract_miscore(value):
    """Extract IntAct MI score."""

    if pd.isna(value):
        return None

    for item in str(value).split("|"):
        item = item.strip()

        if item.startswith("intact-miscore:"):
            try:
                return float(item.split(":", 1)[1])
            except ValueError:
                return None

    return None


# ============================================================
# MITAB parser
# ============================================================

def parse_mitab(filepath):
    """Read and clean an IntAct MITAB 2.8 file."""

    print(f"Reading: {filepath.name}")

    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        dtype=str,
        low_memory=False,
    )

    print(f"  Raw records: {len(df)}")

    # MITAB 2.8 column names
    column_names = [
        "idA",
        "idB",
        "altA",
        "altB",
        "aliasA",
        "aliasB",
        "detection_method",
        "first_author",
        "publication",
        "taxidA",
        "taxidB",
        "interaction_type",
        "source_database",
        "interaction_id",
        "confidence",
        "complex_expansion",
        "biological_role_A",
        "biological_role_B",
        "experimental_role_A",
        "experimental_role_B",
        "type_A",
        "type_B",
        "intact_participant_A",
        "intact_participant_B",
        "source",
        "parameter",
        "negative",
        "inferred_by",
        "conflict",
        "pathway",
        "participants",
        "stoichiometry",
        "experimental_data",
        "host_organism",
        "host_taxid",
        "interaction_detection_method",
        "interaction_parameter",
        "curation",
        "cross_references",
        "annotation_A",
        "annotation_B",
        "complex",
    ]

    # Use only the number of columns actually present
    df = df.iloc[:, :len(column_names)]
    df.columns = column_names[:df.shape[1]]

    # Extract useful standardized fields
    df["proteinA"] = df["idA"].apply(extract_uniprot)
    df["proteinB"] = df["idB"].apply(extract_uniprot)

    df["taxidA_clean"] = df["taxidA"].apply(extract_taxid)
    df["taxidB_clean"] = df["taxidB"].apply(extract_taxid)

    df["pubmed"] = df["publication"].apply(extract_pubmed)

    df["intact_miscore"] = df["confidence"].apply(extract_miscore)

    # Keep only human-human interactions
    df = df[
        (df["taxidA_clean"] == "9606") &
        (df["taxidB_clean"] == "9606")
    ].copy()

    print(f"  Human-human records: {len(df)}")

    # Keep only interactions where both proteins have UniProt IDs
    df = df[
        df["proteinA"].notna() &
        df["proteinB"].notna()
    ].copy()

    print(f"  Records with UniProt IDs: {len(df)}")

    # Select final columns
    df = df[
        [
            "proteinA",
            "proteinB",
            "detection_method",
            "first_author",
            "pubmed",
            "interaction_type",
            "source_database",
            "interaction_id",
            "intact_miscore",
        ]
    ]

    return df


# ============================================================
# Main
# ============================================================

def main():

    files = [
        INPUT_DIR / "Q13363.mitab",
        INPUT_DIR / "P56545.mitab",
    ]

    all_data = []

    for filepath in files:

        if not filepath.exists():
            raise FileNotFoundError(
                f"File not found: {filepath}"
            )

        parsed = parse_mitab(filepath)
        all_data.append(parsed)

    # Combine CtBP1 + CtBP2
    combined = pd.concat(
        all_data,
        ignore_index=True
    )

    # Save
    output_file = OUTPUT_DIR / "intact_raw.csv"

    combined.to_csv(
        output_file,
        index=False
    )

    print()
    print("=" * 60)
    print("DONE")
    print("=" * 60)
    print(f"Total records: {len(combined)}")
    print(f"Saved to: {output_file}")


if __name__ == "__main__":
    main()