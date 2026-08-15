# CtBP1/CtBP2 Protein-Protein Interaction Prediction

A computational pipeline for constructing, analyzing, and predicting the protein-protein interaction (PPI) network of **C-terminal-binding protein 1 (CtBP1)** and **C-terminal-binding protein 2 (CtBP2)**.

The project integrates experimentally supported protein interaction data from multiple biological databases and subsequently applies computational prediction and structural validation methods to identify potential novel CtBP interactors.

---

## Project Objectives

The primary objectives of this project are:

1. Construct a high-confidence CtBP1/CtBP2 protein-protein interaction network.
2. Integrate interaction evidence from multiple public databases.
3. Standardize protein identifiers across databases.
4. Remove duplicate interaction records while preserving experimental provenance.
5. Rank known interactions based on the strength of supporting evidence.
6. Predict potential novel CtBP1/CtBP2 interacting proteins.
7. Validate high-confidence predictions using structural and biological evidence.
8. Perform functional and pathway-level analysis of predicted interactors.

---

## Target Proteins

| Protein | Gene | UniProt Accession | Organism |
|---|---|---|---|
| CtBP1 | CTBP1 | Q13363 | Homo sapiens |
| CtBP2 | CTBP2 | P56545 | Homo sapiens |

---

# Project Structure

```text
ctbp_project/
│
├── data/
│   ├── raw/
│   │   ├── uniprot/
│   │   ├── biogrid/
│   │   ├── string/
│   │   └── intact/
│   │
│   ├── interim/
│   └── processed/
│
├── scripts/
│
├── results/
│
├── figures/
│
├── logs/
│
├── environment.yml
├── README.md
└── .gitignore