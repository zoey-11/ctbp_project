# CtBP1/CtBP2 Protein-Protein Interaction Prediction

A computational pipeline for constructing, analyzing, and predicting the protein-protein interaction (PPI) network of **C-terminal-binding protein 1 (CtBP1)** and **C-terminal-binding protein 2 (CtBP2)**.

The project integrates experimentally supported protein interaction data from multiple biological databases and subsequently applies computational prediction and structural validation methods to identify potential novel CtBP interactors.

---

## Project Objectives

The primary objectives of this project are:

1. Construct a high-confidence protein–protein interaction network for CtBP1 and CtBP2.
2. Integrate interaction evidence for CtBP1 and CtBP2 from multiple public protein–protein interaction databases, including STRING, BioGRID, and IntAct.
3. Standardize protein identifiers and interaction attributes across the different databases to enable consistent integration.
4. Identify and consolidate duplicate interaction records while preserving their experimental and database provenance.
5. Develop a composite confidence-scoring framework to rank known CtBP1/CtBP2 interactions based on the strength and consistency of supporting evidence.
6. Predict potential novel CtBP1/CtBP2 interacting proteins using sequence-based and structure-informed PPI prediction approaches, including SPIRNT and PrePPI.
7. Validate high-confidence novel interaction predictions using structural modeling with AlphaFold-Multimer and supporting biological evidence.
8. Perform functional enrichment and pathway-level analysis of predicted CtBP1/CtBP2 interactors to characterize their potential biological roles.
9. Visualize and compare the integrated CtBP1 and CtBP2 interaction networks to identify shared, unique, and highly supported interaction partners.


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