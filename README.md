# CtBP1/CtBP2 Protein-Protein Interaction Prediction

A computational pipeline for constructing, analyzing, and predicting the protein-protein interaction (PPI) network of **C-terminal-binding protein 1 (CtBP1)** and **C-terminal-binding protein 2 (CtBP2)**.

The project integrates experimentally supported protein interaction data from multiple biological databases and subsequently applies computational prediction and structural validation methods to identify potential novel CtBP interactors.

---

## Project Objectives

The known interaction network is not the output of this project. It is the
mechanistic reference from which CtBP-specific features are derived, and the
benchmark against which predictive methods are calibrated. The output is a
ranked shortlist of **unreported** candidate partners worth testing
experimentally. The claim made is prioritisation rather than accuracy: that
testing the shortlist is a better use of experimental effort than testing an
equivalent number of randomly chosen nuclear proteins.

### Phase 1 — Integrated reference network *(complete)*

1. Construct a high-confidence protein–protein interaction network for CtBP1 and CtBP2.
2. Integrate interaction evidence for CtBP1 and CtBP2 from multiple public protein–protein interaction databases, including STRING, BioGRID, and IntAct.
3. Standardize protein identifiers and interaction attributes across the different databases to enable consistent integration.
4. Identify and consolidate duplicate interaction records while preserving their experimental and database provenance.
5. Develop a composite confidence-scoring framework to rank known CtBP1/CtBP2 interactions based on the strength and consistency of supporting evidence.

### Phase 2 — Benchmark design *(prerequisite for all prediction work)*

6. Resolve every PMID in the evidence records to a publication year and define a
   temporal holdout: CtBP partners first reported after a fixed cutoff form a
   blind test set, as they cannot have been present in the training data of the
   predictive methods used.
7. Construct a prior-matched negative set, sampled at the realistic
   positive:negative ratio of a proteome-wide screen rather than balanced 1:1.
8. Evaluate all predictive methods by precision–recall (PR-AUC, precision@k) and
   by enrichment over a random baseline. Accuracy and AUC-ROC are not reported,
   as both are uninformative at this class imbalance.

### Phase 3 — Prediction

9. Predict potential novel CtBP1/CtBP2 interacting proteins using sequence-based
   and structure-informed PPI prediction approaches, including SPRINT and PrePPI.
10. Develop a CtBP-specific motif predictor: scan the reviewed human proteome for
    PxDLS/PLDLS variants, filtered for nuclear localisation and for occurrence in
    accessible or disordered sequence context.
11. Define an explicit rank-aggregation rule combining sequence, structural,
    motif, and localisation evidence into a single consensus ranking, requiring
    support from at least two orthogonal lines of evidence.

### Phase 4 — Validation and interpretation

12. Validate top-ranked novel predictions using structural modelling with
    AlphaFold-Multimer, with the confidence threshold calibrated on known CtBP
    partners and sequence-shuffled decoys rather than taken as a default.
13. Apply a degree-matched control to confirm that predicted partners are not
    simply high-degree hub proteins recovered by any generic PPI predictor.
14. Perform functional enrichment and pathway-level analysis of predicted CtBP1/CtBP2 interactors to characterize their potential biological roles.
15. Visualize and compare the integrated CtBP1 and CtBP2 interaction networks to identify shared, unique, and highly supported interaction partners.

---

## Methodological Constraints

Four constraints govern how predictions are evaluated, and determine the design
of Phase 2.

**Training-data leakage.** SPRINT is trained on human PPI reference sets derived
from BioGRID/IntAct, and PrePPI's Bayesian framework incorporates existing
interaction-database evidence. Known CtBP interactions are therefore likely to be
present in the training data of the methods being evaluated, so recovering them
measures memorisation rather than prediction. Addressed by the temporal holdout
(objective 6). Where a method is retrainable, CtBP pairs are additionally
excluded from training; PrePPI is distributed as a precomputed resource and
cannot be retrained, so the temporal split is the only available control and this
limitation is reported explicitly.

**Recall is not precision.** Recovery of known partners measures sensitivity. The
quantity of interest is precision on novel calls, which degrades sharply under
proteome-scale class imbalance: at roughly one true partner per seventy candidate
proteins, 95% specificity still yields far more false positives than true ones.

**Absence of evidence is not evidence of absence.** A protein not recorded as a
CtBP partner is unstudied, not demonstrated to be a non-interactor. Measured
precision is therefore a lower bound, and a prediction is not counted as
incorrect solely because it is absent from the databases.

**Study and degree bias.** Heavily studied and highly connected proteins
accumulate more recorded interactions, and generic predictors can learn
promiscuity rather than CtBP-specific binding. Addressed by objective 13.

The PxDLS motif enrichment established in Phase 1 (PLDLS, odds ratio 86.7,
p = 4.6e-17) is mechanistically specific to CtBP and independent of the PPI
databases, making it the one evidence line carrying no leakage risk. Its recall
is limited — 78% of high-confidence partners carry no PxDLS match — so it serves
as a high-precision confirming filter rather than a primary screen.


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