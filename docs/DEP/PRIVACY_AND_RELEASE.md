# Privacy and release controls

1. Raw clinical, biomarker, genotype and SOMAscan files are configured locally and never embedded in code.
2. `data_private/`, `result/` and `publication_candidate/` are ignored by Git.
3. Participant-set exports are routed to `result/private/participant_sets/`.
4. PCA matrices and score tables containing direct identifiers are routed to `result/private/pca/`; a privacy-safe copy without stable IDs is generated for reporting.
5. The Source Data generator stops when a direct or stable identifier column is detected.
6. Publication workbooks are candidates only. They must be audited independently before submission or Zenodo deposition.
7. The code does not provide access to restricted ReDLat data and does not reconstruct participant identities.
