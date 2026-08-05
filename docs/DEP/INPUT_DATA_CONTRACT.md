# DEP input data contract

## Clinical metadata

The metadata CSV must contain, at minimum:

- `SampleId`: local study identifier used only for joining private inputs;
- `SampleGroup`: `CN` or `AD` for the primary analysis;
- `Sex`;
- `Age`;
- `Country`;
- `Education`;
- `ApoE` or the APOE-derived variable expected by the original pipeline.

Additional variables are detected for sensitivity analyses, including CDR-SB, p-tau217, Aβ42/40, NfL and other clinical traits. Their exact accepted aliases remain documented in the relevant canonical script.

## SOMAscan ADAT

The ADAT file must contain `SampleId`, assay values and the internal analyte annotation required to derive `AptName`, gene symbols, Entrez IDs, UniProt identifiers, organism and target names.

## Governance

Inputs must be stored in a restricted local location. The pipeline does not copy raw inputs into the Git repository or publication-candidate directory. Any transfer, access or reuse remains governed by ReDLat approvals and data-use agreements.
