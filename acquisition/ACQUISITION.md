# Acquire the original sources independently

No data downloads or raw-source archives are supplied here. Consult source-specific terms before obtaining or using the data. The source manifest is an identity/field map, not a license or a frozen execution authorization.

| Source | DOI / official catalog | Locally required files |
|---|---|---|
| NLR simulated wind | https://doi.org/10.7799/3007835 ; https://data.nlr.gov/submissions/305 | Accepted wind archives in `nlr_archives`; characterization archive in `characterization` (also listed by wave catalog) |
| NLR simulated wave | https://doi.org/10.7799/3007836 ; https://data.nlr.gov/submissions/306 | Accepted wave archives in `nlr_archives` |
| NLR historical wind | https://doi.org/10.7799/3024096 ; https://data.nlr.gov/submissions/316 | Combined historical wind CSV |
| NLR simulated tidal | https://doi.org/10.7799/3364185 ; https://data.nlr.gov/submissions/317 | Combined simulated tidal CSV |
| NLR historical PV | https://doi.org/10.7799/3372313 ; https://data.nlr.gov/submissions/318 | Combined solar PV CSV |
| Jülich 100 kW V1.0 | https://doi.org/10.26165/JUELICH-DATA/FEKKSS | 112 development and 22 evaluation measurement workbooks listed by public file identifier |

NLR dataset citation: Abel, Riley; Schwarz, Marty; Leighton, Daniel; Gevorgian, Vahan; Wimer, Nicholas (2025, source-specific records), *Public Reference Data for Megawatt-Scale Hydrogen Electrolysis*. National Laboratory of the Rockies Data Catalog. Cite each appropriate DOI and use the current catalog citation. The characterization archive's original local name is used in the source map; providers may present the same resource without the numeric download prefix. Rename the local file to the manifest filename only if its bytes match. Do not regenerate the combined CSVs or change timestamp formats to satisfy the checks. Some provider files may have been updated or replaced; stop if exact historical bytes cannot be obtained. The source-specific production license pages could not be retrieved at the 2026-11-25 review. No NLR data, derived numerical results or licensing assertion is included.

Jülich citation: Keller, Roger; Rauls, Edward; Hehemann, Michael; Müller, Martin (2025). *Operational Data of a 100 kW Proton Exchange Membrane (PEM) Electrolyzer: Collection 2018–2023*, V1.0. Jülich DATA. https://doi.org/10.26165/JUELICH-DATA/FEKKSS. Upstream CC0 1.0 as indicated by the official indexed dataset record. Download the exact measurement resource identified by `file_id`, directory and sheet in `SOURCE_SPEC.json`; save the copy using `<file_id>.xlsx` as specified (the archival local names). In the Dataverse interface, check file IDs and SHA256; the DOI landing page and its resource links are the authoritative routes. Do not download the paired setpoint workbook for this workflow. Keep independent source directories for development and evaluation. The 112/22 counts are original candidate measurement-file counts, not counts of eligible replay files. No author project path or credential is required.

Hash verification does not certify instrument calibration, semantic equivalence of a newer provider revision, or redistribution rights. Cite the sources and respect their license conditions.
