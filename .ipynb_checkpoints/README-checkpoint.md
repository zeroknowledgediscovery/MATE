# MATE: Major Adverse Transplant Episode

This repository contains working definitions and analysis materials for evaluating the **Major Adverse Transplant Episode (MATE)** as a patient-focused outcome after lung transplantation, with current emphasis on operationalizing the phenotype in Epic COSMOS.

## Original MATE definition

The proposed clinical MATE endpoint defines a poor first-year post-lung-transplant outcome as the occurrence of **any** of the following:

1. **Graft failure within 1 year**, defined as death or retransplantation.
2. **Index transplant hospitalization >90 days**.
3. **Dialysis after day 90** post-transplant.
4. **>3 readmissions during the first year** after transplantation (operationally, >=4 readmissions).

The intent is to distinguish patients who are merely alive at 1 year from patients who are **alive and well**, using severe outcomes that are clinically meaningful, scalable, and suitable for validation against patient-reported quality-of-life measures.

## Current Epic COSMOS analysis

Epic COSMOS is being used to establish scalable EHR ascertainment and to determine how the non-death components of MATE behave in a national lung-transplant cohort. **Death is not included in the current COSMOS working phenotype because death ascertainment is less reliable in COSMOS than in transplant registry/LTC data.** Death will remain part of the definitive clinical MATE endpoint and will be evaluated using data sources with more reliable mortality ascertainment.

The current COSMOS analysis should therefore be treated as **MATE-EHR / MATE-proximal**, not as the final validated MATE definition.

### Components to recompute in COSMOS

- Retransplantation within 365 days of index lung transplant.
- Index transplant hospitalization length of stay >90 days.
- Dialysis occurring after day 90.
- Unplanned inpatient readmissions after discharge from the index transplant hospitalization.
- Oxygen dependence after day 90 as a **sensitivity/QOL candidate only**, not a primary MATE component.

### Candidate non-death MATE-EHR definitions

**MATE-EHR v1**  
Retransplantation OR index LOS >90 days OR >=1 dialysis event on days 91-365 OR >=4 unplanned readmissions.

**MATE-EHR v2**  
Retransplantation OR index LOS >90 days OR persistent dialysis after day 90 (>=2 dialysis events on distinct dates, preferably >=7 days apart) OR >=4 unplanned readmissions.

**MATE-EHR v3**  
Retransplantation OR index LOS >90 days OR persistent dialysis after day 90 OR >=5 unplanned readmissions.

The preferred primary definition should be selected based on clinical validity and component behavior, not by forcing a particular prevalence. As a working plausibility check, we expect the **non-death COSMOS MATE-EHR phenotype to fall roughly in the 25-30% range**, which would be compatible with a final clinical MATE prevalence around **30-35% after death is incorporated and overlap among components is accounted for**. These are evaluation ranges, not prevalence-calibration targets.

### Secondary burden phenotype

A secondary **high-burden MATE** phenotype can be defined as >=2 non-death MATE components. This should be used as a severity/burden analysis rather than replacing the primary MATE definition, because a single event such as retransplantation or an index hospitalization >90 days is independently sufficient to represent a major adverse transplant outcome.

## Files

- `COSMOS_MATE_Definition_Worksheet.docx` - shareable worksheet for entering COSMOS counts and prevalence estimates.
- `COSMOS_MATE_Definition_Worksheet.md` - Markdown version of the same computation plan, including overlap/burden checks.

## Immediate analysis goal

Fill the worksheet from COSMOS, compare MATE-EHR v1-v3, inspect component overlap, and determine whether the clinically natural original threshold of **>3 (>=4) readmissions** is appropriate or whether the EHR implementation requires a more specific dialysis definition and/or a higher readmission threshold.
