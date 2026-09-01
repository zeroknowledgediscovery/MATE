# COSMOS MATE Definition Worksheet

Fill `N` and `Prevalence (%)` after rerunning the clinically aligned outcome definitions in Epic COSMOS.

**Working plausibility check:** approximately 25-30% non-death MATE-EHR prevalence in COSMOS; the final clinical MATE is expected to be approximately 30-35% after death is incorporated in LTC/STAR and overlap is accounted for. These ranges are for plausibility assessment, not prevalence calibration.

| Candidate component / definition | Exact COSMOS operational definition | N | Prevalence (%) | Decision / notes |
|---|---|---:|---:|---|
| Retransplantation | Retransplantation within 365 days of index lung transplant. |  |  | Core MATE component; retain. |
| Index hospitalization >90 d | Index transplant admission LOS >90 days. |  |  | Core MATE component; add to COSMOS proxy. |
| Dialysis after day 90 - any | >=1 dialysis/hemodialysis event on days 91-365. |  |  | Closest implementation of proposed MATE. |
| Dialysis after day 90 - persistent | >=2 dialysis events on distinct dates during days 91-365, preferably >=7 days apart. |  |  | Specificity check for persistent dialysis. |
| >=3 unplanned readmissions | >=3 unplanned inpatient readmissions after index discharge through day 365. |  |  | Current preliminary threshold; likely permissive. |
| >=4 unplanned readmissions | >=4 unplanned inpatient readmissions after index discharge through day 365. |  |  | Matches proposal's >3 readmissions; test first. |
| >=5 unplanned readmissions | >=5 unplanned inpatient readmissions after index discharge through day 365. |  |  | Main tightening option if >=4 is too prevalent. |
| >=6 unplanned readmissions | >=6 unplanned inpatient readmissions after index discharge through day 365. |  |  | Sensitivity analysis. |
| Oxygen dependence after day 90 | Oxygen dependence / oxygen DME evidence on or after day 90. |  |  | Sensitivity/QOL candidate; exclude from primary MATE. |
| MATE-EHR v1 | Retransplant OR index LOS >90 d OR any dialysis after day 90 OR >=4 unplanned readmissions. |  |  | Primary definition to evaluate first. |
| MATE-EHR v2 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=4 unplanned readmissions. |  |  | Cleaner EHR-specific candidate. |
| MATE-EHR v3 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=5 unplanned readmissions. |  |  | Use if v1/v2 remain materially >30%. |
| High-burden MATE | >=2 non-death MATE components. |  |  | Secondary severity construct, not primary MATE. |

## Overlap / burden checks

| Overlap / burden calculation | N | % of cohort |
|---|---:|---:|
| Exactly 0 MATE-EHR components |  |  |
| Exactly 1 component |  |  |
| Exactly 2 components |  |  |
| Exactly 3 components |  |  |
| All 4 components |  |  |
| Retransplantation + any other component |  |  |
| Index LOS >90 d + dialysis after day 90 |  |  |
| Index LOS >90 d + high readmission burden |  |  |
| Dialysis after day 90 + high readmission burden |  |  |

## Decision rule

If MATE-EHR v1 is approximately 25-30%, retain the clinically natural >=4 readmission threshold. If it remains materially above 30%, evaluate v2 and then v3. Do not require >=2 components for primary MATE; use that only as a secondary high-burden phenotype.
