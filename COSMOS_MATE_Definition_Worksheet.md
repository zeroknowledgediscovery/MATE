# COSMOS MATE Definition Worksheet

Fill `N` and `Prevalence (%)` after rerunning the clinically aligned outcome definitions in Epic COSMOS.

**Working plausibility check:** approximately 25-30% non-death MATE-EHR prevalence in COSMOS; the final clinical MATE is expected to be approximately 30-35% after death is incorporated in LTC/STAR and overlap is accounted for. These ranges are for plausibility assessment, not prevalence calibration.

| Candidate component / definition | Exact COSMOS operational definition | N | Prevalence (%) | Decision / notes |
|---|---|---:|---:|---|
| Retransplantation | Retransplantation within 365 days of index lung transplant. | 730 | 2.658 | Core MATE component; retain. |
| Index hospitalization >90 d | Index transplant admission LOS >90 days. | 525 | 1.911 | Core MATE component; add to COSMOS proxy. |
| Dialysis after day 90 - any | >=1 dialysis/hemodialysis event on days 91-365. | 4,113 | 14.973 | Closest implementation of proposed MATE. |
| Dialysis after day 90 - persistent | >=2 dialysis events on distinct dates during days 91-365, preferably >=7 days apart. | 1,066 | 3.881 | Specificity check for persistent dialysis. |
| >=3 unplanned readmissions | >=3 unplanned inpatient readmissions after index discharge through day 365. | 6,570 | 23.918 | Current preliminary threshold; likely permissive. |
| >=4 unplanned readmissions | >=4 unplanned inpatient readmissions after index discharge through day 365. | 4,032 | 14.678 | Matches proposal's >3 readmissions; test first. |
| >=5 unplanned readmissions | >=5 unplanned inpatient readmissions after index discharge through day 365. | 2,392 | 8.708 | Main tightening option if >=4 is too prevalent. |
| >=6 unplanned readmissions | >=6 unplanned inpatient readmissions after index discharge through day 365. | 1,446 | 5.264 | Sensitivity analysis. |
| Oxygen dependence after day 90 | Oxygen dependence / oxygen DME evidence on or after day 90. | 3,833 | 13.954 | Sensitivity/QOL candidate; exclude from primary MATE. |
| MATE-EHR v1 | Retransplant OR index LOS >90 d OR any dialysis after day 90 OR >=4 unplanned readmissions. | 7,250 | 26.393 | Primary definition to evaluate first. |
| MATE-EHR v2 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=4 unplanned readmissions. | 5,412 | 19.702 | Cleaner EHR-specific candidate. |
| MATE-EHR v3 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=5 unplanned readmissions. | 3,996 | 14.547 | Use if v1/v2 remain materially >30%. |
| High-burden MATE | >=2 non-death MATE components. | 1,996 | 7.266 | Secondary severity construct, not primary MATE. |

## Overlap / burden checks

| Overlap / burden calculation | N | % of cohort |
|---|---:|---:|
| Exactly 0 MATE-EHR components | 20,219 | 73.607 |
| Exactly 1 component | 5,254 | 19.127 |
| Exactly 2 components | 1,850 | 6.735 |
| Exactly 3 components | 138 | 0.502 |
| All 4 components | 8 | 0.029 |
| Retransplantation + any other component | 307 | 1.118 |
| Index LOS >90 d + dialysis after day 90 | 294 | 1.070 |
| Index LOS >90 d + high readmission burden | 97 | 0.353 |
| Dialysis after day 90 + high readmission burden | 1,505 | 5.479 |

## Decision rule

If MATE-EHR v1 is approximately 25-30%, retain the clinically natural >=4 readmission threshold. If it remains materially above 30%, evaluate v2 and then v3. Do not require >=2 components for primary MATE; use that only as a secondary high-burden phenotype.
