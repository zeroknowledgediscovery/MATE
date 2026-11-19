# COSMOS MATE Definition Worksheet

### **Cohort summary**
* Date of query: November 17, 2026
* Eligible patients: 24,213
* Eligible index dates range: October 07, 2005 to November 16, 2025.

**Working plausibility check:** approximately 25-30% non-death MATE-EHR prevalence in COSMOS; the final clinical MATE is expected to be approximately 30-35% after death is incorporated in LTC/STAR and overlap is accounted for. These ranges are for plausibility assessment, not prevalence calibration.

| Candidate component / definition | Exact COSMOS operational definition | N | Prevalence (%) | Decision / notes |
|---|---|---:|---:|---|
| Retransplantation | Retransplantation within 365 days of index lung transplant. | 303 | 1.251394 | Core MATE component; retain. |
| Index hospitalization >90 d | Index transplant admission LOS >90 days. | 650 | 2.684508 | Core MATE component; add to COSMOS proxy. |
| Dialysis after day 90 - any | >=1 dialysis/hemodialysis event on days 91-365. | 3,968 | 16.387891 | Closest implementation of proposed MATE. |
| Dialysis after day 90 - persistent | >=2 dialysis events on distinct dates during days 91-365, preferably >=7 days apart. | 2,127 | 8.784537 | Dialysis records on >= 3 different days during 90-365 days post-index |
| >=3 unplanned readmissions | >=3 unplanned inpatient readmissions after index discharge through day 365. | 5,137 | 21.215876 | Current preliminary threshold; likely permissive. |
| >=4 unplanned readmissions | >=4 unplanned inpatient readmissions after index discharge through day 365. | 3,109 | 12.840210 | Matches proposal's >3 readmissions; test first. |
| >=5 unplanned readmissions | >=5 unplanned inpatient readmissions after index discharge through day 365. | 1,852 | 7.648784 | Main tightening option if >=4 is too prevalent. |
| >=6 unplanned readmissions | >=6 unplanned inpatient readmissions after index discharge through day 365. | 1,126 | 4.650394 | Sensitivity analysis. |
| Oxygen dependence after day 90 | Oxygen dependence / oxygen DME evidence on or after day 90. | 5,952 | 24.581836 | Sensitivity/QOL candidate; exclude from primary MATE. |
| MATE-EHR v1 | Retransplant OR index LOS >90 d OR any dialysis after day 90 OR >=4 unplanned readmissions. | 6,269 | 25.891050 | Primary definition to evaluate first. |
| MATE-EHR v1+ | MATE-EHR v1 OR oxygen dependence after day 90. | 9,799 | 40.469995 | Expanded v1 sensitivity definition including oxygen dependence. |
| MATE-EHR v2 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=4 unplanned readmissions. | 5,040 | 20.815265 | Cleaner EHR-specific candidate. |
| MATE-EHR v3 | Retransplant OR index LOS >90 d OR persistent dialysis after day 90 OR >=5 unplanned readmissions. | 4,043 | 16.697642 | Use if v1/v2 remain materially >30%. |
| High-burden MATE | >=2 non-death MATE components. | 1,635 | 6.752571 | Secondary severity construct, not primary MATE. |

## Overlap / burden checks

| Overlap / burden calculation | N | % of cohort |
|---|---:|---:|
| Exactly 0 MATE-EHR components | 15,304 | 63.205716 |
| Exactly 1 component | 4,203 | 17.358444 |
| Exactly 2 components | 1,546 | 6.385000 |
| Exactly 3 components | 86 | 0.355181 |
| All 4 components | 3 | 0.012390 |
| Retransplantation + any other component | 126 | 0.520382 |
| Index LOS >90 d + dialysis after day 90 | 371 | 1.532235 |
| Index LOS >90 d + high readmission burden | 111 | 0.458431 |
| Dialysis after day 90 + high readmission burden | 1,218 | 5.030356 |

## Decision rule

If MATE-EHR v1 is approximately 25-30%, retain the clinically natural >=4 readmission threshold. If it remains materially above 30%, evaluate v2 and then v3. Do not require >=2 components for primary MATE; use that only as a secondary high-burden phenotype.
