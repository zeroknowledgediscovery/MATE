# MATE definition used in the manuscript

The current manuscript uses one operational definition of **Major Adverse Transplant Episode (MATE)** in Epic COSMOS.

A recipient meets **MATE** if **any** of the following occurs within the first year after lung transplantation:

| Component | COSMOS operational definition |
|---|---|
| Retransplantation | Retransplantation within 365 days of the index lung transplant |
| Prolonged index hospitalization | Index transplant hospitalization length of stay >90 days |
| Dialysis after day 90 | At least one dialysis/hemodialysis event during days 91–365 |
| Recurrent unplanned readmission | At least 4 unplanned inpatient readmissions through day 365 |
| Persistent oxygen dependence | Oxygen dependence / oxygen DME evidence on or after day 90 |

Thus:

```text
MATE =
    retransplantation within 365 days
    OR index LOS > 90 days
    OR any dialysis/hemodialysis during days 91–365
    OR >= 4 unplanned inpatient readmissions through day 365
    OR oxygen dependence / oxygen DME evidence on or after day 90
```

This is the five-component, oxygen-inclusive phenotype previously labeled `MATE-EHR v1+` during development. The manuscript now refers to it simply as **MATE**.

Death is not included in the present COSMOS EHR phenotype because mortality ascertainment is not sufficiently reliable in COSMOS. The manuscript treats incorporation of death into the broader clinical construct as a future validation step rather than as a second competing MATE definition.
