# Variant Population Frequency Analysis

R-based pipeline for analyzing population-specific variant frequencies
using gnomAD data, with statistical enrichment testing and visualization.

---

## What it does
- Parses and filters gnomAD population frequency data
- Merges with disease-specific variant lists and ACMG classifications
- Calculates Minor Allele Frequency (MAF) across 11 global populations
- Performs Fisher's Exact Test with Bonferroni correction per population
- Generates publication-quality bubble plots showing population-specific enrichment

---

## Populations Analyzed

| Code | Population |
|---|---|
| AFR | African |
| AMR | American |
| ASJ | Ashkenazi Jewish |
| EAS | East Asian |
| FIN | Finnish |
| MID | Middle Eastern |
| NFE | Non-Finnish European |
| AMI | Amish |
| SAS | South Asian |
| OTH | Other |

---

## Input Files

### 1. `all_variant.txt`
gnomAD variant table with population-level AC and AN columns.

Expected columns:
```
variant_id  gene_symbol  pos  hgvs  hgvsp  consequence
joint_ac  joint_an  joint_afr_ac  joint_afr_an ...
```

### 2. `variant_list.txt`
Disease-specific variant list with ACMG classifications.

Expected columns:
```
VCF  AC  AN  ACMG
```

---

## Workflow

```
gnomAD variant table (all_variant.txt)
        ↓
Filter AC and AN columns
        ↓
Merge with disease variant list (variant_list.txt)
        ↓
Calculate population-specific MAF
        ↓
Fisher's Exact Test + Bonferroni correction
        ↓
Bubble plot visualization (Figure.png)
```

---

## How to Run

### 1. Install R dependencies
```r
install.packages(c("dplyr", "ggplot2"))
```

### 2. Update working directory in script
```r
setwd("/path/to/your/data")
```

### 3. Run the script
```r
Rscript maf_analysis.R
```

---

## Output

- Population-specific MAF values per variant
- Fisher's Exact Test p-values (Bonferroni corrected)
- Bubble plot PNG:
  - Bubble size → MAF
  - Color → ACMG classification
  - Red circle → Statistically significant enrichment (p < 0.05)

---

## Plot Description

| Visual Element | Meaning |
|---|---|
| X-axis | Variants |
| Y-axis | Populations |
| Bubble size | Minor Allele Frequency (MAF) |
| Color (pink) | Likely Pathogenic |
| Color (teal) | Pathogenic |
| Red circle | Bonferroni-corrected p < 0.05 |

---

## Requirements

- R ≥ 4.0
- R packages: `dplyr`, `ggplot2`
- gnomAD variant data (v2/v3/v4)

---

## Notes

- Population AN absent counts calculated automatically
- Variants with zero AN excluded from Fisher's Exact Test
- Follows ACMG/AMP variant classification standards
- Script designed for hg38 variant coordinates

---

## Author

Treesa Issen | CSIR-IGIB, New Delhi
[ORCID](https://orcid.org/0009-0001-9056-7101)
