# Output Guide - Understanding Your Results

This guide explains all the outputs from the Drosophila TE Analysis Pipeline and how to interpret them.

## 📋 Overview

The pipeline generates several types of outputs:
- **Quality Control Reports** - Assessment of data quality
- **Count Matrices** - Gene and TE expression quantification
- **Differential Expression Results** - Statistical analysis of changes
- **Alignment Statistics** - Mapping and processing metrics
- **Summary Reports** - Comprehensive analysis overview

## 📁 Directory Structure

After a successful run, your results directory will look like this:

```
results/
├── qc/                              # Quality control outputs
│   ├── multiqc_report.html          # Comprehensive QC summary
│   ├── sample1_fastqc.html          # Per-sample FastQC reports
│   ├── sample1_fastqc.zip           # FastQC data files
│   └── sample1_flagstat.txt         # Alignment statistics
├── tetoolkit/                       # TEToolkit analysis results
│   ├── gene_te_counts_unstranded.txt    # Combined count matrix (unstranded)
│   ├── gene_te_counts_sense.txt         # Sense strand counts
│   ├── gene_te_counts_antisense.txt     # Antisense strand counts
│   ├── sample_info.txt                  # Sample metadata
│   ├── treated_vs_control_DESeq2.R      # DESeq2 analysis script
│   ├── treated_vs_control_gene_TE_analysis.txt   # DE results
│   ├── treated_vs_control_sigdiff_gene_TE.txt    # Significant changes
│   └── treated_vs_control.cntTable      # Raw TEcount output
├── differential_expression/         # Additional DE analysis
└── ANALYSIS_SUMMARY.md             # Complete pipeline summary
```

## 📊 Quality Control Reports

### MultiQC Report (`multiqc_report.html`)

**What it contains:**
- Aggregated quality metrics from all samples
- FastQC results summary
- Alignment statistics overview
- Adapter content analysis
- Sequence duplication levels

**How to interpret:**
1. **Open in web browser** for interactive exploration
2. **General Statistics table** shows key metrics per sample
3. **FastQC sections** highlight quality issues
4. **Status checks** indicate pass/warn/fail for each metric

**Key metrics to check:**
- **% Dups:** <30% is good, >50% may indicate issues
- **% GC:** Should be ~42% for Drosophila
- **M Seqs:** Number of read pairs (millions)
- **Length:** Read length after trimming

### Individual FastQC Reports

**Per-sample HTML reports showing:**
- **Per base sequence quality** - Should be >28 (green zone)
- **Per sequence quality scores** - Peak should be >28
- **Per base N content** - Should be near 0%
- **Sequence length distribution** - Should be uniform after trimming
- **Adapter content** - Should be <1% after trimming
- **Overrepresented sequences** - May indicate contamination

### Alignment Statistics (`*_flagstat.txt`)

**Example output:**
```
1000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
850000 + 0 mapped (85.00% : N/A)
1000000 + 0 paired in sequencing
500000 + 0 read1
500000 + 0 read2
800000 + 0 properly paired (80.00% : N/A)
840000 + 0 with mate mapped to different chr
20000 + 0 with mate mapped to different chr (mapQ>=5)
```

**Interpretation:**
- **Mapped %:** >80% is good for Drosophila
- **Properly paired %:** >75% indicates good library quality
- **Different chr:** High values may indicate genomic rearrangements

## 🧮 Count Matrices

### Combined Gene+TE Count Matrix

**File:** `gene_te_counts_unstranded.txt` (or strand-specific versions)

**Format:**
```
gene_id                sample1    sample2    sample3    sample4
FBgn0000001           1234       1456       1123       1367
FBgn0000002           567        678        567        623
FBgn0000003           89         112        95         108
...
roo                   45         23         67         34
BEL                   12         8          15         9
HMS-Beagle           78         45         89         56
...
```

**Understanding the data:**
- **Rows:** Gene IDs (FBgn...) and TE family names
- **Columns:** Your samples
- **Values:** Raw read counts
- **Gene section:** Annotated protein-coding genes (~17,000)
- **TE section:** Transposable element families (~100+)

### Strand-Specific Count Matrices

**Purpose:** Some TE insertions may be orientation-specific

**Files:**
- `gene_te_counts_sense.txt` - Reads mapping to sense strand
- `gene_te_counts_antisense.txt` - Reads mapping to antisense strand

**When to use:**
- If your library prep is strand-specific
- For studying TE insertion orientation effects
- For antisense transcription analysis

### Sample Information File

**File:** `sample_info.txt`

**Content:**
```
sample          condition    timepoint    replicate
control_rep1    control      NA           1
control_rep2    control      NA           2
treated_rep1    treated      NA           1
treated_rep2    treated      NA           2
```

## 📈 Differential Expression Results

### DESeq2 Analysis Script

**File:** `*_DESeq2.R`

**What it does:**
- Loads count data and sample information
- Performs DESeq2 normalization
- Tests for differential expression
- Generates plots and statistics

**To run manually:**
```bash
cd results/tetoolkit/
Rscript treated_vs_control_DESeq2.R
```

### Main Results Table

**File:** `*_gene_TE_analysis.txt`

**Columns explained:**
```
gene_id         - Gene or TE identifier
baseMean        - Average expression across all samples
log2FoldChange  - Log2 fold change (treatment/control)
lfcSE          - Standard error of log2 fold change
stat           - Wald test statistic
pvalue         - Raw p-value
padj           - Adjusted p-value (FDR corrected)
```

**Example:**
```
gene_id         baseMean    log2FoldChange  lfcSE      stat       pvalue      padj
FBgn0000001     1245.6      -0.45          0.12       -3.75      1.8e-04     0.023
roo             89.3        2.34           0.45       5.20       2.0e-07     1.2e-05
BEL             156.8       -1.67          0.38       -4.39      1.1e-05     0.0089
```

**Interpretation:**
- **log2FoldChange > 0:** Upregulated in treatment
- **log2FoldChange < 0:** Downregulated in treatment
- **padj < 0.05:** Statistically significant (typical cutoff)
- **|log2FoldChange| > 1:** Biologically significant (2-fold change)

### Significant Results

**File:** `*_sigdiff_gene_TE.txt`

**Pre-filtered results showing only:**
- **padj < 0.05** (statistically significant)
- **|log2FoldChange| > 1** (biologically significant)
- **baseMean > 10** (sufficiently expressed)

### Raw TEcount Output

**File:** `*.cntTable`

**Direct output from TEToolkit before DESeq2 processing**
- Same format as count matrices
- Includes additional TEToolkit metadata
- Useful for custom analysis workflows

## 📋 Analysis Summary

### Main Summary Report

**File:** `ANALYSIS_SUMMARY.md`

**Contains:**
- Pipeline configuration details
- Sample processing status
- Key output file locations
- Analysis approach summary
- Quality metrics overview

**Example sections:**
```markdown
## Sample Information
| Sample | Condition | Status |
|--------|-----------|--------|
| ctrl_1 | control   | ✓ Completed |
| treat_1| treated   | ✓ Completed |

## Key Results
- **Total genes analyzed:** 17,234
- **Total TEs analyzed:** 145
- **Significantly changed genes:** 1,247
- **Significantly changed TEs:** 23
```

## 🔍 Interpreting Results

### Gene Expression Changes

**Typical observations:**
- **Stress response genes** often upregulated
- **Metabolic genes** may show complex patterns
- **Developmental genes** context-dependent

**Follow-up analysis:**
- Gene ontology enrichment
- Pathway analysis
- Expression clustering

### TE Expression Changes

**Common patterns:**
- **Heat shock:** Many TEs activated
- **Developmental transitions:** Specific TE families
- **Stress conditions:** DNA transposons vs retrotransposons

**TE-specific considerations:**
- **Recent insertions** vs ancient copies
- **Family-level** vs subfamily analysis
- **Genomic location** effects

### Quality Assessment

**Good quality indicators:**
- **Alignment rate >80%**
- **Even coverage across gene bodies**
- **Biological replicates cluster together**
- **Expected number of expressed genes (~12,000-15,000)**

**Potential issues:**
- **Low alignment rates** - Check reference genome
- **High duplication** - Library complexity issues
- **Outlier samples** - Technical problems
- **Unexpected TE activation** - Stress or contamination

## 📊 Downstream Analysis

### Using Count Matrices

**Load into R:**
```r
# Load count data
counts <- read.table("gene_te_counts_unstranded.txt", 
                    header=TRUE, row.names=1, sep="\t")

# Separate genes and TEs
gene_counts <- counts[grepl("^FBgn", rownames(counts)), ]
te_counts <- counts[!grepl("^FBgn", rownames(counts)), ]

# Basic statistics
summary(rowSums(gene_counts))
summary(rowSums(te_counts))
```

**Load into Python:**
```python
import pandas as pd

# Load count data
counts = pd.read_csv("gene_te_counts_unstranded.txt", 
                    sep="\t", index_col=0)

# Separate genes and TEs
gene_counts = counts[counts.index.str.startswith('FBgn')]
te_counts = counts[~counts.index.str.startswith('FBgn')]

# Basic statistics
print(gene_counts.sum(axis=1).describe())
print(te_counts.sum(axis=1).describe())
```

### Custom Visualizations

**R plotting examples:**
```r
library(ggplot2)
library(DESeq2)

# MA plot
plotMA(dds, alpha=0.05)

# PCA plot
plotPCA(vst(dds), intgroup="condition")

# Heatmap of top TEs
library(pheatmap)
top_tes <- head(te_results[order(te_results$padj), ], 20)
pheatmap(te_counts[rownames(top_tes), ], scale="row")
```

### Functional Analysis

**Gene ontology enrichment:**
```r
library(clusterProfiler)
library(org.Dm.eg.db)

# Get significant genes
sig_genes <- rownames(results[results$padj < 0.05 & 
                             abs(results$log2FoldChange) > 1, ])

# GO enrichment
ego <- enrichGO(gene = sig_genes,
                OrgDb = org.Dm.eg.db,
                keyType = 'FLYBASE',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)
```

## 🚨 Troubleshooting Results

### Low Gene Detection

**Possible causes:**
- Poor RNA quality
- Low sequencing depth
- Alignment issues
- Wrong strand specificity setting

**Solutions:**
- Check FastQC reports
- Verify alignment rates
- Try different TEToolkit strand settings

### Unexpected TE Activation

**Possible causes:**
- Heat shock during sample processing
- Stress conditions
- Contamination
- Age-related activation

**Investigation:**
- Check sample processing notes
- Compare to controls
- Examine specific TE families
- Validate with qPCR

### No Differential Expression

**Possible causes:**
- Insufficient biological effect
- Low statistical power
- Technical variation
- Wrong comparison groups

**Solutions:**
- Check experimental design
- Examine raw count distributions
- Consider different cutoffs
- Validate key targets

### High Variability

**Possible causes:**
- Batch effects
- Sample quality differences
- Developmental timing
- Environmental factors

**Solutions:**
- Check sample clustering
- Consider batch correction
- Examine individual sample QC
- Review experimental protocol

## 📚 Additional Resources

### File Format References
- **GTF format:** [Ensembl GTF specification](http://www.ensembl.org/info/website/upload/gff.html)
- **SAM/BAM:** [SAMtools documentation](http://samtools.github.io/hts-specs/)
- **FastQC:** [FastQC help](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

### Statistical Methods
- **DESeq2:** [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- **TEToolkit:** [TEToolkit paper](https://academic.oup.com/bioinformatics/article/31/22/3593/240926)

### Drosophila Resources
- **FlyBase:** [Gene information](http://flybase.org/)
- **modENCODE:** [Expression data](http://www.modencode.org/)
- **Bergman Lab:** [TE resources](http://bergmanlab.genetics.uga.edu/)

---

This guide should help you understand and interpret all the outputs from your TE analysis pipeline. For specific questions about your results, consider the biological context of your experiment and consult the relevant literature for similar studies.
