# User Guide - Drosophila TE Analysis Pipeline

## Table of Contents

1. [Getting Started](#getting-started)
2. [Pre-Analysis Setup](#pre-analysis-setup)  
3. [Sample Configuration](#sample-configuration)
4. [Pipeline Configuration](#pipeline-configuration)
5. [Running the Pipeline](#running-the-pipeline)
6. [Understanding Results](#understanding-results)
7. [Advanced Usage](#advanced-usage)
8. [Troubleshooting](#troubleshooting)

## Getting Started

### System Requirements

**Minimum Requirements:**
- 16 CPU cores
- 128GB RAM (for STAR index building)
- 500GB available storage
- SLURM job scheduler

**Recommended Requirements:**
- 32+ CPU cores  
- 256GB+ RAM
- High-speed SSD storage
- Dedicated compute node access

### Environment Setup

1. **Load required modules** (example for common HPC systems):
```bash
module purge
module load star/2.7.11b
module load samtools/1.9  
module load R/4.2.2
module load fastqc/0.12.1
module load cutadapt/3.5
module load bedtools/2.29.0
module load MultiQC/1.9
module load TEToolkit/2.0.3
```

2. **Set up R environment**:
```bash
# Create personal R library
mkdir -p ~/R/library
export R_LIBS_USER="$HOME/R/library"

# Install required packages
Rscript -e "install.packages(c('data.table', 'ggplot2', 'dplyr'), repos='https://cran.r-project.org', lib='~/R/library')"
```

## Pre-Analysis Setup

### 1. Reference Files Preparation

Create a reference directory structure:

```bash
mkdir -p /path/to/your/references/dm6_bergman
cd /path/to/your/references/dm6_bergman
```

**Required files:**
- `DM6_genome.fasta` - *D. melanogaster* genome (DM6/BDGP6)
- `dmel_genes.gtf` - Gene annotations 
- `DM6_genome.fasta.out` - RepeatMasker output
- `D_mel_transposon_sequence_set.fa` - Bergman TE library
- `dm6_rmsk_TE.gtf` - TE annotations in GTF format
- `create_te_gtf.py` - TE GTF creation script (if needed)
- `seperate_counts.R` - Results processing script

### 2. Download Reference Files

```bash
# DM6 genome from FlyBase
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2024_01/fasta/dmel-all-chromosome-r6.55.fasta.gz
gunzip dmel-all-chromosome-r6.55.fasta.gz
mv dmel-all-chromosome-r6.55.fasta DM6_genome.fasta

# Gene annotations
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2024_01/gtf/dmel-all-r6.55.gtf.gz  
gunzip dmel-all-r6.55.gtf.gz
mv dmel-all-r6.55.gtf dmel_genes.gtf

# Bergman TE library (from Bergman Lab)
wget https://github.com/bergmanlab/transposons/raw/master/releases/D_mel_transposon_sequence_set_v10.1.fa
mv D_mel_transposon_sequence_set_v10.1.fa D_mel_transposon_sequence_set.fa
```

### 3. Project Directory Setup

```bash
# Clone the pipeline
git clone https://github.com/yourusername/drosophila-te-pipeline.git
cd drosophila-te-pipeline

# Create project-specific directory
mkdir -p projects/my_experiment
cd projects/my_experiment

# Copy pipeline script
cp ../../scripts/tetoolkit_pipeline.sh .
```

## Sample Configuration

### Sample Information Format

The pipeline supports two sample format styles:

**Simple format** (4 fields):
```bash
"sample_name:condition:read1_path:read2_path"
```

**Extended format** (6 fields):
```bash  
"sample_name:condition:timepoint:replicate:read1_path:read2_path"
```

### Example Configurations

**Simple comparison (Control vs Treatment):**
```bash
SAMPLES=(
    "ctrl_rep1:control:/data/fastq/ctrl_rep1_R1.fastq.gz:/data/fastq/ctrl_rep1_R2.fastq.gz"
    "ctrl_rep2:control:/data/fastq/ctrl_rep2_R1.fastq.gz:/data/fastq/ctrl_rep2_R2.fastq.gz"
    "ctrl_rep3:control:/data/fastq/ctrl_rep3_R1.fastq.gz:/data/fastq/ctrl_rep3_R2.fastq.gz"
    "treat_rep1:treated:/data/fastq/treat_rep1_R1.fastq.gz:/data/fastq/treat_rep1_R2.fastq.gz"
    "treat_rep2:treated:/data/fastq/treat_rep2_R1.fastq.gz:/data/fastq/treat_rep2_R2.fastq.gz"
    "treat_rep3:treated:/data/fastq/treat_rep3_R1.fastq.gz:/data/fastq/treat_rep3_R2.fastq.gz"
)
```

**Time course experiment:**
```bash
SAMPLES=(
    "ctrl_0h_rep1:control:0h:1:/data/fastq/ctrl_0h_rep1_R1.fastq.gz:/data/fastq/ctrl_0h_rep1_R2.fastq.gz"
    "ctrl_0h_rep2:control:0h:2:/data/fastq/ctrl_0h_rep2_R1.fastq.gz:/data/fastq/ctrl_0h_rep2_R2.fastq.gz"
    "ctrl_4h_rep1:control:4h:1:/data/fastq/ctrl_4h_rep1_R1.fastq.gz:/data/fastq/ctrl_4h_rep1_R2.fastq.gz"
    "ctrl_4h_rep2:control:4h:2:/data/fastq/ctrl_4h_rep2_R1.fastq.gz:/data/fastq/ctrl_4h_rep2_R2.fastq.gz"
    "treat_0h_rep1:treated:0h:1:/data/fastq/treat_0h_rep1_R1.fastq.gz:/data/fastq/treat_0h_rep1_R2.fastq.gz"
    "treat_0h_rep2:treated:0h:2:/data/fastq/treat_0h_rep2_R1.fastq.gz:/data/fastq/treat_0h_rep2_R2.fastq.gz"
    "treat_4h_rep1:treated:4h:1:/data/fastq/treat_4h_rep1_R1.fastq.gz:/data/fastq/treat_4h_rep1_R2.fastq.gz"
    "treat_4h_rep2:treated:4h:2:/data/fastq/treat_4h_rep2_R1.fastq.gz:/data/fastq/treat_4h_rep2_R2.fastq.gz"
)
```

## Pipeline Configuration

### 1. Edit Main Configuration Section - Note that reference files can be changed, but you must thoroughly check compatibility. GTF annotations and reference genomes must match, and annotation format must be consistent with that use by Casey Bergam Labs (TEtoolkit compatible)
Repeat Masker out file not necessary when using pre-computed TE-GTF.

Open `tetoolkit_pipeline.sh` and modify the USER CONFIGURATION section:

```bash
# Project settings
PROJECT_NAME="my_te_experiment"

# Reference files - UPDATE THESE PATHS
REFERENCE_DIR="/path/to/your/references/dm6_bergman"
GENOME_FILE="$REFERENCE_DIR/DM6_genome.fasta"
GTF_FILE="$REFERENCE_DIR/dmel_genes.gtf"
REPEATMASKER_FILE="$REFERENCE_DIR/DM6_genome.fasta.out"
BERGMAN_TE_LIBRARY="$REFERENCE_DIR/D_mel_transposon_sequence_set.fa"
PYTHON_SCRIPT="$REFERENCE_DIR/create_te_gtf.py"
R_SCRIPT="$REFERENCE_DIR/seperate_counts.R"
```

### 2. Configure SLURM Resources

Adjust SLURM parameters based on your system:

```bash
#SBATCH --job-name=my_te_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1  
#SBATCH --cpus-per-task=32     # Increase for faster processing
#SBATCH --mem=256G             # Increase if you have memory issues
#SBATCH --time=24:00:00        # Adjust based on sample count
#SBATCH --partition=general    # Use your HPC partition name
#SBATCH --qos=general          # Use your QOS settings
```

### 3. TEToolkit Options

Configure TEToolkit behavior:

```bash
# TEToolkit specific options
TETOOLKIT_MODE="multi"      # multi, uniq - how to handle multi-mapping reads
TETOOLKIT_FORMAT="BAM"      # BAM or SAM
TETOOLKIT_STRANDED="yes"    # no, yes, reverse - match your library prep
```

**Strand specificity guide:**
- `no` - Unstranded libraries (most common)
- `yes` - Forward stranded (dUTP method)  
- `reverse` - Reverse stranded (typical Illumina TruSeq)

### 4. Step Control Toggles

Control which parts of the pipeline to run:

```bash
RUN_MODULES=1        # Load modules and verify tools
RUN_VALIDATE=1       # Validate input files
RUN_SETUP=1          # Setup project directories  
RUN_REFERENCES=1     # Prepare reference files and build indices
RUN_SAMPLES=1        # Process samples (QC, trimming, alignment)
RUN_TETOOLKIT=1      # Run TEToolkit analysis
RUN_RESULTS=1        # Process TEToolkit results
RUN_REPORTS=1        # Generate final reports
```

## Running the Pipeline

### 1. Initial Validation Run

First, run a validation-only check:

```bash
# Set only validation steps
RUN_MODULES=1
RUN_VALIDATE=1  
RUN_SETUP=1
# Set all others to 0

sbatch tetoolkit_pipeline.sh
```

### 2. Reference Preparation (One-time)

Build references (this takes the longest):

```bash
RUN_MODULES=1
RUN_VALIDATE=0  # Skip if already validated
RUN_SETUP=0     # Skip if already done  
RUN_REFERENCES=1
# Set all others to 0

sbatch tetoolkit_pipeline.sh
```

### 3. Full Pipeline Run

After references are built, run the complete analysis:

```bash
# Set all toggles to 1, or selectively run remaining steps
RUN_MODULES=0       # Skip if modules already loaded
RUN_VALIDATE=0      # Skip if already validated
RUN_SETUP=0         # Skip if already setup
RUN_REFERENCES=0    # Skip if references already built
RUN_SAMPLES=1       # Run sample processing
RUN_TETOOLKIT=1     # Run TEToolkit analysis  
RUN_RESULTS=1       # Process results
RUN_REPORTS=1       # Generate reports

sbatch tetoolkit_pipeline.sh
```

### 4. Monitor Progress

```bash
# Check job status
squeue -u $USER

# Monitor log files
tail -f te_analysis_[JOBID].out
tail -f te_analysis_[JOBID].err

# Check specific step logs
tail -f logs/*_cutadapt.log
tail -f logs/*_star.log
```

## Understanding Results

### Primary Output Files

1. **Quality Control Reports**
   - `results/qc/multiqc_report.html` - Comprehensive QC summary
   - `results/qc/*_fastqc.html` - Individual sample QC

2. **Count Matrices** 
   - `results/tetoolkit/gene_te_counts_unstranded.txt` - Combined gene+TE counts
   - `results/tetoolkit/gene_te_counts_sense.txt` - Sense strand specific
   - `results/tetoolkit/gene_te_counts_antisense.txt` - Antisense strand specific

3. **Differential Expression**
   - `results/tetoolkit/*_vs_*_DESeq2_*` - DESeq2 analysis results
   - `results/tetoolkit/*_sigdiff_*` - Significantly changed genes/TEs

4. **Summary**
   - `results/ANALYSIS_SUMMARY.md` - Complete pipeline summary

### Count Matrix Format

The TEToolkit count matrices contain both genes and TEs:

```
gene_id/TE_id    sample1    sample2    sample3    ...
FBgn0000001      1234       1456       1123       ...
FBgn0000002      567        678        567        ...
...
ELEMENT1         12         23         15         ...
ELEMENT2         45         67         34         ...
```

### Quality Metrics to Check

1. **Alignment rates** (>80% typically good for Drosophila)
2. **Duplication rates** (<30% preferred)  
3. **Gene body coverage** (should be relatively even)
4. **TE detection rates** (varies by sample type)

## Advanced Usage

### Custom Differential Expression Comparisons

For complex experimental designs, you may want to customize the TEtranscripts comparisons:

1. **Edit the comparison function:**
```bash
# In the script, modify run_tetranscripts_comparisons() function
# Add custom comparison logic based on your experimental design
```

2. **Manual TEtranscripts runs:**
```bash
# After TEcount finishes, run custom comparisons manually
TEtranscripts --format BAM \
             --mode multi \
             --stranded yes \
             --project "custom_comparison" \
             --GTF references/genes.gtf \
             --TE references/te_annotation.gtf \
             -t treatment_sample1.bam,treatment_sample2.bam \
             -c control_sample1.bam,control_sample2.bam
```

### Batch Processing Multiple Experiments

For processing multiple experiments:

1. **Create experiment-specific directories:**
```bash
mkdir -p projects/experiment1 projects/experiment2
```

2. **Use different configuration files:**
```bash
# Copy and modify pipeline for each experiment
cp scripts/tetoolkit_pipeline.sh projects/experiment1/
cp scripts/tetoolkit_pipeline.sh projects/experiment2/
```

3. **Submit as job array:**
```bash
# Create job array script for multiple experiments
sbatch --array=1-2 submit_multiple_experiments.sh
```

### Memory Optimization

For memory-constrained systems:

1. **Reduce STAR index parameters:**
```bash
# In the STAR index building section, add:
--genomeSAindexNbases 12   # Reduce for smaller genomes
--limitGenomeGenerateRAM $((64 * 1024 * 1024 * 1024))  # 64GB limit
```

2. **Process samples in batches:**
```bash
# Split large sample sets into smaller batches
# Run samples 1-5, then 6-10, etc.
```

## Troubleshooting

### Common Error Messages

**"TEcount not found"**
```bash
# Solution: Check module loading
module list
module avail | grep -i tetoolkit
module load TEToolkit/2.0.3
```

**"STAR index error"**
```bash
# Solution: Check available memory and disk space
df -h .  # Check disk space
free -h  # Check available memory
# Increase --mem in SLURM header if needed
```

**"No such file or directory"**
```bash
# Solution: Check all file paths in configuration
ls -la /path/to/your/references/
# Verify all reference files exist and are readable
```

### Performance Issues

**Slow STAR alignment:**
- Increase `--cpus-per-task` 
- Use local/fast storage for reference files
- Consider using STAR's shared memory feature

**High memory usage:**
- Monitor with `top` or `htop`
- Reduce number of concurrent processes
- Use memory-optimized STAR parameters

**Long queue times:**
- Check with your HPC administrator about resource availability
- Consider running during off-peak hours
- Break large jobs into smaller components

### Recovery from Failed Jobs

The pipeline supports resuming from any step:

1. **Identify where the failure occurred**
2. **Set failed step toggle to 1, all previous steps to 0**
3. **Fix the underlying issue**
4. **Resubmit the job**

Example recovery from alignment failure:
```bash
# If alignment failed, set:
RUN_MODULES=0
RUN_VALIDATE=0
RUN_SETUP=0
RUN_REFERENCES=0
RUN_SAMPLES=1      # Restart from sample processing
RUN_TETOOLKIT=1    # Continue with remaining steps
RUN_RESULTS=1
RUN_REPORTS=1
```

### Getting Help

1. **Check log files first:**
   - Main pipeline log: `te_analysis_[JOBID].out` 
   - Error log: `te_analysis_[JOBID].err`
   - Step-specific logs: `logs/`

2. **Verify configuration:**
   - All file paths are correct
   - Sample information is properly formatted
   - SLURM parameters are appropriate for your system

3. **Contact support:**
   - GitHub Issues for bug reports
   - Include relevant log snippets and configuration details
   - Specify your HPC system and module versions

---

This user guide covers the essential aspects of running the Drosophila TE analysis pipeline. For additional questions or advanced customizations, please refer to the GitHub repository documentation or submit an issue.
