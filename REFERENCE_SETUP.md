# Reference Setup Guide

This guide walks you through setting up the reference files needed for the Drosophila TE Analysis Pipeline.

## 📋 Overview

The pipeline requires several reference files:
- **DM6 genome sequence** (FASTA format)
- **Gene annotations** (GTF format)
- **TE annotations** (GTF format)
- **TE consensus sequences** (Bergman library)
- **RepeatMasker output** (optional, for TE GTF creation)

## 🎯 Quick Setup

### Option 1: Automated Download (Recommended)

```bash
# Clone the pipeline
git clone https://github.com/yourusername/drosophila-te-pipeline.git
cd drosophila-te-pipeline

# Run the automated setup script
bash references/download_references.sh /path/to/your/reference/directory
```

### Option 2: Manual Setup

Follow the detailed instructions below for manual setup.

## 🔧 Manual Reference Setup

### Step 1: Create Reference Directory

```bash
# Create your reference directory
mkdir -p /path/to/references/dm6_bergman
cd /path/to/references/dm6_bergman
```

### Step 2: Download DM6 Genome

**From FlyBase (Recommended):**
```bash
# Download DM6 genome (BDGP6)
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2024_01/fasta/dmel-all-chromosome-r6.55.fasta.gz

# Decompress and rename
gunzip dmel-all-chromosome-r6.55.fasta.gz
mv dmel-all-chromosome-r6.55.fasta DM6_genome.fasta
```

**Alternative sources:**
```bash
# From NCBI
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz

# From UCSC
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
```

### Step 3: Download Gene Annotations

**From FlyBase:**
```bash
# Download gene annotations (GTF format)
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2024_01/gtf/dmel-all-r6.55.gtf.gz

# Decompress and rename
gunzip dmel-all-r6.55.gtf.gz
mv dmel-all-r6.55.gtf dmel_genes.gtf
```

### Step 4: Download Bergman TE Library

**From Bergman Lab GitHub:**
```bash
# Download TE consensus sequences
wget https://github.com/bergmanlab/transposons/raw/master/releases/D_mel_transposon_sequence_set_v10.1.fa

# Rename for consistency
mv D_mel_transposon_sequence_set_v10.1.fa D_mel_transposon_sequence_set.fa
```

### Step 5: Get TE Annotations

**Option A: Download pre-built TE GTF (Recommended):**
```bash
# Download from UCSC RepeatMasker track
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.out.gz
gunzip dm6.fa.out.gz
mv dm6.fa.out DM6_genome.fasta.out

# Download pre-built TE GTF
wget https://github.com/bergmanlab/mcclintock/raw/master/install/dm6/dm6_rmsk_TE.gtf
```

**Option B: Create TE GTF from RepeatMasker output:**
```bash
# If you have RepeatMasker output, convert to GTF
python references/create_te_gtf.py DM6_genome.fasta.out dm6_rmsk_TE.gtf
```

### Step 6: Add Processing Scripts

```bash
# Copy the processing scripts
cp path/to/pipeline/src/python/create_te_gtf.py .
cp path/to/pipeline/src/R/separate_counts.R .
```

## 🔍 Verify Reference Files

### Check File Integrity

```bash
# Verify all required files exist
ls -la
# Should see:
# DM6_genome.fasta
# dmel_genes.gtf
# dm6_rmsk_TE.gtf
# D_mel_transposon_sequence_set.fa
# DM6_genome.fasta.out (optional)

# Check file sizes (approximate)
du -h *.fasta *.gtf *.fa
# DM6_genome.fasta should be ~140MB
# dmel_genes.gtf should be ~50MB
# dm6_rmsk_TE.gtf should be ~15-20MB
# D_mel_transposon_sequence_set.fa should be ~1MB
```

### Validate File Formats

```bash
# Check FASTA format
head -5 DM6_genome.fasta
# Should start with: >2L or >X (chromosome names)

# Check GTF format
head -5 dmel_genes.gtf
# Should have 9 tab-separated columns

head -5 dm6_rmsk_TE.gtf
# Should have standard GTF format with TE information

# Count sequences
grep "^>" DM6_genome.fasta | wc -l
# Should be 7-8 sequences (chromosomes + mitochondrial)

grep "^>" D_mel_transposon_sequence_set.fa | wc -l
# Should be 100+ TE consensus sequences
```

### Check Chromosome Naming Consistency

```bash
# Extract chromosome names from genome
grep "^>" DM6_genome.fasta | head -10

# Extract chromosome names from GTF
cut -f1 dmel_genes.gtf | sort | uniq | head -10

# Extract chromosome names from TE GTF
cut -f1 dm6_rmsk_TE.gtf | sort | uniq | head -10
```

**Important:** Ensure chromosome names match between files!
- Common formats: `2L`, `2R`, `3L`, `3R`, `4`, `X`, `Y`
- Some files may have `chr` prefix: `chr2L`, `chr2R`, etc.

## 🛠 Fix Common Issues

### Chromosome Name Mismatches

**Remove 'chr' prefix if needed:**
```bash
# If TE GTF has 'chr' prefix but genome doesn't
sed 's/^chr//' dm6_rmsk_TE.gtf > dm6_rmsk_TE_fixed.gtf
mv dm6_rmsk_TE_fixed.gtf dm6_rmsk_TE.gtf

# If gene GTF has 'chr' prefix but genome doesn't
sed 's/^chr//' dmel_genes.gtf > dmel_genes_fixed.gtf
mv dmel_genes_fixed.gtf dmel_genes.gtf
```

**Add 'chr' prefix if needed:**
```bash
# If genome has 'chr' prefix but GTFs don't
sed 's/^/chr/' dm6_rmsk_TE.gtf > dm6_rmsk_TE_fixed.gtf
mv dm6_rmsk_TE_fixed.gtf dm6_rmsk_TE.gtf
```

### File Permission Issues

```bash
# Make sure files are readable
chmod 644 *.fasta *.gtf *.fa *.out
```

### Missing Features in GTF

**Check required GTF features:**
```bash
# Gene GTF should have these features
cut -f3 dmel_genes.gtf | sort | uniq -c
# Should include: gene, transcript, exon

# TE GTF should have these features
cut -f3 dm6_rmsk_TE.gtf | sort | uniq -c
# Should include: transposable_element or similar
```

## 📊 Reference Statistics

### Typical File Statistics

**DM6 Genome:**
- Size: ~140 MB
- Sequences: 7 main chromosomes (2L, 2R, 3L, 3R, 4, X, Y) + mitochondrial
- Total length: ~143 million bp

**Gene Annotations:**
- Genes: ~17,000 protein-coding genes
- Transcripts: ~30,000+ isoforms
- Exons: ~60,000+ exons

**TE Annotations:**
- TE copies: ~6,000-8,000 annotated elements
- TE families: ~100+ different families
- Coverage: ~20% of genome

### Generate Reference Summary

```bash
# Create a reference summary
cat > REFERENCE_SUMMARY.txt << EOF
# Reference Files Summary
Generated: $(date)

## Genome Statistics
Sequences: $(grep "^>" DM6_genome.fasta | wc -l)
Total length: $(grep -v "^>" DM6_genome.fasta | tr -d '\n' | wc -c) bp

## Gene Annotations
Genes: $(grep -c "gene" dmel_genes.gtf)
Transcripts: $(grep -c "transcript" dmel_genes.gtf)
Exons: $(grep -c "exon" dmel_genes.gtf)

## TE Annotations
TE elements: $(wc -l < dm6_rmsk_TE.gtf)
TE consensus: $(grep -c "^>" D_mel_transposon_sequence_set.fa)

## File Sizes
$(du -h *.fasta *.gtf *.fa)
EOF
```

## 🔄 Update References

### Keeping References Current

**Check for updates periodically:**
- FlyBase releases new annotations ~2-3 times per year
- Bergman TE library updates less frequently
- Monitor release notes for significant changes

**Version tracking:**
```bash
# Create version file
cat > REFERENCE_VERSIONS.txt << EOF
DM6 Genome: Release 6.55 ($(date))
Gene GTF: dmel-all-r6.55.gtf
TE Library: D_mel_transposon_sequence_set_v10.1.fa
TE GTF: dm6_rmsk_TE.gtf
RepeatMasker: $(date)
EOF
```

## 🎯 Quick Validation Script

Create a validation script to check everything:

```bash
#!/bin/bash
# validate_references.sh

echo "Validating reference files..."

# Required files
REQUIRED_FILES=(
    "DM6_genome.fasta"
    "dmel_genes.gtf"
    "dm6_rmsk_TE.gtf"
    "D_mel_transposon_sequence_set.fa"
)

# Check if files exist
for file in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        echo "✓ $file found"
    else
        echo "✗ $file missing"
        exit 1
    fi
done

# Check file formats
echo "Checking file formats..."

# FASTA check
if head -1 DM6_genome.fasta | grep -q "^>"; then
    echo "✓ Genome FASTA format OK"
else
    echo "✗ Genome FASTA format error"
fi

# GTF check
if head -1 dmel_genes.gtf | awk -F'\t' 'NF==9{exit 0}{exit 1}'; then
    echo "✓ Gene GTF format OK"
else
    echo "✗ Gene GTF format error"
fi

echo "Reference validation complete!"
```

## 🚨 Troubleshooting

### Common Problems

**"No such file or directory"**
- Check file paths in your configuration
- Verify files were downloaded completely
- Check file permissions

**"Chromosome not found"**
- Check chromosome naming consistency
- Use the chromosome fixing commands above

**"Empty alignment"**
- Verify GTF and FASTA chromosome names match
- Check GTF file format (9 columns, tab-separated)

**"Index building failed"**
- Ensure adequate disk space (>50GB free)
- Check memory availability (>64GB for STAR)
- Verify FASTA file integrity

### Getting Help

If you encounter issues:
1. Check the troubleshooting section above
2. Validate your files with the validation script
3. Check the pipeline logs for specific error messages
4. Submit an issue with your error details and file statistics

---

With properly set up reference files, your pipeline should run smoothly and produce accurate results!
