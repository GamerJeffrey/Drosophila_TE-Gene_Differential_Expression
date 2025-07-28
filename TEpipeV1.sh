#!/bin/bash
#SBATCH --job-name=drosophila_te_tetoolkit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=te_analysis_%j.out
#SBATCH --error=te_analysis_%j.err

# Drosophila TE Analysis Pipeline with TEToolkit
# Uses Casey Bergman's TE library + TEToolkit for integrated TE/gene analysis

set -euo pipefail

# =============================================================================
# USER CONFIGURATION - EDIT THESE PATHS
# =============================================================================

# Project settings
PROJECT_NAME="tetoolkit_te_analysis"

# Reference files - DM6 with Bergman library (local directory)
REFERENCE_DIR="/home/FCAM/minaba/DM6_RNAseq_TE_pipeline_data"
GENOME_FILE="$REFERENCE_DIR/DM6_genome.fasta"
GTF_FILE="$REFERENCE_DIR/dmel_genes.gtf"
REPEATMASKER_FILE="$REFERENCE_DIR/DM6_genome.fasta.out"
BERGMAN_TE_LIBRARY="$REFERENCE_DIR/D_mel_transposon_sequence_set.fa"
PYTHON_SCRIPT="$REFERENCE_DIR/create_te_gtf.py"
R_SCRIPT="$REFERENCE_DIR/seperate_counts.R"
PLOTS_SCRIPT="$REFERENCE_DIR/create_deseq2_plots.R"

# Local reference directory
LOCAL_REF_DIR="references"

# Sample definitions (MODIFIED FOR YOUR SAMPLES)
# Format: "sample_name:condition:timepoint:replicate:read1_path:read2_path"
SAMPLES=(
    "1E_rep1:treated:baseline:1:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/1E_S1154_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/1E_S1154_L004_R2_001.fastq.gz"
    "2E_rep2:treated:baseline:2:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/2E_S1155_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/2E_S1155_L004_R2_001.fastq.gz"
    "3E_rep3:treated:baseline:3:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/3E_S1156_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/3E_S1156_L004_R2_001.fastq.gz"
    "4C_rep1:control:baseline:1:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/4C_S1157_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/4C_S1157_L004_R2_001.fastq.gz"
    "5C_rep2:control:baseline:2:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/5C_S1158_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/5C_S1158_L004_R2_001.fastq.gz"
    "6C_rep3:control:baseline:3:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/6C_S1159_L004_R1_001.fastq.gz:/home/FCAM/minaba/JeffGamer/RNAseqStella/Inaba_9mRNA_08Nov2024/6C_S1159_L004_R2_001.fastq.gz"
)

# =============================================================================
# STEP TOGGLES - SET TO 1 TO RUN, 0 TO SKIP
# =============================================================================
RUN_MODULES=1        # Load modules and verify tools
RUN_VALIDATE=0       # Validate input files
RUN_SETUP=0          # Setup project directories
RUN_REFERENCES=1     # Prepare reference files and build indices
RUN_SAMPLES=0        # Process samples (QC, trimming, alignment)
RUN_TETOOLKIT=0      # Run TEToolkit analysis
RUN_RESULTS=0        # Process TEToolkit results
RUN_REPORTS=0        # Generate final reports

# Pipeline options
THREADS=${SLURM_CPUS_PER_TASK:-16}
MEMORY_GB=128
MEMORY_BYTES=$((MEMORY_GB * 1024 * 1024 * 1024))

# TEToolkit specific options
TETOOLKIT_MODE="multi"  # multi, uniq - how to handle multi-mapping reads
TETOOLKIT_FORMAT="BAM"  # BAM or SAM
TETOOLKIT_STRANDED="reverse"  # no, yes, reverse - match your library prep

# =============================================================================
# SYSTEM SETUP
# =============================================================================

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Check if step should be run
should_run_step() {
    local step=$1
    
    case "$step" in
        "modules")     [[ $RUN_MODULES -eq 1 ]] ;;
        "validate")    [[ $RUN_VALIDATE -eq 1 ]] ;;
        "setup")       [[ $RUN_SETUP -eq 1 ]] ;;
        "references")  [[ $RUN_REFERENCES -eq 1 ]] ;;
        "samples")     [[ $RUN_SAMPLES -eq 1 ]] ;;
        "tetoolkit")   [[ $RUN_TETOOLKIT -eq 1 ]] ;;
        "results")     [[ $RUN_RESULTS -eq 1 ]] ;;
        "reports")     [[ $RUN_REPORTS -eq 1 ]] ;;
        *)             return 1 ;;
    esac
}

# Parse sample information consistently
parse_sample() {
    local sample=$1
    local sample_name_var=$2
    local condition_var=$3
    local timepoint_var=$4
    local replicate_var=$5
    local read1_var=$6
    local read2_var=$7
    
    # Split by colon
    IFS=':' read -r name cond tp_or_r1 rep_or_r2 r1 r2 <<< "$sample"
    
    eval "$sample_name_var='$name'"
    eval "$condition_var='$cond'"
    
    if [[ -z "$r1" ]]; then
        # 4-field format: name:condition:read1:read2
        eval "$timepoint_var='NA'"
        eval "$replicate_var='1'"
        eval "$read1_var='$tp_or_r1'"
        eval "$read2_var='$rep_or_r2'"
    else
        # 6-field format: name:condition:timepoint:replicate:read1:read2
        eval "$timepoint_var='$tp_or_r1'"
        eval "$replicate_var='$rep_or_r2'"
        eval "$read1_var='$r1'"
        eval "$read2_var='$r2'"
    fi
}

# Extract TE-only counts from combined gene/TE matrix
extract_te_only_counts() {
    local input_file=$1
    local output_file=$2
    
    if [[ ! -f "$input_file" ]]; then
        warning "Input file not found: $input_file"
        return 1
    fi
    
    log "Extracting TE-only counts from $(basename "$input_file")..."
    
    # Check if file has header
    if head -1 "$input_file" | grep -q -E "(gene|TE|sample)" || head -1 "$input_file" | grep -q "^[0-9]" --invert-match; then
        # Has header - keep it
        head -1 "$input_file" > "$output_file"
        grep -v "FBgn" "$input_file" | tail -n +2 >> "$output_file"
    else
        # No header - just remove FBgn lines
        grep -v "FBgn" "$input_file" > "$output_file"
    fi
    
    # Report results
    local orig_lines=$(wc -l < "$input_file")
    local te_lines=$(wc -l < "$output_file")
    local removed_lines=$((orig_lines - te_lines))
    
    info "  Original: $orig_lines lines"
    info "  TE-only: $te_lines lines" 
    info "  Removed: $removed_lines FlyBase genes"
    info "  Output: $output_file"
    
    # Show first few TEs
    info "  First 3 TEs:"
    if head -1 "$output_file" | grep -q -E "(gene|TE|sample)" || head -1 "$output_file" | grep -q "^[0-9]" --invert-match; then
        tail -n +2 "$output_file" | head -3 | cut -f1 | sed 's/^/    /'
    else
        head -3 "$output_file" | cut -f1 | sed 's/^/    /'
    fi
    
    return 0
}

skip_step() {
    local step_name=$1
    info "Skipping $step_name (toggle set to 0)"
}

log() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1" >&2
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" >&2
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1" >&2
}

# Load required modules
load_modules() {
    log "Loading HPC modules..."
    if command -v module &> /dev/null; then
        # Purge modules (ignore errors from deactivate commands)
        module purge 2>/dev/null || true
        module load star/2.7.11b
        module load samtools/1.9
        module load R/4.2.2
        module load fastqc/0.12.1
        module load cutadapt/3.5
        module load bedtools/2.29.0
        module load MultiQC/1.9
        module load TEToolkit/2.0.3
        log "Modules loaded successfully"
    else
        warning "'module' command not found. Skipping module loading. Ensure all required tools are in your PATH."
    fi
}

# Verify tools
verify_tools() {
    log "Verifying tool availability..."
    
    local required_tools=("STAR" "samtools" "fastqc" "cutadapt" "multiqc")
    
    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            error "$tool not found. Make sure required modules are loaded."
        fi
    done
    
    # Check TEToolkit tools separately with more detailed error messages
    if ! command -v "TEcount" &> /dev/null; then
        error "TEcount not found. Check if TEToolkit/2.0.3 module is properly loaded. Try: module list"
    fi
    
    if ! command -v "TEtranscripts" &> /dev/null; then
        error "TEtranscripts not found. Check if TEToolkit/2.0.3 module is properly loaded."
    fi
    
    # Check R packages
    log "Checking R package dependencies..."
    
    # Set up personal R library path
    export R_LIBS_USER="$HOME/R/library"
    
    if ! Rscript -e "library(data.table)" >/dev/null 2>&1; then
        warning "R package 'data.table' not found."
        warning "This is needed for results processing (separate_counts.R script)."
        warning "Install manually with:"
        warning "  mkdir -p ~/R/library"
        warning "  Rscript -e \"install.packages('data.table', repos='https://cran.r-project.org', lib='~/R/library')\""
        warning "Pipeline will continue but results processing may fail."
        warning "You can install it later and re-run with RUN_RESULTS=1"
    else
        log "R package 'data.table' found"
    fi
    
    log "All tools verified"
}

# Validate input files
validate_inputs() {
    log "Validating input files..."
    
    for file in "$GENOME_FILE" "$GTF_FILE" "$REPEATMASKER_FILE" "$BERGMAN_TE_LIBRARY" "$PYTHON_SCRIPT" "$R_SCRIPT"; do
        if [[ ! -f "$file" ]]; then
            error "File not found: $file"
        fi
    done
    
    for sample in "${SAMPLES[@]}"; do
        local sample_name condition timepoint replicate read1 read2
        parse_sample "$sample" sample_name condition timepoint replicate read1 read2
        
        if [[ ! -f "$read1" || ! -f "$read2" ]]; then
            error "FASTQ files not found for sample: $sample_name ($read1, $read2)"
        fi
    done
    
    log "Input validation passed"
}

# Setup project directory
setup_project() {
    log "Setting up project directories..."
    mkdir -p logs
    mkdir -p data/trimmed
    mkdir -p data/aligned
    mkdir -p results/qc
    mkdir -p results/tetoolkit
    mkdir -p results/differential_expression
    mkdir -p scripts
    mkdir -p "$LOCAL_REF_DIR"
}

# Copy and prepare reference files
prepare_references() {
    log "Preparing reference files..."
    
    # Copy reference files locally
    log "Copying reference files..."
    cp "$GENOME_FILE" "$LOCAL_REF_DIR/genome.fasta"
    cp "$GTF_FILE" "$LOCAL_REF_DIR/genes.gtf"
    cp "$REPEATMASKER_FILE" "$LOCAL_REF_DIR/repeatmasker.out"
    cp "$BERGMAN_TE_LIBRARY" "$LOCAL_REF_DIR/te_library.fa"
    cp "$PYTHON_SCRIPT" "$LOCAL_REF_DIR/create_te_gtf.py"
    cp "$R_SCRIPT" "scripts/separate_counts.R"
    
    # Copy plotting script if it exists
    if [[ -f "$PLOTS_SCRIPT" ]]; then
        cp "$PLOTS_SCRIPT" "results/tetoolkit/"
        log "Plotting script copied to results/tetoolkit/"
    else
        warning "Plotting script not found: $PLOTS_SCRIPT"
    fi
    
    # Use the pre-built TE GTF file
    log "Using pre-built dm6_rmsk_TE.gtf file..."
    cp "$REFERENCE_DIR/dm6_rmsk_TE.gtf" "$LOCAL_REF_DIR/te_annotation_raw.gtf"
    
    # Check if chromosome names need fixing (remove chr prefix)
    local first_chr=$(head -1 "$LOCAL_REF_DIR/te_annotation_raw.gtf" | cut -f1)
    local genome_chr=$(grep "^>" "$LOCAL_REF_DIR/genome.fasta" | head -1 | cut -d' ' -f1 | sed 's/^>//')
    
    if [[ "$first_chr" == "chr"* ]] && ! grep -q "^>chr" "$LOCAL_REF_DIR/genome.fasta"; then
        log "Removing 'chr' prefix from TE GTF to match genome..."
        sed 's/^chr//' "$LOCAL_REF_DIR/te_annotation_raw.gtf" > "$LOCAL_REF_DIR/te_annotation.gtf"
    else
        cp "$LOCAL_REF_DIR/te_annotation_raw.gtf" "$LOCAL_REF_DIR/te_annotation.gtf"
    fi
    
    # Check the TE GTF
    local te_count=$(wc -l < "$LOCAL_REF_DIR/te_annotation.gtf")
    log "Using TE GTF with $te_count TE annotations"
    
    # Show first few lines for verification
    log "First few TE annotations:"
    head -3 "$LOCAL_REF_DIR/te_annotation.gtf" | while read line; do
        info "$line"
    done
    
    # Create combined reference for TEToolkit
    log "Creating combined genome + TE library..."
    cat "$LOCAL_REF_DIR/genome.fasta" "$LOCAL_REF_DIR/te_library.fa" > "$LOCAL_REF_DIR/genome_plus_tes.fa"
    
    # Build STAR index for combined reference
    log "Building STAR index for combined reference..."
    mkdir -p "$LOCAL_REF_DIR/star_index"
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$LOCAL_REF_DIR/star_index" \
         --genomeFastaFiles "$LOCAL_REF_DIR/genome_plus_tes.fa" \
         --sjdbGTFfile "$LOCAL_REF_DIR/genes.gtf" \
         --sjdbOverhang 100 \
         --limitGenomeGenerateRAM $MEMORY_BYTES \
         --genomeSAindexNbases 12
    
    log "Reference preparation completed"
}

# Process samples - single alignment approach for TEToolkit
process_samples() {
    log "Processing samples with TEToolkit approach..."
    
    for sample in "${SAMPLES[@]}"; do
        local sample_name condition timepoint replicate read1 read2
        parse_sample "$sample" sample_name condition timepoint replicate read1 read2
        
        log "Processing sample: $sample_name"
        
        # Quality control
        log "Running FastQC for $sample_name..."
        fastqc -o "results/qc" -t 2 "$read1" "$read2"
        
        # Adapter trimming
        log "Trimming adapters for $sample_name..."
        cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o "data/trimmed/${sample_name}_R1_trimmed.fastq.gz" \
            -p "data/trimmed/${sample_name}_R2_trimmed.fastq.gz" \
            --minimum-length 20 \
            --quality-cutoff 20 \
            -j $THREADS \
            "$read1" "$read2" > "logs/${sample_name}_cutadapt.log"
        
        # Alignment with STAR (optimized for TE analysis)
        log "Aligning $sample_name with STAR..."
        STAR --runThreadN $THREADS \
             --genomeDir "$LOCAL_REF_DIR/star_index" \
             --readFilesIn "data/trimmed/${sample_name}_R1_trimmed.fastq.gz" "data/trimmed/${sample_name}_R2_trimmed.fastq.gz" \
             --readFilesCommand zcat \
             --outFileNamePrefix "data/aligned/${sample_name}_" \
             --outSAMtype BAM Unsorted \
             --outSAMstrandField intronMotif \
             --outSAMattributes NH HI AS nM \
             --outFilterMultimapNmax 100 \
             --winAnchorMultimapNmax 100 \
             --outSAMmultNmax -1 \
             --outMultimapperOrder Random \
             --outSAMprimaryFlag OneBestScore \
             --outFilterMismatchNmax 10 \
             --alignIntronMax 1000000 \
             --outBAMsortingThreadN $THREADS
        
        # Sort BAM by queryname (read name) for TEcount compatibility
        log "Sorting BAM by queryname for TEcount..."
        samtools sort -n "data/aligned/${sample_name}_Aligned.out.bam" -o "data/aligned/${sample_name}.bam"
        rm "data/aligned/${sample_name}_Aligned.out.bam"  # Remove unsorted version to save space
        
        # Generate alignment statistics
        samtools flagstat "data/aligned/${sample_name}.bam" > "results/qc/${sample_name}_flagstat.txt"
        
        log "Sample $sample_name processing completed"
    done
}

# Run TEToolkit analysis - UPDATED FOR TREATED VS CONTROL
run_tetoolkit() {
    log "Running TEToolkit analysis for treated vs control..."
    
    # Collect all BAM files
    local bam_files=()
    local sample_names=()
    local conditions=()
    local treated_bams=()
    local control_bams=()
    
    log "Collecting BAM files..."
    for sample in "${SAMPLES[@]}"; do
        local sample_name condition timepoint replicate read1 read2
        parse_sample "$sample" sample_name condition timepoint replicate read1 read2
        
        local bam_path="data/aligned/${sample_name}.bam"
        if [[ -f "$bam_path" ]]; then
            bam_files+=("$bam_path")
            sample_names+=("$sample_name")
            conditions+=("$condition")
            
            # Separate into treatment groups
            if [[ "$condition" == "treated" ]]; then
                treated_bams+=("$bam_path")
                log "  Treated: $bam_path"
            elif [[ "$condition" == "control" ]]; then
                control_bams+=("$bam_path")
                log "  Control: $bam_path"
            fi
        else
            error "BAM file missing: $bam_path"
        fi
    done
    
    log "Total BAM files: ${#bam_files[@]}"
    log "Treated samples: ${#treated_bams[@]} (${treated_bams[*]##*/})"
    log "Control samples: ${#control_bams[@]} (${control_bams[*]##*/})"
    
    # Verify we have samples in both groups
    if [[ ${#treated_bams[@]} -eq 0 || ${#control_bams[@]} -eq 0 ]]; then
        error "Need samples in both treated and control groups for differential expression"
    fi
    
    # Create sample information file (simplified - just sample and condition)
    echo -e "sample\tcondition" > "results/tetoolkit/sample_info.txt"
    for sample in "${SAMPLES[@]}"; do
        local sample_name condition timepoint replicate read1 read2
        parse_sample "$sample" sample_name condition timepoint replicate read1 read2
        echo -e "${sample_name}\t${condition}" >> "results/tetoolkit/sample_info.txt"
    done
    
    log "Sample information saved to: results/tetoolkit/sample_info.txt"
    
    # Run TEcount for basic quantification
    log "Running TEcount for quantification..."
    
    # Create the TEcount command properly with all BAM files as separate arguments
    log "Building TEcount command with all BAM files..."
    
    # Start building the command
    local tecount_args=(
        "TEcount"
        "--format" "BAM"
        "--mode" "$TETOOLKIT_MODE"
        "--stranded" "$TETOOLKIT_STRANDED"
        "--project" "results/tetoolkit/tecount"
        "--GTF" "$LOCAL_REF_DIR/genes.gtf"
        "--TE" "$LOCAL_REF_DIR/te_annotation.gtf"
    )
    
    # Add each BAM file as a separate -b argument
    for bam in "${bam_files[@]}"; do
        tecount_args+=("-b" "$bam")
    done
    
    # Show the full command for debugging
    log "TEcount command: ${tecount_args[*]}"
    
    # Execute the command
    "${tecount_args[@]}"
    
    
    if [[ $? -eq 0 ]]; then
        log "TEcount completed successfully"
        
        # Verify all samples were processed
        if [[ -f "results/tetoolkit/tecount.cntTable" ]]; then
            local header_cols=$(head -1 "results/tetoolkit/tecount.cntTable" | tr '\t' '\n' | wc -l)
            local expected_cols=$((${#bam_files[@]} + 1))  # +1 for gene/TE name column
            log "Count matrix verification:"
            log "  Expected columns: $expected_cols (1 gene column + ${#bam_files[@]} samples)"
            log "  Actual columns: $header_cols"
            
            if [[ $header_cols -eq $expected_cols ]]; then
                log "  All samples processed correctly"
                
                # Show sample columns
                log "  Sample columns:"
                head -1 "results/tetoolkit/tecount.cntTable" | tr '\t' '\n' | tail -n +2 | nl | sed 's/^/    /'
            else
                warning "  Column count mismatch - check TEcount output"
            fi
        fi
    else
        error "TEcount failed"
    fi
    
    # Run TEtranscripts for differential expression
    log "Running TEtranscripts for treated vs control comparison..."
    
    log "Treatment BAMs: ${treated_bams[*]##*/}"
    log "Control BAMs: ${control_bams[*]##*/}"
    
    # Build TEtranscripts command with proper argument structure
    local tetranscripts_args=(
        "TEtranscripts"
        "--format" "BAM"
        "--mode" "$TETOOLKIT_MODE"
        "--stranded" "$TETOOLKIT_STRANDED"
        "--project" "results/tetoolkit/treated_vs_control"
        "--GTF" "$LOCAL_REF_DIR/genes.gtf"
        "--TE" "$LOCAL_REF_DIR/te_annotation.gtf"
    )
    
    # Add treatment BAMs
    tetranscripts_args+=("-t")
    for bam in "${treated_bams[@]}"; do
        tetranscripts_args+=("$bam")
    done
    
    # Add control BAMs
    tetranscripts_args+=("-c")
    for bam in "${control_bams[@]}"; do
        tetranscripts_args+=("$bam")
    done
    
    # Show the full command for debugging
    log "TEtranscripts command: ${tetranscripts_args[*]}"
    
    # Execute the command
    "${tetranscripts_args[@]}"
    
    if [[ $? -eq 0 ]]; then
        log "TEtranscripts completed successfully"
    else
        error "TEtranscripts failed"
    fi
    
    log "TEToolkit analysis completed"
}

# Process TEToolkit results - UPDATED WITH TE EXTRACTION
process_tetoolkit_results() {
    log "Processing TEToolkit results..."
    
    local results_dir="results/tetoolkit"
    
    # Process main TEtranscripts count table
    if [[ -f "$results_dir/treated_vs_control.cntTable" ]]; then
        log "Main count table found: treated_vs_control.cntTable"
        
        # Show dimensions
        local rows=$(tail -n +2 "$results_dir/treated_vs_control.cntTable" | wc -l)
        local cols=$(head -1 "$results_dir/treated_vs_control.cntTable" | tr '\t' '\n' | wc -l)
        log "  Dimensions: $rows features × $((cols-1)) samples"
        
        # Show sample columns
        log "  Sample columns:"
        head -1 "$results_dir/treated_vs_control.cntTable" | tr '\t' '\n' | tail -n +2 | nl -v0 | sed 's/^/    /'
        
        # Extract TE-only counts (only if there are FlyBase genes to remove)
        log "Checking if TE extraction is needed..."
        local fbgn_count=$(grep -c "FBgn" "$results_dir/treated_vs_control.cntTable" || echo "0")
        
        if [[ $fbgn_count -gt 0 ]]; then
            log "Found $fbgn_count FlyBase genes - creating TE-specific count matrix..."
            extract_te_only_counts "$results_dir/treated_vs_control.cntTable" "$results_dir/TE_only_treated_vs_control.cntTable"
        else
            log "No FlyBase genes found - treated_vs_control.cntTable appears to be TE-focused already"
            log "Skipping TE extraction (would produce identical file)"
        fi
        
    else
        warning "Main count table not found: $results_dir/treated_vs_control.cntTable"
    fi
    
    # Process basic TEcount output if it exists (but check if extraction is needed)
    if [[ -f "$results_dir/tecount.cntTable" ]]; then
        log "Basic TEcount output found: tecount.cntTable"
        local rows=$(tail -n +2 "$results_dir/tecount.cntTable" | wc -l)
        local fbgn_count=$(grep -c "FBgn" "$results_dir/tecount.cntTable" || echo "0")
        log "  Total features: $rows, FlyBase genes: $fbgn_count"
        
        # Only extract TEs if there are actually FlyBase genes to remove
        if [[ $fbgn_count -gt 0 ]]; then
            log "  Creating TE-only version (removing $fbgn_count genes)..."
            extract_te_only_counts "$results_dir/tecount.cntTable" "$results_dir/TE_only_tecount.cntTable"
        else
            log "  No FlyBase genes found - tecount.cntTable is already TE-focused"
        fi
    fi
    
    # Check for differential expression results
    if [[ -f "$results_dir/treated_vs_control_DESeq2.R" ]]; then
        log "DESeq2 R script generated: treated_vs_control_DESeq2.R"
        
        # Run the DESeq2 analysis
        log "Running DESeq2 differential expression analysis..."
        cd "$results_dir"
        
        if Rscript "treated_vs_control_DESeq2.R"; then
            log "DESeq2 analysis completed successfully"
            
            # Check for output files and extract TE-specific results
            local de_files=(treated_vs_control_DESeq2_*.txt)
            if [[ -f "${de_files[0]}" ]]; then
                log "Differential expression results generated:"
                for file in treated_vs_control_DESeq2_*.txt; do
                    if [[ -f "$file" ]]; then
                        local count=$(wc -l < "$file")
                        log "    $file ($((count-1)) features)"
                        
                        # Only extract TE-specific results if there are genes to remove
                        local file_fbgn_count=$(grep -c "FBgn" "$file" || echo "0")
                        if [[ $file_fbgn_count -gt 0 ]]; then
                            local te_file="${file%.txt}_TE_only.txt"
                            if extract_te_only_counts "$file" "$te_file"; then
                                local te_count=$(wc -l < "$te_file")
                                log "    $te_file ($((te_count-1)) TEs)"
                            fi
                        else
                            log "    No genes to remove - file appears TE-focused already"
                        fi
                    fi
                done
            fi
            
        else
            warning "DESeq2 analysis failed - check R dependencies and log files"
        fi
        
        cd ../..
    else
        warning "DESeq2 R script not found: $results_dir/treated_vs_control_DESeq2.R"
    fi
    
    # Copy plotting script to results directory for manual use
    if [[ -f "$PLOTS_SCRIPT" ]]; then
        cp "$PLOTS_SCRIPT" "$results_dir/"
        log "Plotting script copied to $results_dir/ for manual use"
        log "To generate plots manually, run: cd $results_dir && Rscript create_deseq2_plots.R"
    fi
    
    log "Results processing completed"
}

# Generate comprehensive reports - UPDATED WITH TE EXTRACTION INFO
generate_reports() {
    log "Generating reports..."
    
    # Run MultiQC
    if command -v multiqc &> /dev/null; then
        multiqc . -o results/qc/ -n multiqc_report.html 2>/dev/null || warning "MultiQC failed - continuing anyway"
    fi
    
    # Create analysis summary
    cat > results/ANALYSIS_SUMMARY.md << EOF
# Drosophila TE Analysis: Treated vs Control with TEToolkit

## Project: $PROJECT_NAME
**Analysis completed:** $(date)
**Comparison:** Treated (E samples) vs Control (C samples)
**TEToolkit mode:** $TETOOLKIT_MODE
**Strand specificity:** $TETOOLKIT_STRANDED

## Sample Groups
**Treated samples (3):** 1E_rep1, 2E_rep2, 3E_rep3
**Control samples (3):** 4C_rep1, 5C_rep2, 6C_rep3

## Key Output Files

### Main Results
- **Count matrix (all)**: \`results/tetoolkit/treated_vs_control.cntTable\`
- **Count matrix (TE-only)**: \`results/tetoolkit/TE_only_treated_vs_control.cntTable\`
- **Differential expression (all)**: \`results/tetoolkit/treated_vs_control_DESeq2_*.txt\`
- **Differential expression (TE-only)**: \`results/tetoolkit/treated_vs_control_DESeq2_*_TE_only.txt\`
- **Sample information**: \`results/tetoolkit/sample_info.txt\`

### Quality Control
- **MultiQC report**: \`results/qc/multiqc_report.html\`
- **FastQC reports**: \`results/qc/*_fastqc.html\`
- **Alignment stats**: \`results/qc/*_flagstat.txt\`

### Visualization (Manual)
- **Plot script**: \`results/tetoolkit/create_deseq2_plots.R\`
- **To generate plots**: \`cd results/tetoolkit && Rscript create_deseq2_plots.R\`

### Reference Files
- **Combined genome**: \`$LOCAL_REF_DIR/genome_plus_tes.fa\`
- **Gene GTF**: \`$LOCAL_REF_DIR/genes.gtf\`
- **TE GTF**: \`$LOCAL_REF_DIR/te_annotation.gtf\`
- **STAR index**: \`$LOCAL_REF_DIR/star_index/\`

## Analysis Approach

This pipeline uses TEToolkit's integrated approach:
1. **Single alignment** to combined genome + TE library
2. **TEcount** for simultaneous gene and TE quantification
3. **TEtranscripts** for differential expression analysis
4. **Smart TE extraction** removes FlyBase genes when present
5. **Proper handling** of multi-mapping reads \(mode: $TETOOLKIT_MODE\)

## Key Advantages of TEToolkit Approach
- Unified quantification of genes and TEs
- Specialized handling of repetitive elements
- Built-in statistical analysis
- Reduced analysis complexity

## Sample Information
EOF

    # Add sample information to summary
    echo "| Sample | Condition | Status |" >> results/ANALYSIS_SUMMARY.md
    echo "|--------|-----------|--------|" >> results/ANALYSIS_SUMMARY.md
    
    for sample in "${SAMPLES[@]}"; do
        local sample_name condition timepoint replicate read1 read2
        parse_sample "$sample" sample_name condition timepoint replicate read1 read2
        
        if [[ -f "data/aligned/${sample_name}.bam" ]]; then
            status="Completed"
        else
            status="Failed"
        fi
        echo "| $sample_name | $condition | $status |" >> results/ANALYSIS_SUMMARY.md
    done

    cat >> results/ANALYSIS_SUMMARY.md << EOF

## Next Steps
1. **Check results**: 
   - All features: \`treated_vs_control.cntTable\`
   - TE-only: \`TE_only_treated_vs_control.cntTable\` (if genes were present)
2. **Differential expression**: 
   - All features: \`treated_vs_control_DESeq2_*.txt\`
   - TE-only: \`treated_vs_control_DESeq2_*_TE_only.txt\` (if genes were present)
3. **Generate plots manually**: 
   - \`cd results/tetoolkit && Rscript create_deseq2_plots.R\`
4. **Quality control**: Review MultiQC report and alignment statistics

## File Descriptions
- **treated_vs_control.cntTable**: Combined gene and TE count matrix from TEtranscripts
- **TE_only_treated_vs_control.cntTable**: TE-only count matrix (FlyBase genes removed)
- **treated_vs_control_DESeq2_gene.txt**: Differential gene expression results
- **treated_vs_control_DESeq2_TE.txt**: Differential TE expression results  
- **treated_vs_control_DESeq2_*_TE_only.txt**: TE-only differential expression results
- **treated_vs_control_DESeq2.R**: R script used for analysis (can be rerun/modified)
- **create_deseq2_plots.R**: R script for generating visualizations (run manually)

## TE Analysis Notes
- **Smart extraction**: Only creates TE-only files when FlyBase genes are actually present
- **FlyBase gene removal**: All features starting with "FBgn" are considered genes
- **TE identification**: Remaining features are transposable elements from the Bergman library
- **Focus on TEs**: Use TE_only_* files for transposon-specific analyses when available
EOF

    log "Analysis summary created: results/ANALYSIS_SUMMARY.md"
}

# Main execution
main() {
    log "Starting Drosophila TE Analysis with TEToolkit"
    info "Project: $PROJECT_NAME"
    info "Samples: ${#SAMPLES[@]}"
    info "TEToolkit mode: $TETOOLKIT_MODE"
    
    # Show which steps will run
    log "Step execution plan:"
    info "  Load modules: $([ $RUN_MODULES -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Validate inputs: $([ $RUN_VALIDATE -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Setup directories: $([ $RUN_SETUP -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Prepare references: $([ $RUN_REFERENCES -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Process samples: $([ $RUN_SAMPLES -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Run TEToolkit: $([ $RUN_TETOOLKIT -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Process results: $([ $RUN_RESULTS -eq 1 ] && echo "RUN" || echo "SKIP")"
    info "  Generate reports: $([ $RUN_REPORTS -eq 1 ] && echo "RUN" || echo "SKIP")"

    # Trap for error handling
    trap 'warning "Pipeline exited early. Check outputs before proceeding to next step."' ERR
    
    # Step 1: Load modules
    if should_run_step "modules"; then
        load_modules
        verify_tools
    else
        skip_step "module loading and tool verification"
    fi
    
    # Step 2: Validate inputs
    if should_run_step "validate"; then
        validate_inputs
    else
        skip_step "input validation"
    fi
    
    # Step 3: Setup project
    if should_run_step "setup"; then
        setup_project
    else
        skip_step "project setup"
    fi
    
    # Step 4: Prepare references
    if should_run_step "references"; then
        prepare_references
    else
        skip_step "reference preparation"
    fi
    
    # Step 5: Process samples
    if should_run_step "samples"; then
        process_samples
    else
        skip_step "sample processing"
    fi
    
    # Step 6: Run TEToolkit
    if should_run_step "tetoolkit"; then
        run_tetoolkit
    else
        skip_step "TEToolkit analysis"
    fi
    
    # Step 7: Process results
    if should_run_step "results"; then
        process_tetoolkit_results
    else
        skip_step "results processing"
    fi
    
    # Step 8: Generate reports
    if should_run_step "reports"; then
        generate_reports
    else
        skip_step "report generation"
    fi

    log "Pipeline completed successfully!"
    log "Results summary: results/ANALYSIS_SUMMARY.md"
    log "TEToolkit outputs: results/tetoolkit/"
    log "To generate plots: cd results/tetoolkit && Rscript create_deseq2_plots.R"
}

# Execute main function
if [[ "${1:-}" == "help" || "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    echo "Drosophila TE Analysis Pipeline with TEToolkit"
    echo "=============================================="
    echo ""
    echo "Usage: $0"
    echo ""
    echo "Step Control:"
    echo "Edit the step toggles at the top of this script to control which steps run:"
    echo "  RUN_MODULES=1      # Load modules and verify tools"
    echo "  RUN_VALIDATE=1     # Validate input files"
    echo "  RUN_SETUP=1        # Setup project directories"
    echo "  RUN_REFERENCES=1   # Prepare reference files and build indices"
    echo "  RUN_SAMPLES=1      # Process samples (QC, trimming, alignment)"
    echo "  RUN_TETOOLKIT=1    # Run TEToolkit analysis"
    echo "  RUN_RESULTS=1      # Process TEToolkit results"
    echo "  RUN_REPORTS=1      # Generate final reports"
    echo ""
    echo "Set any toggle to 0 to skip that step, 1 to run it."
    echo ""
    echo "Recommended workflow:"
    echo "1. Run modules/validate/setup first"
    echo "2. Run references (takes longest - building STAR index)"
    echo "3. Run samples (alignment step)"
    echo "4. Run tetoolkit (quantification)"
    echo "5. Run results/reports (final processing)"
    echo "6. Generate plots manually: cd results/tetoolkit && Rscript create_deseq2_plots.R"
    echo ""
    echo "NEW: Smart TE extraction automatically removes FlyBase genes when present"
    echo "Creates TE-only versions of count matrices and differential expression results"
    exit 0
fi

main "$@"
