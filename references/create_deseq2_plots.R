#!/usr/bin/env Rscript

# Standalone DESeq2 Visualization Script for TEToolkit Results
# Auto-detects experimental conditions and creates comprehensive plots
# Now with automatic dependency installation

cat("=== DESeq2 Visualization Script ===\n")
cat("Setting up dependencies...\n")

# Function to install and load packages
install_and_load <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            cat("Installing missing package:", pkg, "\n")
            
            # Check if it's a Bioconductor package
            if (pkg %in% c("DESeq2", "BiocManager")) {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager", repos = "https://cran.r-project.org")
                }
                BiocManager::install(pkg, ask = FALSE, update = FALSE)
            } else {
                # CRAN package
                install.packages(pkg, repos = "https://cran.r-project.org")
            }
        }
        
        # Load the package
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        cat("✓", pkg, "loaded\n")
    }
}

# Required packages
required_packages <- c("DESeq2", "ggplot2", "pheatmap")

# Install and load packages
tryCatch({
    install_and_load(required_packages)
    cat("All dependencies loaded successfully!\n\n")
}, error = function(e) {
    cat("ERROR installing packages:", e$message, "\n")
    cat("Please try installing manually:\n")
    cat("install.packages(c('ggplot2', 'pheatmap'))\n")
    cat("BiocManager::install('DESeq2')\n")
    stop("Dependency installation failed")
})

cat("Auto-detecting TEtranscripts results...\n\n")

# Function to find DESeq2 R script
find_deseq2_script <- function() {
    # Look for TEtranscripts generated R script
    possible_files <- c(
        "treated_vs_control_DESeq2.R",
        list.files(pattern = "*_DESeq2\\.R$")
    )
    
    for (file in possible_files) {
        if (file.exists(file)) {
            cat("Found DESeq2 script:", file, "\n")
            return(file)
        }
    }
    
    stop("ERROR: No DESeq2 R script found. Make sure TEtranscripts completed successfully.")
}

# Function to find differential expression results
find_de_results <- function() {
    # Look for differential expression results
    possible_files <- c(
        "treated_vs_control_sigdiff_gene_TE.txt",
        list.files(pattern = "*sigdiff_gene_TE\\.txt$")
    )
    
    found_files <- c()
    for (file in possible_files) {
        if (file.exists(file)) {
            found_files <- c(found_files, file)
        }
    }
    
    return(found_files)
}

# Load DESeq2 results
deseq2_script <- find_deseq2_script()
cat("Loading DESeq2 results from:", deseq2_script, "\n")

# Function to fix paths in DESeq2 script
fix_deseq2_paths <- function(script_file) {
    # Read the script content
    script_content <- readLines(script_file)
    
    # Look for common problematic paths and fix them
    script_content <- gsub("results/tetoolkit/", "", script_content)
    script_content <- gsub("'results/tetoolkit/", "'", script_content)
    script_content <- gsub('"results/tetoolkit/', '"', script_content)
    
    # Create a temporary fixed script
    temp_script <- paste0(tempfile(), ".R")
    writeLines(script_content, temp_script)
    
    return(temp_script)
}

tryCatch({
    # Check if we need to fix paths
    script_content <- readLines(deseq2_script)
    needs_path_fix <- any(grepl("results/tetoolkit/", script_content))
    
    if (needs_path_fix) {
        cat("Detected path issues, fixing paths in DESeq2 script...\n")
        fixed_script <- fix_deseq2_paths(deseq2_script)
        source(fixed_script)
        # Clean up temp file
        unlink(fixed_script)
    } else {
        source(deseq2_script)
    }
    
    cat("✓ DESeq2 data loaded successfully\n")
}, error = function(e) {
    # Try alternative approach: change working directory temporarily
    cat("Path fix didn't work, trying directory change approach...\n")
    
    # Look for the actual data files in subdirectories
    possible_dirs <- c("results/tetoolkit", "tetoolkit", ".")
    found_dir <- NULL
    
    for (dir in possible_dirs) {
        if (dir.exists(dir)) {
            # Check if the expected files exist in this directory
            expected_files <- c("treated_vs_control.cntTable", 
                              list.files(pattern = "\\.cntTable$"))
            
            for (file in expected_files) {
                if (file.exists(file.path(dir, file))) {
                    found_dir <- dir
                    break
                }
            }
            if (!is.null(found_dir)) break
        }
    }
    
    if (!is.null(found_dir) && found_dir != ".") {
        cat("Found data files in:", found_dir, "\n")
        cat("Changing working directory temporarily...\n")
        
        old_wd <- getwd()
        setwd(found_dir)
        
        tryCatch({
            source(file.path(old_wd, deseq2_script))
            setwd(old_wd)  # Change back
            cat("✓ DESeq2 data loaded successfully\n")
        }, error = function(e2) {
            setwd(old_wd)  # Make sure to change back
            stop("ERROR loading DESeq2 script even after directory change: ", e2$message)
        })
    } else {
        # List available files to help troubleshoot
        cat("Available files in current directory:\n")
        files <- list.files(pattern = "\\.(cntTable|txt|R)$", recursive = TRUE)
        if (length(files) > 0) {
            for (f in files) cat("  -", f, "\n")
        } else {
            cat("  No relevant files found\n")
        }
        
        stop("ERROR loading DESeq2 script: ", e$message, 
             "\nCould not find the required data files. Please check that TEtranscripts completed successfully.")
    }
})

# Verify required objects exist and find the DESeqDataSet object
cds <- NULL
possible_names <- c("cds", "dds", "deseq_data", "ds", "dataset")

# Check for DESeqDataSet objects in the environment
all_objects <- ls()
deseq_objects <- c()

for (obj_name in all_objects) {
    obj <- get(obj_name)
    if (inherits(obj, "DESeqDataSet")) {
        deseq_objects <- c(deseq_objects, obj_name)
    }
}

if (length(deseq_objects) == 0) {
    cat("Available objects in environment:\n")
    for (obj in all_objects) {
        cat("  -", obj, ":", class(get(obj))[1], "\n")
    }
    stop("ERROR: No DESeqDataSet object found in the environment")
} else if (length(deseq_objects) == 1) {
    cds_name <- deseq_objects[1]
    cds <- get(cds_name)
    cat("Found DESeqDataSet object:", cds_name, "\n")
    if (cds_name != "cds") {
        # Create an alias for consistency
        assign("cds", cds)
    }
} else {
    cat("Multiple DESeqDataSet objects found:\n")
    for (obj in deseq_objects) {
        cat("  -", obj, "\n")
    }
    # Use the first one found
    cds_name <- deseq_objects[1]
    cds <- get(cds_name)
    cat("Using:", cds_name, "\n")
    if (cds_name != "cds") {
        assign("cds", cds)
    }
}

# Also check for results object
res_obj <- NULL
if (exists("res")) {
    res_obj <- res
    cat("Found results object: res\n")
} else if (exists("results")) {
    res_obj <- results
    cat("Found results object: results\n")
} else {
    cat("No pre-computed results object found, will generate from DESeqDataSet\n")
}

# Get sample information and auto-detect conditions
sample_info <- as.data.frame(colData(cds))

# Handle the specific column naming from TEtranscripts
if ("groups" %in% colnames(sample_info)) {
    sample_info$condition <- sample_info$groups
} else if ("condition" %in% colnames(sample_info)) {
    # Already has condition column
} else {
    # Use the first factor column as condition
    factor_cols <- sapply(sample_info, is.factor)
    if (any(factor_cols)) {
        condition_col <- names(factor_cols)[factor_cols][1]
        sample_info$condition <- sample_info[[condition_col]]
        cat("Using column", condition_col, "as condition\n")
    } else {
        stop("Could not identify condition column in sample data")
    }
}

conditions <- unique(sample_info$condition)
n_samples <- nrow(sample_info)

cat("Sample Information Detected:\n")
cat("  Total samples:", n_samples, "\n")
cat("  Conditions:", paste(conditions, collapse = ", "), "\n")

# Print sample table
cat("  Sample details:\n")
for (i in 1:nrow(sample_info)) {
    cat("    ", rownames(sample_info)[i], ":", sample_info$condition[i], "\n")
}
cat("\n")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots/combined", showWarnings = FALSE)
dir.create("plots/te_only", showWarnings = FALSE)

# Function to check for TE-only count table
find_te_only_data <- function() {
    possible_files <- c(
        "TE_only_treated_vs_control.cntTable",
        "results/tetoolkit/TE_only_treated_vs_control.cntTable",
        list.files(pattern = "TE_only.*\\.cntTable$", recursive = TRUE)
    )
    
    for (file in possible_files) {
        if (file.exists(file)) {
            cat("Found TE-only count table:", file, "\n")
            return(file)
        }
    }
    
    return(NULL)
}

# Load TE-only data if available
te_only_file <- find_te_only_data()
te_cds <- NULL

if (!is.null(te_only_file)) {
    cat("Loading TE-only data for separate analysis...\n")
    tryCatch({
        # Load TE-only count data
        te_data <- read.table(te_only_file, header = TRUE, row.names = 1)
        
        # Apply same filtering as original script
        min_read <- 1
        te_data <- te_data[apply(te_data, 1, function(x){max(x)}) > min_read, ]
        
        # Create sample info (should match the original)
        te_groups <- factor(c(rep("TGroup", 3), rep("CGroup", 3)))
        te_sampleInfo <- data.frame(groups = te_groups, row.names = colnames(te_data))
        
        # Create DESeq2 dataset for TEs
        te_cds <- DESeqDataSetFromMatrix(countData = te_data, 
                                        colData = te_sampleInfo, 
                                        design = ~ groups)
        te_cds$condition <- relevel(te_cds$groups, "CGroup")
        te_cds <- DESeq(te_cds)
        te_res <- results(te_cds, independentFiltering = FALSE)
        
        cat("✓ TE-only DESeq2 analysis completed\n")
        cat("  TE features:", nrow(te_data), "\n")
        
    }, error = function(e) {
        cat("✗ Error loading TE-only data:", e$message, "\n")
        te_cds <- NULL
    })
}

# Find differential expression files
de_files <- find_de_results()
has_combined <- any(grepl("sigdiff_gene_TE\\.txt$", de_files))
has_te_cds <- !is.null(te_cds)

cat("Analysis options available:\n")
cat("  Combined (genes + TEs):", ifelse(has_combined, "✓", "✗"), "\n")
cat("  TE-only analysis:", ifelse(has_te_cds, "✓", "✗"), "\n")

if (length(de_files) > 0) {
    cat("  Differential expression files found:\n")
    for (file in de_files) {
        cat("    -", file, "\n")
    }
}
cat("\n")

# ============================================================================
# COMBINED PLOTS (Genes + TEs)
# ============================================================================

if (has_combined) {
    cat("Creating combined plots (genes + TEs)...\n")
    
    # 1. PCA Plot
    cat("  1. PCA plot...")
    tryCatch({
        rld <- rlog(cds, blind = FALSE)
        
        pca_plot <- plotPCA(rld, intgroup = "condition") + 
            ggtitle("PCA: All Features (Genes + TEs)") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        ggsave("plots/combined/pca_plot.pdf", pca_plot, width = 8, height = 6)
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 2. MA Plot
    cat("  2. MA plot...")
    tryCatch({
        # Use the results object we found or generate one
        if (!is.null(res_obj)) {
            plot_res <- res_obj
        } else if (exists("res")) {
            plot_res <- res
        } else {
            plot_res <- results(cds)
        }
        
        pdf("plots/combined/ma_plot.pdf", width = 8, height = 6)
        plotMA(plot_res, main = "MA Plot: All Features (Genes + TEs)", ylim = c(-4, 4))
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 3. Volcano Plot
    cat("  3. Volcano plot...")
    tryCatch({
        combined_file <- de_files[grepl("sigdiff_gene_TE\\.txt$", de_files)][1]
        de_data <- read.table(combined_file, header = TRUE, row.names = 1)
        de_data <- de_data[!is.na(de_data$padj) & is.finite(de_data$padj), ]
        de_data$log10padj <- -log10(de_data$padj)
        
        # Count significant features
        sig_up <- sum(de_data$padj < 0.05 & de_data$log2FoldChange > 1, na.rm = TRUE)
        sig_down <- sum(de_data$padj < 0.05 & de_data$log2FoldChange < -1, na.rm = TRUE)
        
        pdf("plots/combined/volcano_plot.pdf", width = 10, height = 8)
        plot(de_data$log2FoldChange, de_data$log10padj,
             xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value",
             main = paste("Volcano Plot: All Features\nUp-regulated:", sig_up, "| Down-regulated:", sig_down),
             pch = 20, col = "gray60", cex = 0.8)
        
        # Highlight significant points
        sig_idx <- de_data$padj < 0.05 & abs(de_data$log2FoldChange) > 1
        points(de_data$log2FoldChange[sig_idx], de_data$log10padj[sig_idx], 
               col = "red", pch = 20, cex = 0.8)
        
        # Add threshold lines
        abline(h = -log10(0.05), lty = 2, col = "blue", alpha = 0.7)
        abline(v = c(-1, 1), lty = 2, col = "blue", alpha = 0.7)
        
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 4. Heatmap
    cat("  4. Heatmap...")
    tryCatch({
        if (!exists("rld")) rld <- rlog(cds, blind = FALSE)
        
        # Get top 50 variable features
        rv <- rowVars(assay(rld))
        select <- order(rv, decreasing = TRUE)[1:min(50, length(rv))]
        
        pdf("plots/combined/heatmap.pdf", width = 10, height = 12)
        pheatmap(assay(rld)[select, ], 
                show_rownames = FALSE,
                main = "Heatmap: Top 50 Variable Features (Genes + TEs)",
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE)
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 5. Count Plots
    cat("  5. Count plots...")
    tryCatch({
        combined_file <- de_files[grepl("sigdiff_gene_TE\\.txt$", de_files)][1]
        de_data <- read.table(combined_file, header = TRUE, row.names = 1)
        de_data <- de_data[!is.na(de_data$padj), ]
        
        # Get top 6 significant features
        top_features <- head(rownames(de_data[order(de_data$padj), ]), 6)
        valid_features <- top_features[top_features %in% rownames(cds)]
        
        if (length(valid_features) > 0) {
            pdf("plots/combined/count_plots.pdf", width = 12, height = 8)
            par(mfrow = c(2, 3))
            
            for (feature in valid_features[1:min(6, length(valid_features))]) {
                plotCounts(cds, gene = feature, intgroup = "condition", 
                          main = paste0(feature, "\n(padj = ", 
                                      formatC(de_data[feature, "padj"], format = "e", digits = 2), ")"))
            }
            dev.off()
        }
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
}

# ============================================================================
# TE-ONLY PLOTS (from separate TE-only analysis)
# ============================================================================

if (has_te_cds) {
    cat("\nCreating TE-only plots from separate analysis...\n")
    
    # 1. TE-only PCA Plot
    cat("  1. TE PCA plot...")
    tryCatch({
        te_rld <- rlog(te_cds, blind = FALSE)
        
        te_pca_plot <- plotPCA(te_rld, intgroup = "condition") + 
            ggtitle("PCA: Transposable Elements Only") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        
        ggsave("plots/te_only/pca_plot_TE_only.pdf", te_pca_plot, width = 8, height = 6)
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 2. TE-only MA Plot
    cat("  2. TE MA plot...")
    tryCatch({
        pdf("plots/te_only/ma_plot_TE_only.pdf", width = 8, height = 6)
        plotMA(te_res, main = "MA Plot: Transposable Elements Only", ylim = c(-4, 4))
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 3. TE-only Volcano Plot
    cat("  3. TE volcano plot...")
    tryCatch({
        te_data_plot <- as.data.frame(te_res)
        te_data_plot <- te_data_plot[!is.na(te_data_plot$padj) & is.finite(te_data_plot$padj), ]
        te_data_plot$log10padj <- -log10(te_data_plot$padj)
        
        # Count significant TEs
        te_sig_up <- sum(te_data_plot$padj < 0.05 & te_data_plot$log2FoldChange > 1, na.rm = TRUE)
        te_sig_down <- sum(te_data_plot$padj < 0.05 & te_data_plot$log2FoldChange < -1, na.rm = TRUE)
        
        pdf("plots/te_only/volcano_plot_TE_only.pdf", width = 10, height = 8)
        plot(te_data_plot$log2FoldChange, te_data_plot$log10padj,
             xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value",
             main = paste("Volcano Plot: Transposable Elements Only\nUp-regulated:", te_sig_up, "| Down-regulated:", te_sig_down),
             pch = 20, col = "gray60", cex = 0.8)
        
        # Highlight significant TEs
        te_sig_idx <- te_data_plot$padj < 0.05 & abs(te_data_plot$log2FoldChange) > 1
        points(te_data_plot$log2FoldChange[te_sig_idx], te_data_plot$log10padj[te_sig_idx], 
               col = "darkred", pch = 20, cex = 1)
        
        # Add threshold lines
        abline(h = -log10(0.05), lty = 2, col = "blue", alpha = 0.7)
        abline(v = c(-1, 1), lty = 2, col = "blue", alpha = 0.7)
        
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 4. TE-only Heatmap
    cat("  4. TE heatmap...")
    tryCatch({
        if (!exists("te_rld")) te_rld <- rlog(te_cds, blind = FALSE)
        
        # Get top 50 variable TEs
        te_rv <- rowVars(assay(te_rld))
        te_select <- order(te_rv, decreasing = TRUE)[1:min(50, length(te_rv))]
        
        pdf("plots/te_only/heatmap_TE_only.pdf", width = 10, height = 12)
        pheatmap(assay(te_rld)[te_select, ], 
                show_rownames = FALSE,
                main = "Heatmap: Top 50 Variable Transposable Elements",
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE)
        dev.off()
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
    
    # 5. TE-only Count Plots
    cat("  5. TE count plots...")
    tryCatch({
        te_data_sig <- as.data.frame(te_res)
        te_data_sig <- te_data_sig[!is.na(te_data_sig$padj), ]
        
        # Get top 6 significant TEs
        top_tes <- head(rownames(te_data_sig[order(te_data_sig$padj), ]), 6)
        valid_tes <- top_tes[top_tes %in% rownames(te_cds)]
        
        if (length(valid_tes) > 0) {
            pdf("plots/te_only/count_plots_TE_only.pdf", width = 12, height = 8)
            par(mfrow = c(2, 3))
            
            for (te in valid_tes[1:min(6, length(valid_tes))]) {
                plotCounts(te_cds, gene = te, intgroup = "condition", 
                          main = paste0(te, "\n(padj = ", 
                                      formatC(te_data_sig[te, "padj"], format = "e", digits = 2), ")"))
            }
            dev.off()
        }
        cat(" ✓\n")
    }, error = function(e) { cat(" ✗ (", e$message, ")\n") })
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== VISUALIZATION SUMMARY ===\n")

# Count generated files
combined_plots <- list.files("plots/combined", pattern = "\\.pdf$", full.names = FALSE)
te_plots <- list.files("plots/te_only", pattern = "\\.pdf$", full.names = FALSE)

cat("Generated plots:\n")
cat("  Combined (genes + TEs):", length(combined_plots), "files\n")
if (length(combined_plots) > 0) {
    for (plot in combined_plots) {
        cat("    - plots/combined/", plot, "\n", sep = "")
    }
}

cat("  TE-only:", length(te_plots), "files\n")
if (length(te_plots) > 0) {
    for (plot in te_plots) {
        cat("    - plots/te_only/", plot, "\n", sep = "")
    }
}

total_plots <- length(combined_plots) + length(te_plots)
cat("\nTotal plots created:", total_plots, "\n")

if (total_plots > 0) {
    cat("\n✓ Visualization completed successfully!\n")
    cat("Check the plots/ directory for all generated visualizations.\n")
} else {
    cat("\n⚠ No plots were created. Check that DESeq2 analysis completed successfully.\n")
}
