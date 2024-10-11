
#### load libraries & utility function 
library("Seurat")
# library("data.table")
library("ggplot2")
library("dplyr")

# source utility functions
# source("workflow/scripts/utils.R")
### Utility Functions
library("data.table")

# save Seurat object results
save_seurat_object <- function (seurat_obj, result_dir){
    
    # make result directory if not exist
    if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

    # save seurat object
    saveRDS(seurat_obj, file=file.path(result_dir, "object.rds"))
    # save metadata
    fwrite(as.data.frame(seurat_obj[[]]), file = file.path(result_dir, "metadata.csv"), row.names=TRUE)
    # save stats
    stats <- paste0("cells: ",ncol(seurat_obj),"\nfeatures: ",nrow(seurat_obj))
    write(stats, file=file.path(result_dir, "stats.txt"))
}


# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5){
    # make result directory if not exist
    if (!dir.exists(results_path)){
        dir.create(results_path, recursive = TRUE)
    }
    
    # failsafe for large plots
    width <- min(100, width)
    height <- min(100, height)
    
    for (format in c('png')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = 300,
            width = width,
            height = height,
          limitsize = FALSE,
            units="in"
        )
    }
}

# helper function to check if all values in a group are the same
all_equal <- function(x) {
  if(length(unique(x)) == 1) unique(x) else NA
}

# inputs
filtered_object_path <- snakemake@input[["filtered_object"]]

# outputs
pseudobulk_counts_path <- snakemake@output[["pseudobulk_counts"]]
metadata_path <- snakemake@output[["metadata"]]
cell_count_plot_path <- snakemake@output[["cell_count_plot"]]

# parameters
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]

pseudobulk_by <- snakemake@config[["pseudobulk"]][["by"]]
pseudobulk_method <- match.fun(snakemake@config[["pseudobulk"]][["method"]]) # can be "sum", "mean", or "median"
pseudobulk_th <- snakemake@config[["pseudobulk"]][["cell_count_th"]]

### load filtered data
seurat_object <- readRDS(file = file.path(filtered_object_path))

# get metadata
metadata <- as.data.frame(seurat_object[[]])
metadata$ID <- rownames(metadata)

# generate aggregated metadata sheet
metadata_aggregated <- metadata %>%
  group_by(across(all_of(pseudobulk_by))) %>%
  summarise(across(everything(), all_equal), cell_count = n(), .groups = "drop") %>% # retain all columns that have the same value within a group
  select(-where(~ any(is.na(.)))) # Remove columns that have NAs
# format metadata
metadata_aggregated <- as.data.frame(metadata_aggregated)
rownames(metadata_aggregated) <- apply(metadata_aggregated[, pseudobulk_by], 1, function(x) paste(x, collapse = "_"))
# filter by cell_count_th
metadata_aggregated <- metadata_aggregated[metadata_aggregated$cell_count>=pseudobulk_th, ]


# pseudobulk per modality (RNA, ...)
for (modality in c("RNA", ab_flag, crispr_flag, custom_flag)){
    if (modality==""){
        next
    }

    # extract data
    tmp_data <- as.data.frame(t(as.data.frame(GetAssayData(object = seurat_object, slot = "counts", assay = modality))))
    features <- colnames(tmp_data)
    tmp_data$ID <- rownames(tmp_data)

    # join with metadata
    tmp_data <- inner_join(tmp_data, metadata[,c("ID",pseudobulk_by)], by = "ID")

    # pseudobulk by method
    tmp_pseudobulk <- tmp_data %>%
      group_by(across(all_of(pseudobulk_by))) %>%
      summarise(across(all_of(features), pseudobulk_method, .names = "{.col}"), .groups = "drop")

    # convert to integers
    numeric_cols <- sapply(tmp_pseudobulk, is.numeric)
    tmp_pseudobulk[, numeric_cols] <- lapply(tmp_pseudobulk[, numeric_cols], function(x) as.integer(round(x)))

    # reformat df
    tmp_pseudobulk <- as.data.frame(tmp_pseudobulk)
    rownames(tmp_pseudobulk) <- apply(tmp_pseudobulk[, pseudobulk_by], 1, function(x) paste(x, collapse = "_"))
    tmp_pseudobulk[,pseudobulk_by] <- NULL
    tmp_pseudobulk <- as.data.frame(t(tmp_pseudobulk))

    # filter by cell_count_th
    tmp_pseudobulk <- tmp_pseudobulk[,rownames(metadata_aggregated)]

    # add total counts for modality to metadata
    metadata_aggregated[[paste0("total_counts_",modality)]] <- colSums(tmp_pseudobulk)[rownames(metadata_aggregated)]

    # save pseudobulked data frame
    fwrite(as.data.frame(tmp_pseudobulk), file=file.path(dirname(pseudobulk_counts_path),paste0(modality,".csv")), row.names=TRUE)
}

# save metadata
fwrite(as.data.frame(metadata_aggregated), file=file.path(metadata_path), row.names=TRUE)

# visualize pseudobulk cell counts    
cell_count_plot <- ggplot(metadata_aggregated, aes(x = cell_count)) +
                                      geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "blue", color = NA, alpha = 0.5) +
                                      geom_density(color = "red") +
                                      theme_minimal() +
                                      labs(title = "Histogram and Density Plot of Pseudobulk Cell Counts", x = "Cell counts", y = "Density/Frequency")+
  theme(plot.title = element_text(size = 10), # Reduce title font size
        axis.title = element_text(size = 10), # Reduce axis titles font size
        axis.text = element_text(size = 10)) # Reduce axis text font size
                                      
# save plot
ggsave_new(filename=sub("\\.[[:alnum:]]+$", "", basename(cell_count_plot_path)),
           results_path=dirname(cell_count_plot_path), 
           plot=cell_count_plot, 
           width=4, 
           height=4
          )
