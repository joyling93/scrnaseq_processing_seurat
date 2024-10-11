
#### load libraries & utility function 
library("Seurat")

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

# inputs
merged_object_path <- snakemake@input[["merged_object"]]

# outputs
split_object_path <- snakemake@output[["split_object"]]

# parameters
split_by <- snakemake@wildcards[["split"]]
# 'flags' for modalities
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]

### load merged data
merged_object <- readRDS(file = file.path(merged_object_path))
metadata <- merged_object[[]]


### split data
split_list <- unlist(regmatches(split_by, regexpr("__", split_by), invert = TRUE))
split <- split_list[1]
cat <- split_list[2]

split_expr <- parse(text=paste0(split,'==',deparse(cat)))

tmp_metadata <- subset(metadata, subset= eval(split_expr))
tmp_object <- merged_object[,rownames(tmp_metadata)]

### save data
save_seurat_object(seurat_obj=tmp_object, result_dir=dirname(split_object_path))
