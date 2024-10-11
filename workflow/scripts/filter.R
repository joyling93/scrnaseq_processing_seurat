

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
raw_object_path <- snakemake@input[["raw_object"]]

# outputs
filtered_object_path <- snakemake@output[["filtered_object"]]

# parameters
filter_expression <- snakemake@params[["filter_expression"]]
# 'flags' for modalities
ab_flag <- snakemake@params[["ab_flag"]]#'AB'
crispr_flag <- snakemake@params[["crispr_flag"]]#'gRNA'
custom_flag <- snakemake@params[["custom_flag"]]#'HTO'

### load raw data
raw_object <- readRDS(file = file.path(raw_object_path))
metadata <- raw_object[[]]

### filter data
if (filter_expression!=""){
    filter_expr <- parse(text=filter_expression)

    filtered_metadata <- subset(metadata, subset= eval(filter_expr))
    filtered_object <- raw_object[,rownames(filtered_metadata)]
}else{
    filtered_metadata <- metadata
    filtered_object <- raw_object
}


### save data
save_seurat_object(seurat_obj=filtered_object,
                   result_dir=dirname(file.path(filtered_object_path))
                  )
