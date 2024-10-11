#### load libraries
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

#### configs

# inputs
sample_object_paths <- snakemake@input
extra_metadata_path <- snakemake@config[["extra_metadata"]]

# outputs
merged_object <- snakemake@output[["merged_object"]]

# parameters
project_name <- snakemake@config[["project_name"]] #"test" 
result_dir <- dirname(merged_object)

# 'flags' for modalities
ab_flag <- snakemake@config[["modality_flags"]][['Antibody_Capture']]
crispr_flag <- snakemake@config[["modality_flags"]][['CRISPR_Guide_Capture']]
custom_flag <- snakemake@config[["modality_flags"]][['Custom']]


### load data
sample_objects = c()
for (sample_path in sample_object_paths){
        sample_objects <- append(sample_objects, readRDS(file = file.path(sample_path)))
}

# if more than one sample merge into one object, otherwise just rename
if (length(sample_objects)>1){
    merged_data <- merge(sample_objects[[1]],
                         y = sample_objects[2:length(sample_objects)],
                         project = project_name,
                         add.cell.ids = unlist(lapply(sample_objects, Project))
              )
}else{
    merged_data <- sample_objects[[1]]
}

# if extra metadata is provided add it to the merged object
if (extra_metadata_path != ""){
    print("extra metadata is added")
#     extra_metadata <- read.csv(extra_metadata_path, row.names = 1, header= TRUE)
    extra_metadata <- data.frame(fread(file.path(extra_metadata_path), header=TRUE), row.names=1)
    
    for (col in colnames(extra_metadata)){
        metadata_tmp <- data.frame(matrix(nrow=ncol(merged_data), ncol=1, dimnames=list(colnames(merged_data), c(col))))
        metadata_tmp[rownames(extra_metadata),col] <- extra_metadata[rownames(extra_metadata),col]
        
        # check if categorical (heuristic: less than 50 unique values) and make factor
        if (length(unique(metadata_tmp[[col]]))<50){
            metadata_tmp[, col] <- as.factor(metadata_tmp[, col])
            
            # check if any entry is an empty string i.e., "" and replace with "unknown"
            metadata_tmp[[col]][metadata_tmp[[col]] == ""] <- "unknown"
        }
        
        merged_data[[col]] <- metadata_tmp[colnames(merged_data), col]
    }
}

### save merged data
save_seurat_object(seurat_obj=merged_data,result_dir=dirname(merged_object))
