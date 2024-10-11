
#### load libraries & utility function 
library("inspectdf")
library("tibble")
library("ggplot2")

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
metadata_path <- snakemake@input[["metadata"]]

# outputs
result_dir <- snakemake@output[["metadata_plots"]]
stats_dir <- snakemake@output[["metadata_stats"]]

# load data
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)
metadata <- as_tibble(metadata)

### plot metadata types
metadata_type_stats <- inspect_types(metadata)
tmp_plot <- show_plot(metadata_type_stats)

# plots specs
width <- 5
height <- 5

# save plot
ggsave_new(filename="types", 
           results_path=result_dir, 
           plot=tmp_plot, 
           width=width, 
           height=height
          )

### categorical data
# inspect categorical data
metadata_cat_stats <- inspect_cat(metadata)

# plot stats
tmp_plot <- metadata_cat_stats %>% show_plot(label_thresh=0.01)

# plots specs
height <- nrow(metadata_cat_stats)*0.5
width <- 10

# save plot
ggsave_new(filename="categorical", 
           results_path=result_dir, 
           plot=tmp_plot, 
           width=width, 
           height=height
          )


### numerical data
# inspect numerical data
metadata_num_stats <- inspect_num(metadata, breaks=100)

# plot stats
tmp_plot <- metadata_num_stats %>% show_plot()

# plots specs
height <- nrow(metadata_num_stats)/3*1.5

# save plot
ggsave_new(filename="numerical",
           results_path=result_dir,
           plot=tmp_plot, 
           width=width, 
           height=height
          )


### save all statistics as CSV files
dir.create(stats_dir, recursive = TRUE)

for (cat in names(metadata_cat_stats$levels)){
#     write.csv(metadata_cat_stats$levels[[cat]], file=file.path(stats_dir, paste0(step,"_metadata_",cat,".csv")), row.names=FALSE)
    fwrite(as.data.frame(metadata_cat_stats$levels[[cat]]), file=file.path(stats_dir, paste0("metadata_",cat,".csv")), row.names=FALSE)
}

# write.csv(metadata_num_stats[,-ncol(metadata_num_stats)], file=file.path(stats_dir, paste0(step,"_metadata_","numerical",".csv")), row.names=FALSE)
fwrite(as.data.frame(metadata_num_stats[,-ncol(metadata_num_stats)]), file=file.path(stats_dir, "metadata_numerical.csv"), row.names=FALSE)
