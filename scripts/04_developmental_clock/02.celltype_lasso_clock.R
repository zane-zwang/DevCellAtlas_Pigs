# Purpose: celltype lasso clock.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

setwd(file.path(analysis_dir, "clock/adipose"))
library(Seurat)
library(dplyr)
library(tidyverse)
library(glmnet)
library(viridis)
library(ggpubr)

set.seed(26) # set main seed
num_repeats <- 5 
lobo_results <- list() 
all_models <- list() 

# input data
clock_input_rds <- Sys.getenv("CLOCK_INPUT_RDS", "path/to/input/adipose_SCT.rds")
obj <- readRDS(clock_input_rds)
meta <- obj@meta.data %>%
  select(celltype, age) %>%
  mutate(cell_id = rownames(obj@meta.data))

# Retain cell types with more than 100 nuclei per represented stage in at least four stages.
eligible_groups <- meta %>%
  filter(!is.na(celltype), !is.na(age)) %>%
  count(celltype, age, name = "n_nuclei") %>%
  filter(n_nuclei > 100)
eligible_celltypes <- eligible_groups %>%
  count(celltype, name = "n_stages") %>%
  filter(n_stages >= 4) %>%
  pull(celltype)
eligible_groups <- eligible_groups %>% filter(celltype %in% eligible_celltypes)
eligible_cells <- meta %>%
  inner_join(eligible_groups, by = c("celltype", "age")) %>%
  pull(cell_id)
obj <- subset(obj, cells = eligible_cells)
meta <- obj@meta.data %>%
  select(celltype, age) %>%
  mutate(cell_id = rownames(obj@meta.data))

# define functions
Convert_to_Dataframe <- function(svz) {
    DefaultAssay(svz) <- "RNA"
    svz[["SCT"]] <- NULL
    meta <- svz@meta.data
    meta <- meta[, c("sample", "age", "celltype", "orig.ident")]
    raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
    raw_counts <- raw_counts[, colSums(raw_counts) > 0]
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
    pseudocells <- c()
    # If dynamic then only sample with replacement if required due to shortage of cells.
    if (replace == "dynamic") {
        if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
    }
    for (i in c(1:n)) {
        batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
        pseudocells <- rbind(pseudocells, colSums(batch))
    }
    colnames(pseudocells) <- colnames(df)
    return(as_tibble(pseudocells))
}


for (i in 1:num_repeats) {
  # set seeds
  current_seed <- 26 + i
  set.seed(current_seed)

  ## Step 1: construct pseudocells
  print(paste("Repeat", i, "- Generating pseudocells"))
  batches <- meta %>%
  group_by(celltype, age) %>%
  mutate(orig.ident = sample(c("batch1", "batch2", "batch3", "batch4", "batch5"),
                             size = n(),
                             replace = TRUE,
                             prob = c(1/5, 1/5, 1/5, 1/5, 1/5))) %>%
  ungroup()
  table(batches$orig.ident)
  obj@meta.data$orig.ident <- batches$orig.ident
  svz <- obj

  df <- Convert_to_Dataframe(svz) %>% 
    group_by(celltype, age, orig.ident, sample) %>% 
    nest()

  df <- df %>% mutate(pseudocell_all = map(data, ~ bootstrap.pseudocells(.x, size = 15, n = 100)))

  df$data <- NULL
  saveRDS(df, paste0("models_rep5/bootstrap_pseudocell_repeat_", i, ".rds"))

  ## Step 2: lognorm
  print(paste("Repeat", i, "- Log-normalizing counts"))
  lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
  }

  df <- df %>% mutate(lognorm = map(pseudocell_all, lognorm))
  df$pseudocell_all <- NULL
  df <- unnest(df, lognorm)

  ## Step 3: LOBO test
  print(paste("Repeat", i, "- Performing Leave-One-Batch-Out predictions"))
 full_df <- c()
 batches <- unique(df$orig.ident)
 Celltypes <- unique(df$celltype)

 for (CT in Celltypes) {
   print(paste("Celltype:", CT))
   df_CT <- filter(df, celltype == CT)

   for (batch in batches) {
     print(paste("  Batch:", batch))
     df_fmo <- filter(df_CT, orig.ident != batch)
     df_ss <- filter(df_CT, orig.ident == batch)

     if (length(unique(df_fmo$age)) < 2) {
       warning(paste("  Skipped:", CT, "-", batch, "- only one unique age in training set"))
       next
     }

     # train model
     model <- tryCatch({
       cv.glmnet(x = as.matrix(df_fmo[, -c(1:4)]),
                 y = as.matrix(df_fmo[, 2]),
                 type.measure = "mae", standardize = FALSE, relax = FALSE, nfolds = 5)
     }, error = function(e) {
       message(paste("  ERROR:", CT, "-", batch, "->", e$message))
       return(NULL)
     })

     # skip error
     if (is.null(model)) {
       next
     }

     # prediction
     testPredictions <- predict(model, newx = as.matrix(df_ss[, -c(1:4)]), s = "lambda.min")
     df_ss$Predictions <- testPredictions[, 1]
     test_df <- df_ss %>% select("sample", "age", "celltype", "orig.ident", "Predictions")
     colnames(test_df) <- c("sample", "age", "Celltype", "Batch", "Prediction")
     full_df <- rbind(full_df, test_df)
   }
 }

  # save LOBO results
  saveRDS(full_df, paste0("models_rep5/lobo_predictions_repeat_", i, ".rds"))
  lobo_results[[i]] <- full_df

  ## Step 4: train all data
  print(paste("Repeat", i, "- Training models by celltype"))
  df <- df %>% ungroup %>%
            select(-c(sample, orig.ident))
  by_celltype <- df %>% group_by(celltype) %>% nest()
  colnames(by_celltype)[2] <- "lognormalized"
  
  celltype_model <- function(input) {
    if (nrow(input) < 2) stop("Not enough rows for model training")
    if (ncol(input) < 2) stop("Not enough columns for model training")
    if (any(is.na(input))) stop("Input contains NA values")
    if (all(colSums(input[, -1]) == 0)) stop("All columns are zero")
    cv.glmnet(x = as.matrix(input[, -1]), y = as.matrix(input[, 1]),
              type.measure = "mae", standardize = F, relax = F, nfolds = 5)
  }

  models <- by_celltype %>% mutate(model = map(lognormalized, celltype_model))
  
  # save all data model results
  saveRDS(models, paste0("models_rep5/models_repeat_", i, ".rds"))
  all_models[[i]] <- models
}

# save data
saveRDS(all_models, "models_rep5/all_repeats_models.rds")
saveRDS(lobo_results, "models_rep5/all_repeats_lobo_results.rds")
