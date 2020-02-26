# This script includes all custom functions for the IMPUTATION_demo.R

# This function uses a cvFolds assignment to set matrix values to missing and retains the original values
# folds = a cvFolds object that split all genotype indices into different groups
# genotypes = a numeric genotype matrix (rows = taxa, columns = markers)
# fold_ct = the number of folds to use for genotype subsetting (if NULL, use all folds)
get_gen_cv_sets <- function(folds, genotypes, fold_ct = NULL)
{
  cv_genotype_list <- list()
  
  if (!is.null(fold_ct))
  {
    iterator <- 1:fold_ct
  } else
  {
    iterator <- unique(folds$which)
  }
  
  for (i in iterator)
  {
    fold_genotypes <- unlist(genotypes) # Decompose the matrix to simplify cross-validation index reference
    fold_genotypes[folds$subsets[folds$which == i]] <- NA # Set the values to missing
    fold_genotypes <- data.frame(matrix(fold_genotypes, ncol = ncol(genotypes), nrow = nrow(genotypes))) # Recompose the genotype matrix
    names(fold_genotypes) <- names(genotypes)
    row.names(fold_genotypes) <- row.names(genotypes)
    
    fold_list <- list(fold_genotypes) 
    names(fold_list) <- c('fold_genotypes')
    
    cv_genotype_list[[i]] <- fold_list
    names(cv_genotype_list)[i] <- paste('fold_', i)
  }
  
  return(cv_genotype_list)
}

# Perform mean imputation for a given genotype table
# genotypes
impute_mean <- function(genotypes)
{
  for (i in 1:ncol(genotypes))
  {
    na_indices <- is.na(genotypes[,i])
    genotypes[na_indices,i] <- mean(genotypes[,i], na.rm = T)
  }
  return(genotypes)
}

# Calculate imptuation accuracy
calc_impute_acc <- function(genotypes, genotypes_na, genotypes_imp)
{
  na_indices <- which(is.na(genotypes_na), arr.ind = T) # Get the indices of missing values
  
  orig_values <- genotypes[na_indices]
  imp_values <- genotypes_imp[na_indices]
  
  accuracy <- cor(orig_values, imp_values)
  
  return(accuracy)
}

# This function will test the accuracy of mean imputation across a set of folds
test_mean_impute_cv <- function(genotypes, genotypes_na_folds)
{
  accuracy_vect <- vector(mode = 'numeric', length = length(genotypes_na_folds))
  
  for (i in 1:length(accuracy_vect))
  {
    print(paste("Testing fold: ", i, sep = ''))
    genotypes_na <- genotypes_na_folds[[i]][[1]] # Extract the genotypes with missing values
    genotypes_imp <- impute_mean(genotypes_na) # Perform imputation
    accuracy_vect[i] <- calc_impute_acc(genotypes, genotypes_na, genotypes_imp) # Calculate an accuracy
  }
  
  return(accuracy_vect)
}