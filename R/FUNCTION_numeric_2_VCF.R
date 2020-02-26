# This function will convert a GAPIT formatted numeric genotypes with a marker map into a BEAGLE compatible VCF file (https://faculty.washington.edu/browning/beagle/beagle_5.1_08Nov19.pdf)

#' Convert numeric genotypes into BEAGLE compatible VCF format
#' 
#' @param genotypes A numeric matrix or dataframe coded as c(0,1,2). Columns = markers, rows = samples/taxa
#' @param marker_map A dataframe with the same order as the genotypes. Column names = c("rs", chr", "position")
#' 
#' @return A dataframe in VCF format with marker metadata in the first 9 columns followed by genotype calls
numeric_2_VCF <- function(genotypes, marker_map)
{
  # Set constants
  ref_default <- 'A' # Default reference base, set to A arbitrarily
  alt_default <- 'T' # Default alternate base, set to T arbitrarily
  qual_default <- 60 # Default qualtiy score, set to a very high value (in case quality matters...)
  meta_cols <- 9 # There are 9 marker data columns before the actual genotypes
  n_markers <- ncol(genotypes)
  n_samples <- nrow(genotypes)
  
  # Format the genotype table to match VCF specifications
  vcf_table <- data.frame(matrix(ncol = meta_cols+n_samples, nrow = n_markers)) # Create a data frame to store both marker metadata and genotypes
  names(vcf_table) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", row.names(genotypes))
  
  # Convert values into VCF string format
  vcf_genotypes <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes))
  vcf_genotypes[genotypes == 0] <- "0/0"
  vcf_genotypes[genotypes == 1] <- "0/1"
  vcf_genotypes[genotypes == 2] <- "1/1"
  vcf_genotypes[is.na(genotypes)] <- "./."
  
  # Transpose the genotype matrix to match the VCF
  vcf_genotypes <- t(vcf_genotypes) 
  
  # Insert the VCF calls into the VCF table
  vcf_table[,(meta_cols + 1):ncol(vcf_table)] <- vcf_genotypes
  
  # Add in metadata from the map table
  vcf_table$CHROM <- marker_map$chr
  vcf_table$ID <- marker_map$rs
  vcf_table$POS <- marker_map$pos
  vcf_table$REF <- ref_default
  vcf_table$ALT <- alt_default
  vcf_table$QUAL <- qual_default
  vcf_table$FILTER <- "PASS"
  vcf_table$INFO <- "None"
  vcf_table$FORMAT <- "GT"
  
  return(vcf_table)
}

#' This function will write a VCF to a text file in the proper format
#' 
#' @param vcf_genotypes A dataframe in VCF format
#' @param outfile The name of the output .vcf file
#' 
#' @return 0 if file was created, 1 if file was not created
write_vcf <- function(vcf_genotypes, outfile)
{
  # Set constants
  VCF_version = 4.2 # (This might be incorporated as a feature in the future for different methods for different versions)
  
  sink(paste(outfile, ".vcf", sep = "")) # Redirect all R output to the specified file
  cat(paste("##fileformat=VCFv", VCF_version, "\n", sep = '')) # Write the first line representing theversion
  cat('#')
  sink() # Close redirection
  
  write.table(vcf_genotypes, 
              file = paste(outfile, ".vcf", sep = ""), 
              append = T, 
              quote = F, 
              sep = "\t", 
              na = ".", 
              row.names = F, 
              col.names = T)
  
  if (file.exists(paste(outfile, ".vcf", sep = ""))) {
    return(0)
  } else
  {
    return(1)
  }
}

#' This function will write a VCF to a text file in the proper format
#' 
#' @param infile The name of the .vcf file to import
#' @param skip The number of metadata lines in the .vcf
#' 
#' @return A VCF dataframe
read_vcf <- function(infile, skip = 8)
{
  vcf_table <- read.table(infile, 
                          header = T, 
                          skip = skip, 
                          stringsAsFactors = F, 
                          comment.char = "", 
                          na = ".")
  names(vcf_table)[1] <- "CHROM"
  
  return(vcf_table)
}

