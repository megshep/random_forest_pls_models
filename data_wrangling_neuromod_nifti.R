library(pls)
library(RNifti)
library(dplyr)
library(ggplot2)


#DATA WRANGLING
#first thing that I need to do is load in the data into a table format
setwd("P:/pls/data")

getwd()


#load in the ELS induced amygdala change file
els_low_file <- readNifti('rbilateral_amygdala.nii')
els_low <- as.numeric(els_low_file)

# List of gene names (or identifiers)
gene_names <- c('ADRA1A', 'ADRA1B', 'ADRA1D', 'ADRA2A', 'ADRA2B', 'ADRA2C', 'ADRB1', 'ADRB2', 'ADRB3',
                '5HT1A', '5HT1B', '5HT1D', '5HT1E', '5HT1F', '5HT2A', '5HT2B', '5HT2C', '5HT3A', '5HT4', '5HT5A', '5HT6', '5HT7',
                'D1', 'D2', 'D3', 'D4', 'D5',
                'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5',
                'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6', 'CHRNA7', 'CHRNA9', 'CHRNA10',
                'CHRNB2', 'CHRNB3')

# Loop through each neuromodulator and process the corresponding NIfTI file
for (gene in gene_names) {
  # Construct the file name dynamically
  file_name <- paste0(gene, '_mirr_mRNA.nii')
  
  # Check if the file exists
  if (file.exists(file_name)) {
    # Read the NIfTI file
    nifti_data <- readNifti(file_name)
    
    # Convert the data to numeric
    assign(paste0(gene, "_map"), as.numeric(nifti_data))
    
    # Optionally, print a message confirming successful processing
    print(paste(gene, "data successfully loaded and converted to numeric."))
  } else {
    # If the file does not exist, print an error message
    print(paste(file_name, "does not exist."))
  }
}

#need to change serotonin map names to remove it starting with a 5

# Get the variable names
var_names <- ls()

# Filter only variables that start with "5HT"
matching_vars <- grep("^5HT", var_names, value = TRUE)

# Remove the '5' from the variable names and add 'HT'
new_var_names <- gsub("^5HT", "HT", matching_vars)

# Rename the variables
for(i in seq_along(matching_vars)) {
  assign(new_var_names[i], get(matching_vars[i]))
  rm(list = matching_vars[i])
}

# Check the new variable names
ls()


exp_1 <- data.frame(CHRNA2_map, CHRNA3_map, CHRNA4_map, CHRNA5_map, CHRNA6_map, CHRNA7_map, CHRNA9_map, CHRNA10_map, 
                    D1_map, D2_map, D3_map, D4_map, D5_map, 
                    ADRA1A_map, ADRA1B_map, ADRA1D_map, ADRA2A_map, ADRA2B_map, ADRA2C_map,ADRB1_map, ADRB2_map, ADRB3_map, 
                    CHRNB2_map, CHRNB3_map, 
                    HT1A_map, HT1B_map, HT1D_map, HT1E_map, HT1F_map, HT2A_map, HT2B_map, HT2C_map,
                    HT3A_map, HT4_map, HT5A_map, HT6_map, HT7_map, CHRM1_map, CHRM2_map, CHRM3_map, CHRM4_map, CHRM5_map,
                    els_low)


exp_low <- data.frame(exp_1)

#there should now be a file which has all of the ELS volume change files and the neuromod values - checking this
head(exp_low)

