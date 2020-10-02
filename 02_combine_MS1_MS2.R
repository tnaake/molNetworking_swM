## combine MS2 and MS1 
## script to take MS2 spectra and trace back MS1 precursor ions
## extract the entries of these precursor ions in MS1 and write a truncated 
## peaklist 
source("01_ms2_similarity_functions.R")

## read MS1 peaklist 
setwd("/home/thomas/Projects/molNetworking_swM/results_MS1/")
## neg
pl_neg_swM <- read.table("peaklist_neg_log_tic_batch_final_swM.txt", sep = "\t")
## pos
pl_pos_swM <- read.table("peaklist_pos_log_tic_batch_final_swM.txt", sep = "\t")

## load Spectrum2 and Spectra objects
setwd("/home/thomas/Projects/molNetworking_swM/results_MS2/")
load("ms2_spectra_swM_cid_neg.RData")
load("ms2_spectra_swM_hcd_neg.RData")
load("ms2_spectra_swM_cid_pos.RData")
load("ms2_spectra_swM_hcd_pos.RData")

## sweet maize CID, negative
pl_cid_neg_map <- map_precursor(spectra_swM_cid_neg, pl_neg_swM, ppm = 20, 
    rt_tol = 0.20, rt_ms1 = "seconds")
## sweet maize HCD, negative
pl_hcd_neg_map <- map_precursor(spectra_swM_hcd_neg, pl_neg_swM, ppm = 20, 
    rt_tol = 0.20, rt_ms1 = "seconds")
## take for positive mode a higher rt_tol since retention time deviation
## occured for one sample
## sweet maize CID, positive
pl_cid_pos_map <- map_precursor(spectra_swM_cid_pos, pl_pos_swM, ppm = 20, 
    rt_tol = 0.35, rt_ms1 = "seconds")
## sweet maize HCD, positive
pl_hcd_pos_map <- map_precursor(spectra_swM_hcd_pos, pl_pos_swM, ppm = 20, 
    rt_tol = 0.35, rt_ms1 = "seconds")
    
## write the peaklist table to a file
setwd("/home/thomas/Projects/molNetworking_swM/results_MS1/")
write.table(pl_cid_neg_map, file = "peaklist_neg_swM_cid_map.txt", sep = "\t",
    quote = FALSE, dec = ".")
write.table(pl_hcd_neg_map, file = "peaklist_neg_swM_hcd_map.txt", sep = "\t",
    quote = FALSE, dec = ".")
write.table(pl_cid_pos_map, file = "peaklist_pos_swM_cid_map.txt", sep = "\t",
    quote = FALSE, dec = ".")
write.table(pl_hcd_pos_map, file = "peaklist_pos_swM_hcd_map.txt", sep = "\t",
        quote = FALSE, dec = ".")
