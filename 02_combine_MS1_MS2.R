## combine MS2 and MS1 
## script to take MS2 spectra and trace back MS1 precursor ions
## extract the entries of these precursor ions in MS1 and write a truncated 
## peaklist 

## read MS1 peaklist 
setwd("/home/thomas/Projects/molNetworking_swM/results_MS1/")
## neg
pl_neg_swM <- read.table("", sep = "\t")
## pos
pl_pos_swM <- read.table("", sep = "\t")

peaklist <- read.table("peaklist_log_batch_tic_final.txt", sep = "\t")


## load Spectrum2 and Spectra objects
setwd("/home/thomas/Projects/molNetworking_swM/results_MS2/")
load("ms2_spectra_swM_cid_neg.RData")
load("ms2_spectra_swM_hcd_neg.RData")
load("ms2_spectra_swM_cid_pos.RData")
load("ms2_spectra_swM_hcd_pos.RData")

## iterate through each Spectrum2 
spectra <- spectra_swM_cid_neg
#' @name map_precursor
#' @title Map precursors of MS2 to MS1
#' @description The function `map_precursor` takes a `Spectra` objects
#' containing MS2 information and a peaklist `matrix` with MS1 information and 
#' returns those features in the MS1 `matrix` whose features with m/z and 
#' retention time pairs match to the precursor ions of the `Spectra` object. 
#' @param spectra `Spectra` object containing MS2 information
#' @param peaklist `matrix`, containing the columns `"mz"` and `"rt"`
#' @param ppm `numeric`, m/z tolerance in ppm
#' @param rt_tol `numeric`, retention time tolerance in minutes
#' @details If several features in `peaklist` match to a MS2 spectrum, the 
#' function `map_precursor` will return all matching features. 
#' @usage map_precursor(spectra, peaklist, ppm = 20, rt_tol = 0.15)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @return `matrix`
#' @examples 
#' map_precursor(spectra, peaklist, ppm = 20, rt_tol = 0.15)
map_precursor <- function(spectra, peaklist, ppm = 20, rt_tol = 0.15, rt_ms1 = c("seconds", "minutes")) {
    

    if (rt_ms1 == "seconds") {
        peaklist[, "rt"] <- peaklist[, "rt"] / 60
    }
    
    res <- vector("list", nrow(peaklist))
    ## create vector that stores retention time of MS1 
    ms1_rt <- peaklist[, "rt"]
    
    pl_map <- lapply(seq_along(spectra), function(x) {
        
        ## retrieve precursorMz and retention time of Spectra 
        mz <- spectra[[x]]@precursorMz
        rt <- spectra[[x]]@rt
        
        ## restrict the search space based on retention time
        pl <- peaklist[ms1_rt >= (rt - rt_tol) & ms1_rt <= (rt + rt_tol),]
        
        ## calculate lower and upper m/z based on ppm
        upper <- mz / abs(ppm / 10 ^ 6  - 1 ) 
        lower <- mz / abs(ppm / 10 ^ 6  + 1 ) 
        
        ## restrict the search space based on m/z
        ms1_mz <- pl[, "mz"]
        pl <- pl[ms1_mz >= lower & ms1_mz <= upper, ]
        
        ## add information of the mapped Spectrum2
        if (nrow(pl) > 0) pl <- cbind(spectra=names(spectra)[x], pl)
        
        ## return
        return(pl)
    })
    
    ## rbind the list 
    do.call("rbind", pl_map) 
}  

## sweet maize CID, negative
pl_cid_neg_map <- map_precursor(spectra_swM_cid_neg, pl_neg_swM, ppm = 20, 
    rt_tol = 0.15, rt_ms1 = "seconds")
## sweet maize HCD, negative
pl_hcd_neg_map <- map_precursor(spectra_swM_hcd_neg, pl_neg_swM, ppm = 20, 
    rt_tol = 0.15, rt_ms1 = "seconds")
## sweet maize CID, positive
pl_cid_pos_map <- map_precursor(spectra_swM_cid_pos, pl_pos_swM, ppm = 20, 
    rt_tol = 0.15, rt_ms1 = "seconds")
## sweet maize HCD, positive
pl_hcd_pos_map <- map_precursor(spectra_swM_hcd_pos, pl_pos_swM, ppm = 20, 
    rt_tol = 0.15, rt_ms1 = "seconds")
