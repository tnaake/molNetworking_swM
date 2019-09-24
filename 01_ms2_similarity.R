################################################################################
################################## functions ###################################
################################################################################

## functions to create MS/MS features from corresponding files that were 
## collected under different collision energies

#' @name assemblySpectra
#' @title Deconvolute final spectra from different collision energies
#' @description assemblySpectra takes files from different collision energies
#' as input and assembles a final deconvoluted spectra from it, taking the 
#' highest intensity values from corresponding peaks across the different 
#' collision energies. 
#' @usage assemlySpectra(spectra, aln, cols_aln, sample_name, tol)
#' @param spectra list of matrices, each object is from MS-DIAL and contains
#' information on the MS/MS spectrum
#' @param aln matrix, alignment file from MS-DIAL
#' @param cols_aln character, colnames in aln that correspond to the objects in
#' spectra
#' @param sample_name character, add the name to the names of the returned list
#' (Alignment.ID, retention time and m/z of the precursor, taken from aln, will
#' be automatically added to the name)
#' @param tol numeric, tolerance parameter
#' @details assemblySpectra will first bin within each spectra (su)
#' @return list of assembled spectra. Each spectra is a matrix that in the 
#' first column contains m/z values and in the second column contains the 
#' corresponding intensities
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples 
#' sp <- list(MK_hcd30, MK_hcd40, MK_hcd50)
#' cols_aln <- c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50")
#' assemblySpectra(spectra=sp, aln=aln, cols_aln=cols_aln, sample_name="MK_HCD")
assemblySpectra <- function(spectra=list(MK_hcd30, MK_hcd40, MK_hcd50), aln=aln,
    cols_aln=c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD", tol=0.01) {
    
    
    ## create list of vector nrow(aln) that will store the final deconvoluted 
    ## spectra, names refer to the aligned spectra in aln
    res <- vector("list", nrow(aln))
    
    
    ## iterate through aln and assemble final deconvoluted spectra
    res <- lapply(seq_len(nrow(aln)), function(i) {
        
        ## retrieve PeakID from alignment object aln per row, get IDs from all
        ## collision energies and assign as a list to inds_ev
        inds_ev <- lapply(cols_aln, function(x) aln[i, x])
        
        ## retrieve the corresponding row in spectra (for each collision energy)
        ## the entry will store information on the spectrum, use the "PeakID"
        ## as the identifier
        i_spectra <- lapply(seq_along(spectra), function(x) {
            if (inds_ev[[x]] != "-2") {
                sp <- spectra[[x]][spectra[[x]][, "PeakID"] == inds_ev[[x]],]
                ## if there is no matching MS/MS spectrum recorded set to NULL
                if (nrow(sp) == 0) sp <- NULL
            } else sp <- NULL
            return(sp)
        })
        
        ## if not all entries are NULL
        if (!is.null(unlist(i_spectra))) {
        
            ## retrieve the column "MSMS.spectrum" from each entry in i_spectra,
            ## this column contains information on the peaks and the 
            ## corresponding intensities, strsplit the entries and write them 
            ## to a matrix
            i_msms <- lapply(i_spectra, function(x) {
                if (!is.null(x)) {
                    msms <- strsplit(x[, "MSMS.spectrum"], split=" ")[[1]]  
                    msms <- lapply(msms, function(y) strsplit(y, split=":"))
                    msms <- matrix(unlist(msms), ncol=2, byrow=TRUE)
                    mode(msms) <- "numeric"
                    return(msms)
                }
            })
        
            ## bin within each spectra of i_spectra_msms
            ## sum up corresponding intensities from peaks that are  binned
            ## together
            i_msms_bin <- lapply(i_msms, function(x) {
                if (!is.null(x)) {
                    binSpectra(x, tol=tol, fun="sum")   
                }
            })
        
            ## deconvolute final spectra, bin across spectra
            deconv <- deconvolute(i_msms_bin, tol=tol)
        
        } else deconv <- NULL
        
        return(deconv)
    })
    
    ## construct names for res
    names(res) <- paste(sample_name, aln[, "Alignment.ID"], 
        aln[, "Average.Rt.min."], aln[, "Average.Mz"], sep="_")
    
    ## truncate res that it only contains the entries that are not NULL
    inds_keep <- unlist(lapply(res, function(x) !is.null(x)))
    res <- res[inds_keep]
    
    return(res)
}


#' @name binSpectra
#' @title Bin a spectra
#' @description \code{binSpectra} combines peaks that are close to each 
#' other
#' @usage binSpectra(spectra, tol)
#' @param spectra matrix, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#' @param tol numeric, tolerance parameter
#' @param fun character, either "max" or "sum"
#' @details Bin peaks and sum the intensities of corresponding peaks or 
#' take the maximum (depending on \code{fun}). If \code{fun} is equal to 
#' "max" the m/z value of the peak with the highest intensity is 
#' retained togheter with the highest intensity, if \code{fun} is equal to 
#' "sum" the median of m/z values in the bin are retained and intensities 
#' are added. 
#' @return matrix containing in the first column m/z information and in the 
#' second column intensity information
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples
#' binSpectra(spectra, tol=0.01, fun="sum")
binSpectra <- function(spectra, tol=0.01, fun=c("max", "sum")) {
    
    ##if (!is.matrix(spectra)) stop("spectra is not a matrix")
    ## match arguments for fun
    fun <- match.arg(fun)
    
    ## if spectra only contains > 2 do the binning, otherwise return spectra
    if (nrow(spectra) > 2) {
        frag_s <- spectra[,1]
        steps <- (max(frag_s) - min(frag_s)) / tol
        if (steps > 2) {
            bins <- tapply(frag_s, cut(frag_s, steps), mean)
            bins <- bins[!is.na(bins)]
            ## find nearest ones and sum all intensities up 
            inds <- lapply(frag_s, FUN = function(x) which.min(abs(x - bins)))
            inds <- unlist(inds)
    
            ## iterate through duplicated peaks and combine them
            for (j in names(which(table(inds) != 1))) {
                inds_dup <- which(inds == j)
                spectra_dup <- spectra[inds_dup,]
        
                ## either use max or sum the intensities
                if (fun == "max") {
                    spectra[inds_dup, 1] <- spectra_dup[which.max(spectra_dup[, 2]), 1]
                    spectra[inds_dup, 2] <- max(spectra_dup[, 2])    
                    ## set all except the ones with the highest intensity to NA
                    spectra[inds_dup[-which.max(spectra_dup[, 2])], ] <- NA 
                } 
                if (fun == "sum") {
                    spectra[inds_dup, 1] <- median(spectra_dup[, 1])
                    spectra[inds_dup, 2] <- sum(spectra_dup[, 2])   
                    ## set all except the first one to NA
                    spectra[inds_dup[-1],] <- NA 
                }
            }
            ## remove NA values
            spectra <- spectra[!is.na(spectra[,1]), ]
        } else {
            spectra <- matrix(c(frag_s[which.max(spectra[, 2])], max(spectra[, 2])), ncol=2)
        }
    }
    return(spectra)
}

#' @name deconvolute
#' @title Deconvolute a final spectrum
#' @description \code{deconvolute} takes as input a list of MS/MS spectra. It will 
#' bin the m/z values and take the maximum of the corresponding intensity value.
#' @usage deconvolute(spectra, tol=0.01)
#' @param spectra list of matrices, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#' @param tol numeric, tolerance parameter for binning
#' @details Takes a list of matrices that contain spectral information and 
#' bins the m/z values. If several peaks are located in the same bin, 
#' the peak with the highest intensity is retained and others are removed
#' from the spectrum. 
#' @return matrix containing in the first column m/z information and in the 
#' second column intensity information
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples
#' sp <- list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50)
#' deconvolute(spectra=sp, tol=0.01)
deconvolute <- function(spectra=list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50), tol=0.01) {
    ## rbind spectra
    spectra <- do.call("rbind", spectra)
    
    ## order spectra
    spectra <- spectra[order(spectra[,1]),]
    
    ## convert to matrix
    if (!is.matrix(spectra)) spectra <- matrix(spectra, ncol = 2, byrow=TRUE)
    
    ## bin: take the max from the same corresonding peak
    spectra <- binSpectra(spectra, tol=tol, fun="max")
    
    ## return
    return(spectra)
}

## use other function from MetCirc (>v1.15.0) and MSnbase (>v2.11.4) and 
## igraph (>v1.2.4.1)for visualization 
library(MSnbase)
library(igraph)

################################################################################
################################## workflow ####################################
################################################################################
## Maize MS2: create a similarity network
setwd("/home/thomas/Projects/molNetworking_swM/data_MS2/neg")

## 1) 
################################### negative ###################################
## load file that contains alignment information
aln_neg <- read.table("PeakID_0_2019822154.txt", sep="\t", fill=TRUE, skip=4, 
                      header=TRUE, stringsAsFactors=FALSE, quote="\"")
## truncate file
## MK: maize kernel, ML: maize leaf, QC: sweet maize
cols_keep <- c("Alignment.ID", "Average.Rt.min.", "Average.Mz", 
    "Metabolite.name", "Adduct.type", 
    "mixMK.ddms2.CID30", "mixMK.ddms2.CID40", 
    "mixMK.ddms2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50",
    "mixML.ddms2.CID30", "mixML.ddms2.CID40",
    "mixML.ddms2.HCD30", "mixML.ddms2.HCD40", "mixML.ddms2.HCD50",
    "QC.ddms.CID30", "QC.ddms.CID40", 
    "QC.ddms.HCD30", "QC.ddms.HCD40", "QC.ddms.HCD50")

         
aln_neg <- aln_neg[, cols_keep]
## keep only these rows that have alignment information
inds_remove <- apply(aln_neg[, cols_keep[-c(1:5)]], 1, function(x) all(x == "-2"))
aln_neg <- aln_neg[!inds_remove, ]
## remove lines that contain any NA
aln_neg <- aln_neg[!apply(aln_neg, 1, function(x) any(is.na(x))), ]

## 2) 
## load files that contain spectra, do this for each sample individually
## for CID only 30 and 40 eV available, for HCD 30, 40 and 50 eV available

## for kernel CID
MK_cid30_neg <- read.table("mixMK-ddms2-CID30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_cid40_neg <- read.table("mixMK-ddms2-CID40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for kernel HCD
MK_hcd30_neg <- read.table("mixMK-ddms2-HCD30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd40_neg <- read.table("mixMK-ddms2-HCD40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd50_neg <- read.table("mixMK-ddms2-HCD50.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for leaf CID
ML_cid30_neg <- read.table("mixML-ddms2-CID30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_cid40_neg <- read.table("mixML-ddms2-CID40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for leaf HCD
ML_hcd30_neg <- read.table("mixML-ddms2-HCD30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd40_neg <- read.table("mixML-ddms2-HCD40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd50_neg <- read.table("mixML-ddms2-HCD50.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for sweet maize CID
swM_cid30_neg <- read.table("QC-ddms-CID30.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_cid40_neg <- read.table("QC-ddms-CID40.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for sweet maize HCD
swM_hcd30_neg <- read.table("QC-ddms-HCD30.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd40_neg <- read.table("QC-ddms-HCD40.txt", sep="\t", fill=TRUE,
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd50_neg <- read.table("QC-ddms-HCD50.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")

################################### positive ###################################
## load file that contains alignment information
setwd("../pos")
aln_pos <- read.table("PeakID_0_20198221558.txt", sep="\t", fill=TRUE, 
                      skip=4, header=TRUE, stringsAsFactors=FALSE, quote="'")
## truncate file
## MK: maize kernel, ML: maize, QCsm: sweet maize
cols_keep <- c("Alignment.ID", "Average.Rt.min.", "Average.Mz", 
    "Metabolite.name", "Adduct.type",
    "mixMK.ddms2.CID30" ,"mixMK.ddms2.CID40",
    "mixMK.ddms2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50",
    "mixML.pos.ddms2.CID30", "mixML.pos.ddms2.CID40", 
    "mixML.pos.ddms2.HCD30", "mixML.pos.ddms2.HCD40", "mixML.pos.ddms2.HCD50",
    "QCsm.pos.ddms.CID30", "QCsm.pos.ddms.CID40",
    "QCsm.pos.ddms.HCD30", "QCsm.pos.ddms.HCD40", "QCsm.pos.ddms.HCD50")

aln_pos <- aln_pos[, cols_keep]
## keep only these rows that have alignment information
inds_remove <- apply(aln_pos[,cols_keep[-c(1:5)]], 1, function(x) all(x == "-2"))
aln_pos <- aln_pos[!inds_remove, ]
## remove lines that contain any NA
aln_pos <- aln_pos[!apply(aln_pos, 1, function(x) any(is.na(x))), ]

## 2) 
## load files that contain spectra, do this for each sample individually
## for CID only 30 and 40 eV available, for HCD 30, 40 and 50 eV available
## for kernel CID
MK_cid30_pos <- read.table("mixMK-ddms2-CID30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_cid40_pos <- read.table("mixMK-ddms2-CID40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for kernel HCD
MK_hcd30_pos <- read.table("mixMK-ddms2-HCD30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd40_pos <- read.table("mixMK-ddms2-HCD40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
MK_hcd50_pos <- read.table("mixMK-ddms2-HCD50.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for leaf CID
ML_cid30_pos <- read.table("mixML-pos-ddms2-CID30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_cid40_pos <- read.table("mixML-pos-ddms2-CID40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for leaf HCD
ML_hcd30_pos <- read.table("mixML-pos-ddms2-HCD30.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd40_pos <- read.table("mixML-pos-ddms2-HCD40.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")
ML_hcd50_pos <- read.table("mixML-pos-ddms2-HCD50.txt", sep="\t", fill=TRUE, 
                           header=TRUE, stringsAsFactors=FALSE, quote="\"")

## for sweet maize CID
swM_cid30_pos <- read.table("QCsm-pos-ddms-CID30.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_cid40_pos <- read.table("QCsm-pos-ddms-CID40.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
## for sweet maize HCD
swM_hcd30_pos <- read.table("QCsm-pos-ddms-HCD30.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd40_pos <- read.table("QCsm-pos-ddms-HCD40.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")
swM_hcd50_pos <- read.table("QCsm-pos-ddms-HCD50.txt", sep="\t", fill=TRUE, 
                            header=TRUE, stringsAsFactors=FALSE, quote="\"")

## 3) 
## truncate files and keep only columns that are of interest
cols_keep <- c("PeakID", "Title", "RT.min.", "Precursor.m.z", "Height", 
    "Area", "MetaboliteName", "AdductIon", "Isotope", "S.N", "MSMS.spectrum")

## neg
MK_cid30_neg <- MK_cid30_neg[MK_cid30_neg[, "MSMS.spectrum"] != "", cols_keep]
MK_cid40_neg <- MK_cid40_neg[MK_cid40_neg[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd30_neg <- MK_hcd30_neg[MK_hcd30_neg[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd40_neg <- MK_hcd40_neg[MK_hcd40_neg[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd50_neg <- MK_hcd50_neg[MK_hcd50_neg[, "MSMS.spectrum"] != "", cols_keep]
ML_cid30_neg <- ML_cid30_neg[ML_cid30_neg[, "MSMS.spectrum"] != "", cols_keep]
ML_cid40_neg <- ML_cid40_neg[ML_cid40_neg[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd30_neg <- ML_hcd30_neg[ML_hcd30_neg[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd40_neg <- ML_hcd40_neg[ML_hcd40_neg[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd50_neg <- ML_hcd50_neg[ML_hcd50_neg[, "MSMS.spectrum"] != "", cols_keep]
swM_cid30_neg <- swM_cid30_neg[swM_cid30_neg[, "MSMS.spectrum"] != "", cols_keep]
swM_cid40_neg <- swM_cid40_neg[swM_cid40_neg[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd30_neg <- swM_hcd30_neg[swM_hcd30_neg[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd40_neg <- swM_hcd40_neg[swM_hcd40_neg[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd50_neg <- swM_hcd50_neg[swM_hcd50_neg[, "MSMS.spectrum"] != "", cols_keep]

## pos
MK_cid30_pos <- MK_cid30_pos[MK_cid30_pos[, "MSMS.spectrum"] != "", cols_keep]
MK_cid40_pos <- MK_cid40_pos[MK_cid40_pos[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd30_pos <- MK_hcd30_pos[MK_hcd30_pos[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd40_pos <- MK_hcd40_pos[MK_hcd40_pos[, "MSMS.spectrum"] != "", cols_keep]
MK_hcd50_pos <- MK_hcd50_pos[MK_hcd50_pos[, "MSMS.spectrum"] != "", cols_keep]
ML_cid30_pos <- ML_cid30_pos[ML_cid30_pos[, "MSMS.spectrum"] != "", cols_keep]
ML_cid40_pos <- ML_cid40_pos[ML_cid40_pos[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd30_pos <- ML_hcd30_pos[ML_hcd30_pos[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd40_pos <- ML_hcd40_pos[ML_hcd40_pos[, "MSMS.spectrum"] != "", cols_keep]
ML_hcd50_pos <- ML_hcd50_pos[ML_hcd50_pos[, "MSMS.spectrum"] != "", cols_keep]
swM_cid30_pos <- swM_cid30_pos[swM_cid30_pos[, "MSMS.spectrum"] != "", cols_keep]
swM_cid40_pos <- swM_cid40_pos[swM_cid40_pos[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd30_pos <- swM_hcd30_pos[swM_hcd30_pos[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd40_pos <- swM_hcd40_pos[swM_hcd40_pos[, "MSMS.spectrum"] != "", cols_keep]
swM_hcd50_pos <- swM_hcd50_pos[swM_hcd50_pos[, "MSMS.spectrum"] != "", cols_keep]

## 4) 
## truncate files that they only contain the spectra that are aligned
## (spectra that are present in aln)

## neg
## kernel 
MK_cid30_neg <- MK_cid30_neg[MK_cid30_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.CID30"],]
MK_cid40_neg <- MK_cid40_neg[MK_cid40_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.CID40"],]
MK_hcd30_neg <- MK_hcd30_neg[MK_hcd30_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddm2.HCD30"],]
MK_hcd40_neg <- MK_hcd40_neg[MK_hcd40_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.HCD40"],]
MK_hcd50_neg <- MK_hcd50_neg[MK_hcd50_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.HCD50"],]

## leaf 
ML_cid30_neg <- ML_cid30_neg[ML_cid30_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.CID30"],]
ML_cid40_neg <- ML_cid40_neg[ML_cid40_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.CID40"],]
ML_hcd30_neg <- ML_hcd30_neg[ML_hcd30_neg[, "PeakID"] %in% aln_neg[, "mixML.ddm2.HCD30"],]
ML_hcd40_neg <- ML_hcd40_neg[ML_hcd40_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.HCD40"],]
ML_hcd50_neg <- ML_hcd50_neg[ML_hcd50_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.HCD50"],]

## sweet maize
swM_cid30_neg <- swM_cid30_neg[swM_cid30_neg[, "PeakID"] %in% aln_neg[, "QC.ddms.CID30"],]
swM_cid40_neg <- swM_cid40_neg[swM_cid40_neg[, "PeakID"] %in% aln_neg[, "QC.ddms.CID40"],]
swM_hcd30_neg <- swM_hcd30_neg[swM_hcd30_neg[, "PeakID"] %in% aln_neg[, "QC.ddms.HCD30"],]
swM_hcd40_neg <- swM_hcd40_neg[swM_hcd40_neg[, "PeakID"] %in% aln_neg[, "QC.ddms.HCD40"],]
swM_hcd50_neg <- swM_hcd50_neg[swM_hcd50_neg[, "PeakID"] %in% aln_neg[, "QC.ddms.HCD50"],]

## pos
## kernel 
MK_cid30_pos <- MK_cid30_pos[MK_cid30_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.CID30"],]
MK_cid40_pos <- MK_cid40_pos[MK_cid40_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.CID40"],]
MK_hcd30_pos <- MK_hcd30_pos[MK_hcd30_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddm2.HCD30"],]
MK_hcd40_pos <- MK_hcd40_pos[MK_hcd40_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.HCD40"],]
MK_hcd50_pos <- MK_hcd50_pos[MK_hcd50_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.HCD50"],]

## leaf 
ML_cid30_pos <- ML_cid30_pos[ML_cid30_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.CID30"],]
ML_cid40_pos <- ML_cid40_pos[ML_cid40_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.CID40"],]
ML_hcd30_pos <- ML_hcd30_pos[ML_hcd30_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddm2.HCD30"],]
ML_hcd40_pos <- ML_hcd40_pos[ML_hcd40_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.HCD40"],]
ML_hcd50_pos <- ML_hcd50_pos[ML_hcd50_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.HCD50"],]

## sweet maize
swM_cid30_pos <- swM_cid30_pos[swM_cid30_pos[, "PeakID"] %in% aln_pos[, "QCsm.pos.ddms.CID30"],]
swM_cid40_pos <- swM_cid40_pos[swM_cid40_pos[, "PeakID"] %in% aln_pos[, "QCsm.pos.ddms.CID40"],]
swM_hcd30_pos <- swM_hcd30_pos[swM_hcd30_pos[, "PeakID"] %in% aln_pos[, "QCsm.pos.ddms.HCD30"],]
swM_hcd40_pos <- swM_hcd40_pos[swM_hcd40_pos[, "PeakID"] %in% aln_pos[, "QCsm.pos.ddms.HCD40"],]
swM_hcd50_pos <- swM_hcd50_pos[swM_hcd50_pos[, "PeakID"] %in% aln_pos[, "QCsm.pos.ddms.HCD50"],]

## 5) 
## assemble spectra
## neg
assembly_MK_cid_neg <- assemblySpectra(spectra=list(MK_cid30_neg, MK_cid40_neg), 
    aln=aln_neg, cols_aln=c("mixMK.ddm2.CID30", "mixMK.ddms2.CID40"),
    sample_name="MK_CID")
assembly_MK_hcd_neg <- assemblySpectra(spectra=list(MK_hcd30_neg, MK_hcd40_neg, MK_hcd50_neg), 
    aln=aln_neg, cols_aln=c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD")
assembly_ML_cid_neg <- assemblySpectra(spectra=list(ML_cid30_neg, ML_cid40_neg), 
    aln=aln_neg, cols_aln=c("mixML.ddm2.CID30", "mixML.ddms2.CID40"),
    sample_name="ML_CID")
assembly_ML_hcd_neg <- assemblySpectra(spectra=list(ML_hcd30_neg, ML_hcd40_neg, ML_hcd50_neg), 
    aln=aln_neg, cols_aln=c("mixML.ddm2.HCD30", "mixML.ddms2.HCD40", "mixML.ddms2.HCD50"),
    sample_name="ML_HCD")
assembly_swM_cid_neg <- assemblySpectra(spectra=list(swM_cid30_neg, swM_cid40_neg), ## 1359 1319
    aln=aln_neg, cols_aln=c("QC.ddms.CID30", "QC.ddms.CID40"),
    sample_name="swM_CID")
assembly_swM_hcd_neg <- assemblySpectra(spectra=list(swM_hcd30_neg, swM_hcd40_neg, swM_hcd50_neg), 
    aln=aln_neg, cols_aln=c("QC.ddms.HCD30", "QC.ddms.HCD40", "QC.ddms.HCD50"),
    sample_name="swM_HCD")

## save
setwd("/home/thomas/Projects/molNetworking_swM/results_MS2/")
save(assembly_MK_cid_neg, file="ms2_assembly_MK_cid_neg.RData")
save(assembly_MK_hcd_neg, file="ms2_assembly_MK_hcd_neg.RData")
save(assembly_ML_cid_neg, file="ms2_assembly_ML_cid_neg.RData")
save(assembly_ML_hcd_neg, file="ms2_assembly_ML_hcd_neg.RData")
save(assembly_swM_cid_neg, file="ms2_assembly_swM_cid_neg.RData")
save(assembly_swM_hcd_neg, file="ms2_assembly_swM_hcd_neg.RData")

## pos
assembly_MK_cid_pos <- assemblySpectra(spectra=list(MK_cid30_pos, MK_cid40_pos), 
    aln=aln_pos, cols_aln=c("mixMK.ddm2.CID30", "mixMK.ddms2.CID40"),
    sample_name="MK_CID")
assembly_MK_hcd_pos <- assemblySpectra(spectra=list(MK_hcd30_pos, MK_hcd40_pos, MK_hcd50_pos), 
    aln=aln_pos, cols_aln=c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD")
assembly_ML_cid_pos <- assemblySpectra(spectra=list(ML_cid30_pos, ML_cid40_pos), 
    aln=aln_pos, cols_aln=c("mixML.ddm2.CID30", "mixML.ddms2.CID40"),
    sample_name="ML_CID")
assembly_ML_hcd_pos <- assemblySpectra(spectra=list(ML_hcd30_pos, ML_hcd40_pos, ML_hcd50_pos), 
    aln=aln_pos, cols_aln=c("mixML.ddm2.HCD30", "mixML.ddms2.HCD40", "mixML.ddms2.HCD50"),
    sample_name="ML_HCD")
assembly_swM_cid_pos <- assemblySpectra(spectra=list(swM_cid30_pos, swM_cid40_pos), 
    aln=aln_pos, cols_aln=c("QCsm.pos.ddms.CID30", "QCsm.pos.ddms.CID40"),
    sample_name="swM_CID")
assembly_swM_hcd_pos <- assemblySpectra(spectra=list(swM_hcd30_pos, swM_hcd40_pos, swM_hcd50_pos), 
    aln=aln_pos, cols_aln=c("QCsm.pos.ddms.HCD30", "QCsm.pos.ddms.HCD40", "QCsm.pos.ddms.HCD50"),
    sample_name="swM_HCD")

## save
save(assembly_MK_cid_pos, file="ms2_assembly_MK_cid_pos.RData")
save(assembly_MK_hcd_pos, file="ms2_assembly_MK_hcd_pos.RData")
save(assembly_ML_cid_pos, file="ms2_assembly_ML_cid_pos.RData")
save(assembly_ML_hcd_pos, file="ms2_assembly_ML_hcd_pos.RData")
save(assembly_swM_cid_pos, file="ms2_assembly_swM_cid_pos.RData")
save(assembly_swM_hcd_pos, file="ms2_assembly_swM_hcd_pos.RData")

## 6) 
## convert final deconvoluted spectra to Spectrum2 objects from MSnbase
## load the MSnbase package
library(MSnbase)

## retrieve retention time, precursorMZ, m/z of peaks and corresponding 
## intensities from the assembled spectrum

#' @name construct_Spectrum2
#' @description The function `construct_Spectrum2` retrieves the relevant 
#' information from the `assembly` (retention time of precursor, 
#' m/z of precursor and m/z and intensities of the fragments) 
#' and creates a list of `Spectrum2` objects.
#' @param assembly list of matrices
#' @details The function `construct_Spectrum2` takes as an argument a list of
#' matrices. Each matrix contains in its first column m/z values and in the 
#' second column intensity values. 
#' @usage construct_Spectrum2(assembly)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @return `list` of `Spectrum2` objects
#' @examples 
#' construct_Spectrum(assembly)
construct_Spectrum2 <- function(assembly) {
    
    ## get retention time
    rt <- unlist(lapply(strsplit(names(assembly), split="_"), "[", 4))  
    rt <- as.numeric(rt)
    ## get precursor m/z
    prec_mz <- unlist(lapply(strsplit(names(assembly), split="_"), "[", 5))
    prec_mz <- as.numeric(prec_mz)
    ## get m/z values and corresponding intensities
    mz <- lapply(assembly, function(x) x[,1])
    int <- lapply(assembly, function(x) x[,2])
    ## create Spectrum2 objects with the corresponding information
    lapply(1:length(assembly), function(x) new("Spectrum2", 
        rt=rt[x], precursorMz=prec_mz[x], mz=mz[[x]], intensity=int[[x]]))
}

## neg
## kernel 
spl_MK_cid_neg <- construct_Spectrum2(assembly_MK_cid_neg)
spl_MK_hcd_neg <- construct_Spectrum2(assembly_MK_hcd_neg)

## leaf 
spl_ML_cid_neg <- construct_Spectrum2(assembly_ML_cid_neg)
spl_ML_hcd_neg <- construct_Spectrum2(assembly_ML_hcd_neg)

## sweet maize
spl_swM_cid_neg <- construct_Spectrum2(assembly_swM_cid_neg)
spl_swM_hcd_neg <- construct_Spectrum2(assembly_swM_hcd_neg)

## pos
## kernel 
spl_MK_cid_pos <- construct_Spectrum2(assembly_MK_cid_pos)
spN_MK_hcd_pos <- construct_Spectrum2(assembly_MK_hcd_pos)

## leaf 
spl_ML_cid_pos <- construct_Spectrum2(assembly_ML_cid_pos)
spN_ML_hcd_pos <- construct_Spectrum2(assembly_ML_hcd_pos)

## sweet maize
spl_swM_cid_pos <- construct_Spectrum2(assembly_swM_cid_pos)
spl_swM_hcd_pos <- construct_Spectrum2(assembly_swM_hcd_pos)

## 7)
## create Spectra object from the list of Spectrum2 objects

#' @name create_Spectra
#' @description The function `create_Spectra` creates from a list of `Spectrum2`
#' objects a `Spectra` object. 
#' @param spN_l `list` of `Spectrum2` objects
#' @details The function `create_Spectra` creates the `Spectra` object from 
#' a list of `Spectrum2` objects. It binds a `DataFrame` object 
#' as `elementMetadata` with the columns `"precursorMz"` containing 
#' the m/z values of the precursor ion, `"rt"` containing the retention time 
#' and `"show"` set to `TRUE` for all `Spectrum2` entries. 
#' @usage create_Spectra(spN_l)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @return `Spectra` object
#' @examples
#' create_Spectra(spN)
create_Spectra <- function(spN) {
    MSnbase::Spectra(spN, 
        elementMetadata=S4Vectors::DataFrame(
            precursorMz = unlist(lapply(spN, function(x) x@precursorMz)),
            rt = unlist(lapply(spN, function(x) x@rt)), 
            show=rep(TRUE, length(spN)))) 
}

## neg 
spectra_MK_cid_neg <- create_Spectra(spl_MK_cid_neg)
spectra_MK_hcd_neg <- create_Spectra(spN_MK_hcd_l_neg)
spectra_ML_cid_neg <- create_Spectra(spN_ML_cid_l_neg)
spectra_ML_hcd_neg <- create_Spectra(spN_ML_hcd_l_neg)
spectra_swM_cid_neg <- create_Spectra(spN_swM_cid_l_neg)
spectra_swM_hcd_neg <- create_Spectra(spN_swM_hcd_l_neg)

save(spl_MK_cid_neg, spectra_MK_cid_neg, file = "ms2_spectra_MK_cid_neg.RData")
save(spl_MK_hcd_neg, spectra_MK_hcd_neg, file = "ms2_spectra_MK_hcd_neg.RData")
save(spl_ML_cid_neg, spectra_ML_cid_neg, file = "ms2_spectra_ML_cid_neg.RData")
save(spl_ML_hcd_neg, spectra_ML_hcd_neg, file = "ms2_spectra_ML_hcd_neg.RData")
save(spl_swM_cid_neg, spectra_swM_cid_neg, file = "ms2_spectra_swM_cid_neg.RData")
save(spl_swM_hcd_neg, spectra_swM_hcd_neg, file = "ms2_spectra_swM_hcd_neg.RData")

## pos
spectra_MK_cid_pos <- create_Spectra(spN_MK_cid_l_pos)
spectra_MK_hcd_pos <- create_Spectra(spN_MK_hcd_l_pos)
spectra_ML_cid_pos <- create_Spectra(spN_ML_cid_l_pos)
spectra_ML_hcd_pos <- create_Spectra(spN_ML_hcd_l_pos)
spectra_swM_cid_pos <- create_Spectra(spN_swM_cid_l_pos)
spectra_swM_hcd_pos <- create_Spectra(spN_swM_hcd_l_pos)

save(spl_MK_cid_pos, spectra_MK_cid_pos, file = "ms2_spectra_MK_cid_pos.RData")
save(spl_MK_hcd_pos, spectra_MK_hcd_pos, file = "ms2_spectra_MK_hcd_pos.RData")
save(spl_ML_cid_pos, spectra_ML_cid_pos, file = "ms2_spectra_ML_cid_pos.RData")
save(spl_ML_hcd_pos, spectra_ML_hcd_pos, file = "ms2_spectra_ML_hcd_pos.RData")
save(spl_swM_cid_pos, spectra_swM_cid_pos, file = "ms2_spectra_swM_cid_pos.RData")
save(spl_swM_hcd_pos, spectra_swM_hcd_pos, file = "ms2_spectra_swM_hcd_pos.RData")

## 8) 
## calculate pair-wise similarities using the normalizeddotproduct
## neg
library(MetCirc)
similarityMat_MK_cid_neg <- compare_Spectra(spectra_MK_cid_neg, fun=normalizeddotproduct, binSize=0.01)
similarityMat_MK_hcd_neg <- compare_Spectra(spectra_MK_hcd_neg, fun=normalizeddotproduct, binSize=0.01)
similarityMat_ML_cid_neg <- compare_Spectra(spectra_ML_cid_neg, fun=normalizeddotproduct, binSize=0.01)
similarityMat_ML_hcd_neg <- compare_Spectra(spectra_ML_hcd_neg, fun=normalizeddotproduct, binSize=0.01)
similarityMat_swM_cid_neg <- compare_Spectra(spectra_swM_cid_neg, fun=normalizeddotproduct, binSize=0.01)
similarityMat_swM_hcd_neg <- compare_Spectra(spectra_swM_hcd_neg, fun=normalizeddotproduct, binSize=0.01)

save(similarityMat_MK_cid_neg, file = "ms2_similarityMat_MK_cid_neg.RData")
save(similarityMat_MK_hcd_neg, file = "ms2_similarityMat_MK_hcd_neg.RData")
save(similarityMat_ML_cid_neg, file = "ms2_similarityMat_ML_cid_neg.RData")
save(similarityMat_ML_hcd_neg, file = "ms2_similarityMat_ML_hcd_neg.RData")
save(similarityMat_swM_cid_neg, file = "ms2_similarityMat_swM_cid_neg.RData")
save(similarityMat_swM_hcd_neg, file = "ms2_similarityMat_swM_hcd_neg.RData")

## pos 
similarityMat_MK_cid_pos <- compare_Spectra(spectra_MK_cid_pos, fun=normalizeddotproduct, binSize=0.01)
similarityMat_MK_hcd_pos <- compare_Spectra(spectra_MK_hcd_pos, fun=normalizeddotproduct, binSize=0.01)
similarityMat_ML_cid_pos <- compare_Spectra(spectra_ML_cid_pos, fun=normalizeddotproduct, binSize=0.01)
similarityMat_ML_hcd_pos <- compare_Spectra(spectra_ML_hcd_pos, fun=normalizeddotproduct, binSize=0.01)
similarityMat_swM_cid_pos <- compare_Spectra(spectra_swM_cid_pos, fun=normalizeddotproduct, binSize=0.01)
similarityMat_swM_hcd_pos <- compare_Spectra(spectra_swM_hcd_pos, fun=normalizeddotproduct, binSize=0.01)

save(similarityMat_MK_cid_pos, file = "ms2_similarityMat_MK_cid_pos.RData")
save(similarityMat_MK_hcd_pos, file = "ms2_similarityMat_MK_hcd_pos.RData")
save(similarityMat_ML_cid_pos, file = "ms2_similarityMat_ML_cid_pos.RData")
save(similarityMat_ML_hcd_pos, file = "ms2_similarityMat_ML_hcd_pos.RData")
save(similarityMat_swM_cid_pos, file = "ms2_similarityMat_swM_cid_pos.RData")
save(similarityMat_swM_hcd_pos, file = "ms2_similarityMat_swM_hcd_pos.RData")

## 9) 
## visualize network 
library(igraph)

#' @name get_network_plot
#' @description 
#' @param similarityMat `matrix`, square matrix with similarity values
#' @param assembly 
#' @param threshold `numeric`,  set edges below `threshold` to 0
#' @param file `character
#' @details Set `threshold` to `<= 0` if no thresholding should happen
#' @usage get_network_plot(similarityMat, assembly, threshold = 0.05, file)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @return `igraph` object, will save the network plot to `file`
#' @examples 
#' get_network_plot(similarityMat, assembly, threshold = 0.05, file)
get_network_plot <- function(similarityMat, assembly, threshold = 0.05, file) {
    colnames(similarityMat) <- rownames(similarityMat) <- names(assembly)
    similarityMat[similarityMat < threshold] <- 0
    net <- graph_from_adjacency_matrix(similarityMat, weighted = TRUE,
                                                mode = "undirected", diag = F)
    ## plot the network, save to file
    pdf(file)
    plot(net, vertex.label = NA, vertex.size = 0.01)
    dev.off()
    return(net)
}

## neg
get_network_plot(similarityMat_MK_cid_neg, assembly_MK_cid_neg, file = "network_assembly_MK_cid_neg.pdf")
get_network_plot(similarityMat_MK_hcd_neg, assembly_MK_hcd_neg, file = "network_assembly_MK_hcd_neg.pdf")
get_network_plot(similarityMat_ML_cid_neg, assembly_ML_cid_neg, file = "network_assembly_ML_cid_neg.pdf")
get_network_plot(similarityMat_ML_hcd_neg, assembly_ML_hcd_neg, file = "network_assembly_ML_hcd_neg.pdf")
get_network_plot(similarityMat_swM_cid_neg, assembly_swM_cid_neg, file = "network_assembly_swM_cid_neg.pdf")
get_network_plot(similarityMat_swM_hcd_neg, assembly_swM_hcd_neg, file = "network_assembly_swM_hcd_neg.pdf")

## pos
get_network_plot(similarityMat_MK_cid_pos, assembly_MK_cid_pos, file = "network_assembly_MK_cid_pos.pdf")
get_network_plot(similarityMat_MK_hcd_pos, assembly_MK_hcd_pos, file = "network_assembly_MK_hcd_pos.pdf")
get_network_plot(similarityMat_ML_cid_pos, assembly_ML_cid_pos, file = "network_assembly_ML_cid_pos.pdf")
get_network_plot(similarityMat_ML_hcd_pos, assembly_ML_hcd_pos, file = "network_assembly_ML_hcd_pos.pdf")
get_network_plot(similarityMat_swM_cid_pos, assembly_swM_cid_pos, file = "network_assembly_swM_cid_pos.pdf")
get_network_plot(similarityMat_swM_hcd_pos, assembly_swM_hcd_pos, file = "network_assembly_swM_hcd_pos.pdf")

