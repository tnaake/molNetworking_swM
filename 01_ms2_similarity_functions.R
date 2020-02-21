################################################################################
################################## functions ###################################
################################################################################

## functions to create MS/MS features from corresponding files that were 
## collected under different collision energies

#' @name assemblySpectra
#'
#' @title Deconvolute final spectra from different collision energies
#'
#' @description 
#' `assemblySpectra` takes files from different collision energies
#' as input and assembles a final deconvoluted spectra from it, taking the 
#' highest intensity values from corresponding peaks across the different 
#' collision energies.
#'
#' @param spectra list of matrices, each object is from MS-DIAL and contains
#' information on the MS/MS spectrum
#'
#' @param aln matrix, alignment file from MS-DIAL
#'
#' @param cols_aln character, colnames in aln that correspond to the objects in
#' spectra
#'
#' @param sample_name character, add the name to the names of the returned list
#' (Alignment.ID, retention time and m/z of the precursor, taken from aln, will
#' be automatically added to the name)
#'
#' @param tol numeric, tolerance parameter
#'
#' @details `assemblySpectra` will first bin within each spectra (su)
#'
#' @return list of assembled spectra. Each spectra is a matrix that in the 
#' first column contains m/z values and in the second column contains the 
#' corresponding intensities
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
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
                sp <- spectra[[x]][spectra[[x]][, "PeakID"] == inds_ev[[x]], ]
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
                    binSpectra(x, tol=tol, fun = "sum")   
                }
            })
            
            ## deconvolute final spectra, bin across spectra
            deconv <- deconvolute(i_msms_bin, tol = tol)
            
        } else deconv <- NULL
        
        return(deconv)
    })
    
    ## construct names for res
    names(res) <- paste(sample_name, aln[, "Alignment.ID"], 
        aln[, "Average.Rt.min."], aln[, "Average.Mz"], sep = "_")
    
    ## truncate res that it only contains the entries that are not NULL
    inds_keep <- unlist(lapply(res, function(x) !is.null(x)))
    res <- res[inds_keep]
    
    return(res)
}

#' @name createRefSpectra
#' 
#' @title Create reference spectra from alignment file
#' 
#' @description The alignment object `aln` given by MS-DIAL contains per
#' aligned spectrum a reference spectrum. The function `createRefSpectra`
#' extracts these spectra and returns a list of matrices containing the 
#' m/z values and the corresponding intensities. 
#' 
#' @param aln alignment object (`data.frame`) containing reference spectrum in 
#' the column `MS.MS.spectrum`, as given by MS-DIAL
#' 
#' @details `createRefSpectra` assigns names to the list in the following 
#' format `"ID_id_RT_mz"`, where `id` is the ID of the aligned spectrum taken
#' from the column `"Alignment.ID"`, `RT` is the retention time in minutes
#' taken from the column `"Average.Rt.min."` and `mz` is the m/z value taken 
#' from the column `"Average.Mz"`.
#' 
#' @return list of matrices
#' 
#' @author Thomas Naake <thomasnaake@@googlemail.com> 
#' 
#' @examples
#' createRefSpectra(aln_neg)
#' createRefSpectra(aln_pos)
createRefSpectra <- function(aln) {
    ## retrieve the column "MSMS.spectrum" from each entry in i_spectra,
    ## this column contains information on the peaks and the 
    ## corresponding intensities, strsplit the entries and write them 
    ## to a matrix
    msms <- strsplit(aln[, "MS.MS.spectrum"], split = " ")
    msms <- lapply(msms, function(y) {
        tmp <- strsplit(y, split = ":")
        tmp <- do.call("rbind", tmp)
        mode(tmp) <- "numeric"
        tmp
    })
    names(msms) <- paste("Spectrum_ID", aln[, "Alignment.ID"], 
        aln[, "Average.Rt.min."], aln[, "Average.Mz"], sep = "_")
    return(msms)
}

#' @name binSpectra
#'
#' @title Bin a spectra
#'
#' @description 
#' `binSpectra` combines peaks that are close to each other
#'
#' @param spectra matrix, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#'
#' @param tol numeric, tolerance parameter
#'
#' @param fun character, either "max" or "sum"
#'
#' @details
#' Bin peaks and sum the intensities of corresponding peaks or take the maximum
#' (depending on `fun`). If `fun` is equal to "max" the m/z value of the peak 
#' with the highest intensity is retained togheter with the highest intensity, 
#' if `fun` is equal to "sum" the median of m/z values in the bin are retained 
#' and intensities are added.
#'
#' @return
#' matrix containing in the first column m/z information and in the 
#' second column intensity information
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
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
            spectra <- matrix(c(frag_s[which.max(spectra[, 2])], max(spectra[, 2])), 
                ncol = 2)
        }
    }
    return(spectra)
}

#' @name deconvolute
#' 
#' @title Deconvolute a final spectrum
#' 
#' @description 
#' `deconvolute` takes as input a list of MS/MS spectra. It will bin the m/z 
#' values and take the maximum of the corresponding intensity value.
#'
#' @param spectra list of matrices, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#'
#' @param tol numeric, tolerance parameter for binning
#'
#' @details
#' Takes a list of matrices that contain spectral information and 
#' bins the m/z values. If several peaks are located in the same bin, 
#' the peak with the highest intensity is retained and others are removed
#' from the spectrum.
#'
#' @return matrix containing in the first column m/z information and in the 
#' second column intensity information
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' sp <- list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50)
#' deconvolute(spectra=sp, tol=0.01)
deconvolute <- function(spectra=list(i_msms_hcd30, i_msms_hcd40, i_msms_hcd50), 
    tol = 0.01) {
    
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


#' @name construct_Spectrum2
#'
#' @description
#' The function `construct_Spectrum2` retrieves the relevant information from 
#' the `assembly` (retention time of precursor, m/z of precursor and m/z 
#' and intensities of the fragments) and creates a list of `Spectrum2` objects.
#'
#' @param assembly list of matrices
#'
#' @details
#' The function `construct_Spectrum2` takes as an argument a list of
#' matrices. Each matrix contains in its first column m/z values and in the 
#' second column intensity values.
#' 'construct_Spectrum2` retrieves retention time, precursorMZ, m/z of peaks and 
#' corresponding intensities from the assembled spectrum.
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @return `list` of `Spectrum2` objects
#'
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
    mz <- lapply(assembly, function(x) x[, 1])
    int <- lapply(assembly, function(x) x[, 2])
    ## create Spectrum2 objects with the corresponding information
    lapply(1:length(assembly), function(x) new("Spectrum2", 
        rt=rt[x], precursorMz=prec_mz[x], mz=mz[[x]], intensity=int[[x]]))
}


#' @name create_Spectra
#'
#' @description
#' The function `create_Spectra` creates  a `Spectra` object from a list of 
#' `Spectrum2` objects.
#'
#' @param spl `list` of `Spectrum2` objects
#'
#' @details
#' The function `create_Spectra` creates the `Spectra` object from a list of 
#' `Spectrum2` objects. It binds a `DataFrame` object as `elementMetadata` with 
#' the columns `"precursorMz"` containing the m/z values of the precursor ion, 
#' `"rt"` containing the retention time and `"show"` set to `TRUE` for all 
#' `Spectrum2` entries.
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @return `Spectra` object
#'
#' @examples
#' create_Spectra(spl)
create_Spectra <- function(spl) {
    MSnbase::Spectra(spl, 
        elementMetadata=S4Vectors::DataFrame(
            precursorMz = unlist(lapply(spl, function(x) x@precursorMz)),
            rt = unlist(lapply(spl, function(x) x@rt)), 
            show=rep(TRUE, length(spl)))) 
}


#' @name get_network_plot
#'
#' @description 
#'
#' @param similarityMat `matrix`, square matrix with similarity values
#'
#' @param assembly 
#'
#' @param threshold `numeric`,  set edges below `threshold` to 0
#'
#' @param file `character
#'
#' @details 
#' Set `threshold` to `<= 0` if no thresholding should happen
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @return `igraph` object, will save the network plot to `file`
#'
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


