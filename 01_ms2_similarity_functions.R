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
assemblySpectra <- function(spectra = list(MK_hcd30, MK_hcd40, MK_hcd50), 
    aln = aln,
    cols_aln = c("mixMK.ddm2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name = "MK_HCD", tol = 0.01) {
    
    
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
                    msms <- strsplit(x[, "MSMS.spectrum"], split = " ")[[1]]  
                    msms <- lapply(msms, function(y) strsplit(y, split = ":"))
                    msms <- matrix(unlist(msms), ncol = 2, byrow = TRUE)
                    mode(msms) <- "numeric"
                    return(msms)
                }
            })
            
            ## bin within each spectra of i_spectra_msms
            ## sum up corresponding intensities from peaks that are  binned
            ## together
            i_msms_bin <- lapply(i_msms, function(x) {
                if (!is.null(x)) {
                    binSpectra(x, tol = tol, fun = "sum")   
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


#' @name collapse_assembly
#' 
#' @title Collapse assembly
#' 
#' @description 
#' 
#' @param assembly
#' @param rt_dev rt deviance
#' @param mz_dev
#' 
#' @details 
#' 
#' @return list of matrices
#' 
#' @author Thomas Naake <thomasnaake@@googlemail.com> 
#' 
#' @examples
#' collapse_assembly
collapse_assembly <- function(assembly, rt_dev, mz_dev) {
    
    names_assembly <- names(assembly)
    names_assembly <- strsplit(names_assembly, split = "_")
    rt <- lapply(names_assembly, "[", 4)
    rt <- as.numeric(unlist(rt))
    mz <- lapply(names_assembly, "[", 5)
    mz <- as.numeric(unlist(mz))
    
    ## iterate through list of features and calculate deviances
    for (i in 1:length(mz)) {
        mz_devs <- abs(mz[i] - mz)
        rt_devs <- abs(rt[i] - rt)
        inds <- which(mz_devs <= mz_dev & rt_devs <= rt_dev)
        
        ## remove i from the vector of indices that are within the boundaries
        inds <- inds[inds != i]
        
        ## combine all the matching assemblies into one
        assembly_comb <- Reduce(rbind, assembly[c(i, inds)])
        assembly[[i]] <- assembly_comb
        
        ## remove the entries of the assemblies that were matched
        assembly[inds] <- lapply(inds, function(x) matrix(nrow = 0, ncol = 2))
    }
    
    ## remove the empty entries
    assembly_nrows <- unlist(lapply(assembly, nrow))
    assembly <- assembly[assembly_nrows != 0]
    
    ## bin the mz values (this will collapse features with mz values < 0.01
    ## together by aggregating via the sum of intensities)
    assembly <- lapply(assembly, function(x) 
        binAssembly(x, tol = 0.01, fun = "sum"))
    
    return(assembly)
}

#' @name binAssembly
#'
#' @title Bin an assembly
#'
#' @description 
#' `binAssembly` combines peaks that are close to each other
#'
#' @param assembly matrix, the first column contains information on the 
#' m/z and the second column contains information on the intensity
#' @param tol numeric, tolerance parameter
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
#' binAssembly(assembly, tol=0.01, fun="sum")
binAssembly <- function(assembly, tol = 0.01, fun=c("max", "sum")) {
    
    if (!is.matrix(assembly)) stop("assembly is not a matrix")
    ## match arguments for fun
    fun <- match.arg(fun)
    
<<<<<<< Updated upstream
    ## if assembly only contains > 2 do the binning, otherwise return assembly
    if (nrow(assembly) > 2) {
        ## order the assembly according to mz
        assembly <- assembly[order(assembly[, 1]), ]
        frag_s <- assembly[, 1]
=======
    ## if spectra only contains > 2 do the binning, otherwise return spectra
    if (nrow(spectra) > 2) {
        
        frag_s <- spectra[,1]
>>>>>>> Stashed changes
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
                assembly_dup <- assembly[inds_dup,]
                
                ## either use max or sum the intensities
                if (fun == "max") {
                    assembly[inds_dup, 1] <- assembly_dup[which.max(assembly_dup[, 2]), 1]
                    assembly[inds_dup, 2] <- max(assembly_dup[, 2])    
                    ## set all except the ones with the highest intensity to NA
                    assembly[inds_dup[-which.max(assembly_dup[, 2])], ] <- NA 
                } 
                if (fun == "sum") {
                    assembly[inds_dup, 1] <- median(assembly_dup[, 1])
                    assembly[inds_dup, 2] <- sum(assembly_dup[, 2])   
                    ## set all except the first one to NA
                    assembly[inds_dup[-1], ] <- NA 
                }
            }
            ## remove NA values
<<<<<<< Updated upstream
            assembly <- assembly[!is.na(assembly[, 1]), ]
        } else {
            assembly <- matrix(c(frag_s[which.max(assembly[, 2])], max(assembly[, 2])), 
                              ncol = 2)
=======
            spectra <- spectra[!is.na(spectra[,1]), ]
            
        } else {
            spectra <- matrix(
                c(frag_s[which.max(spectra[, 2])], max(spectra[, 2])), 
                ncol = 2)
>>>>>>> Stashed changes
        }
    }
    return(assembly)
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
    if (!is.matrix(spectra)) 
        spectra <- matrix(spectra, ncol = 2, byrow = TRUE)
    
    ## bin: take the max from the same corresonding peak
    spectra <- binSpectra(spectra, tol = tol, fun = "max")
    
    ## return
    return(spectra)
}


#' @name create_Spectra
#'
#' @description
#' The function `create_Spectra` retrieves the relevant information from 
#' the `assembly` (retention time of precursor, m/z of precursor and m/z 
#' and intensities of the fragments) and creates a `Spectra` object from this
#' list.
#'
#' @param assembly list of matrices
#'
#' @details
#' The function `create_Spectra` takes as an argument a list of
#' matrices. Each matrix contains in its first column m/z values and in the 
#' second column intensity values.
#' `create_spectra` retrieves retention time, precursorMZ, m/z of peaks and 
#' corresponding intensities from the assembled spectrum.
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @return `Spectra` object
#'
#' @examples
#' create_Spectra(assembly)
create_Spectra <- function(assembly) {
    
    names_a <- names(assembly)
<<<<<<< Updated upstream
    rt <- unlist(lapply(strsplit(names_a, split="_"), "[", 4))  
    rt <- as.numeric(rt)
    ## get precursor m/z
    prec_mz <- unlist(lapply(strsplit(names_a, split="_"), "[", 5))
=======
    rt <- unlist(lapply(strsplit(names_a, split = "_"), "[", 4))  
    rt <- as.numeric(rt)
    ## get precursor m/z
    prec_mz <- unlist(lapply(strsplit(names_a, split = "_"), "[", 5))
>>>>>>> Stashed changes
    prec_mz <- as.numeric(prec_mz)
    
    ## get m/z values and corresponding intensities
    mz <- lapply(assembly, function(x) x[, 1])
    intensity <- lapply(assembly, function(x) x[, 2])
    
    spd <- DataFrame(msLevel = c(2L), rtime = rt, precursorMz = prec_mz)
    spd$mz <- mz
    spd$intensity <- intensity
    
    data <- Spectra::Spectra(spd)
    spd_names <- lapply(strsplit(names_a, "_"), "[", 1:3)
    spd_names <- unlist(lapply(spd_names, paste, collapse = "_"))
    data$spectrum_id <- spd_names
    
    return(data)
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

#' @name map_precursor
#'
#' @title Map precursors of MS2 to MS1
#'
#' @description The function `map_precursor` takes a `Spectra` object
#' containing MS2 information and a peaklist `matrix` with MS1 information and
#' returns those features in the MS1 `matrix` whose features with m/z and
#' retention time pairs match to the precursor ions of the `Spectra` object.
#'
#' @param spectra `Spectra` object containing MS2 information
#'
#' @param peaklist `matrix`, containing the columns `"mz"` and `"rt"`
#'
#' @param ppm `numeric`, m/z tolerance in ppm
#'
#' @param rt_tol `numeric`, retention time tolerance in minutes
#'
#' @details If several features in `peaklist` match to a MS2 spectrum, the
#' function `map_precursor` will return all matching features.
#'
#' @usage map_precursor(spectra, peaklist, ppm = 20, rt_tol = 0.15)
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @return `matrix`
#'
#' @examples
#' map_precursor(spectra, peaklist, ppm = 20, rt_tol = 0.15)
map_precursor <- function(spectra, peaklist, ppm = 20, rt_tol = 0.15, 
    rt_ms1 = c("seconds", "minutes")) {
    
    
    if (rt_ms1 == "seconds") {
        peaklist[, "rt"] <- peaklist[, "rt"] / 60
    }
    
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
        if (nrow(pl) > 0) 
            pl <- cbind(spectra = names(spectra)[x], pl)
        
        ## return   
        return(pl)
    })
    
    ## rbind the list 
    do.call("rbind", pl_map) 
}

