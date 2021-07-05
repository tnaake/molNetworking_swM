## source the functions file
source("~/Projects/molNetworking_swM/01_ms2_similarity_functions.R")

## use other function from MetCirc (>v1.15.0) and MSnbase (>v2.11.4) and 
## igraph (>v1.2.4.1)for visualization 
library(MSnbase)
library(igraph)

################################################################################
################################## workflow ####################################
################################################################################
## Maize MS2: create a similarity network
setwd("~/Projects/molNetworking_swM/data_MS2/neg")

## 1) 
################################### negative ###################################
## remove manually the first four rows, and the columns "Post curation result",
## "Fill %", "MS/MS assigned", "Reference RT", "Reference m/z", "Formula",
## "Ontology", "INCHIKEY", "SMILES", "MSI level", "Comment", 
## "Mannualy modified", "Isotope tracking parent ID", 
## "Isotope tracking weight number", "Total score", "RT similarity", 
## "Dot product", "Reverse dot product", "Fragment presence %", "S/N average",
## "Spectrum reference file name", "MS1 isotopic spectrum", "MS/MS spectrum"
##
## save the file to PeakID_0_[pos/neg]_cut.txt

## load file that contains alignment information
aln_neg <- read.table("PeakID_0_neg_cut.txt", sep="\t", 
    header=TRUE, stringsAsFactors=FALSE, quote='"')
## truncate file
## mixMK: maize kernel, mixML: maize leaf, QC: sweet maize

## remove qcMK and qcML columns
aln_neg <- aln_neg[, !colnames(aln_neg) %in% c("qcMK.ddms2.CID30", 
        "qcMK.ddms2.HCD40", "qcML.ddms2.CID30", "qcML.ddms2.HCD40")]

## keep only these rows that have alignment information
cols <- colnames(aln_neg)
inds_remove <- apply(aln_neg[, cols[-c(1:5)]], 1, function(x) all(x == "-2"))
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
aln_pos <- read.table("PeakID_0_pos_cut.txt", sep="\t", 
    header=TRUE, stringsAsFactors=FALSE, quote='"')
## truncate file
## MK: maize kernel, ML: maize, QCsm: sweet maize

## remove qcMK and qcML columns
aln_pos <- aln_pos[, !colnames(aln_pos) %in% 
    c("qcMK.ddms2.CID30", "qcMK.ddms2.HCD40")]

## keep only these rows that have alignment information
cols <- colnames(aln_pos)
inds_remove <- apply(aln_pos[, cols[-c(1:5)]], 1, function(x) all(x == "-2"))
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
MK_hcd30_neg <- MK_hcd30_neg[MK_hcd30_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.HCD30"],]
MK_hcd40_neg <- MK_hcd40_neg[MK_hcd40_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.HCD40"],]
MK_hcd50_neg <- MK_hcd50_neg[MK_hcd50_neg[, "PeakID"] %in% aln_neg[, "mixMK.ddms2.HCD50"],]

## leaf 
ML_cid30_neg <- ML_cid30_neg[ML_cid30_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.CID30"],]
ML_cid40_neg <- ML_cid40_neg[ML_cid40_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.CID40"],]
ML_hcd30_neg <- ML_hcd30_neg[ML_hcd30_neg[, "PeakID"] %in% aln_neg[, "mixML.ddms2.HCD30"],]
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
MK_hcd30_pos <- MK_hcd30_pos[MK_hcd30_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.HCD30"],]
MK_hcd40_pos <- MK_hcd40_pos[MK_hcd40_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.HCD40"],]
MK_hcd50_pos <- MK_hcd50_pos[MK_hcd50_pos[, "PeakID"] %in% aln_pos[, "mixMK.ddms2.HCD50"],]

## leaf 
ML_cid30_pos <- ML_cid30_pos[ML_cid30_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.CID30"],]
ML_cid40_pos <- ML_cid40_pos[ML_cid40_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.CID40"],]
ML_hcd30_pos <- ML_hcd30_pos[ML_hcd30_pos[, "PeakID"] %in% aln_pos[, "mixML.pos.ddms2.HCD30"],]
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
    aln=aln_neg, cols_aln=c("mixMK.ddms2.CID30", "mixMK.ddms2.CID40"),
    sample_name="MK_CID")
assembly_MK_hcd_neg <- assemblySpectra(spectra=list(MK_hcd30_neg, MK_hcd40_neg, MK_hcd50_neg), 
    aln=aln_neg, cols_aln=c("mixMK.ddms2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD")
assembly_ML_cid_neg <- assemblySpectra(spectra=list(ML_cid30_neg, ML_cid40_neg), 
    aln=aln_neg, cols_aln=c("mixML.ddms2.CID30", "mixML.ddms2.CID40"),
    sample_name="ML_CID")
assembly_ML_hcd_neg <- assemblySpectra(spectra=list(ML_hcd30_neg, ML_hcd40_neg, ML_hcd50_neg), 
    aln=aln_neg, cols_aln=c("mixML.ddms2.HCD30", "mixML.ddms2.HCD40", "mixML.ddms2.HCD50"),
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

## load
load("ms2_assembly_MK_cid_neg.RData")
load("ms2_assembly_MK_hcd_neg.RData")
load("ms2_assembly_ML_cid_neg.RData")
load("ms2_assembly_ML_hcd_neg.RData")
load("ms2_assembly_swM_cid_neg.RData")
load("ms2_assembly_swM_hcd_neg.RData")


## pos
assembly_MK_cid_pos <- assemblySpectra(spectra=list(MK_cid30_pos, MK_cid40_pos), 
    aln=aln_pos, cols_aln=c("mixMK.ddms2.CID30", "mixMK.ddms2.CID40"),
    sample_name="MK_CID")
assembly_MK_hcd_pos <- assemblySpectra(spectra=list(MK_hcd30_pos, MK_hcd40_pos, MK_hcd50_pos), 
    aln=aln_pos, cols_aln=c("mixMK.ddms2.HCD30", "mixMK.ddms2.HCD40", "mixMK.ddms2.HCD50"),
    sample_name="MK_HCD")
assembly_ML_cid_pos <- assemblySpectra(spectra=list(ML_cid30_pos, ML_cid40_pos),
    aln=aln_pos, cols_aln=c("mixML.pos.ddms2.CID30", "mixML.pos.ddms2.CID40"),
    sample_name="ML_CID")
assembly_ML_hcd_pos <- assemblySpectra(spectra=list(ML_hcd30_pos, ML_hcd40_pos, ML_hcd50_pos),
    aln=aln_pos, cols_aln=c("mixML.pos.ddms2.HCD30", "mixML.pos.ddms2.HCD40", "mixML.pos.ddms2.HCD50"),
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

## load
load("ms2_assembly_MK_cid_pos.RData")
load("ms2_assembly_MK_hcd_pos.RData")
load("ms2_assembly_ML_cid_pos.RData")
load("ms2_assembly_ML_hcd_pos.RData")
load("ms2_assembly_swM_cid_pos.RData")
load("ms2_assembly_swM_hcd_pos.RData")

# 6) collapse the spectra if the individual assemblies have close rt and 
## precursor m/z values (assume they are identical)
## neg
assembly_collapsed_MK_cid_neg <- collapse_assembly(assembly_MK_cid_neg, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_MK_hcd_neg <- collapse_assembly(assembly_MK_hcd_neg, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_ML_cid_neg <- collapse_assembly(assembly_ML_cid_neg, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_ML_hcd_neg <- collapse_assembly(assembly_ML_hcd_neg, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_swM_cid_neg <- collapse_assembly(assembly_swM_cid_neg, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_swM_hcd_neg <- collapse_assembly(assembly_swM_hcd_neg, rt_dev = 0.1, mz_dev = 0.01)

## save
setwd("/home/thomas/Projects/molNetworking_swM/results_MS2/")
save(assembly_collapsed_MK_cid_neg, file="ms2_assembly_collapsed_MK_cid_neg.RData")
save(assembly_collapsed_MK_hcd_neg, file="ms2_assembly_collapsed_MK_hcd_neg.RData")
save(assembly_collapsed_ML_cid_neg, file="ms2_assembly_collapsed_ML_cid_neg.RData")
save(assembly_collapsed_ML_hcd_neg, file="ms2_assembly_collapsed_ML_hcd_neg.RData")
save(assembly_collapsed_swM_cid_neg, file="ms2_assembly_collapsed_swM_cid_neg.RData")
save(assembly_collapsed_swM_hcd_neg, file="ms2_assembly_collapsed_swM_hcd_neg.RData")

## load
load("ms2_assembly_collapsed_MK_cid_neg.RData")
load("ms2_assembly_collapsed_MK_hcd_neg.RData")
load("ms2_assembly_collapsed_ML_cid_neg.RData")
load("ms2_assembly_collapsed_ML_hcd_neg.RData")
load("ms2_assembly_collapsed_swM_cid_neg.RData")
load("ms2_assembly_collapsed_swM_hcd_neg.RData")

## pos
assembly_collapsed_MK_cid_pos <- collapse_assembly(assembly_MK_cid_pos, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_MK_hcd_pos <- collapse_assembly(assembly_MK_hcd_pos, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_ML_cid_pos <- collapse_assembly(assembly_ML_cid_pos, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_ML_hcd_pos <- collapse_assembly(assembly_ML_hcd_pos, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_swM_cid_pos <- collapse_assembly(assembly_swM_cid_pos, rt_dev = 0.1, mz_dev = 0.01)
assembly_collapsed_swM_hcd_pos <- collapse_assembly(assembly_swM_hcd_pos, rt_dev = 0.1, mz_dev = 0.01)

## save
save(assembly_collapsed_MK_cid_pos, file="ms2_assembly_collapsed_MK_cid_pos.RData")
save(assembly_collapsed_MK_hcd_pos, file="ms2_assembly_collapsed_MK_hcd_pos.RData")
save(assembly_collapsed_ML_cid_pos, file="ms2_assembly_collapsed_ML_cid_pos.RData")
save(assembly_collapsed_ML_hcd_pos, file="ms2_assembly_collapsed_ML_hcd_pos.RData")
save(assembly_collapsed_swM_cid_pos, file="ms2_assembly_collapsed_swM_cid_pos.RData")
save(assembly_collapsed_swM_hcd_pos, file="ms2_assembly_collapsed_swM_hcd_pos.RData")

## load
load("ms2_assembly_collapsed_MK_cid_pos.RData")
load("ms2_assembly_collapsed_MK_hcd_pos.RData")
load("ms2_assembly_collapsed_ML_cid_pos.RData")
load("ms2_assembly_collapsed_ML_hcd_pos.RData")
load("ms2_assembly_collapsed_swM_cid_pos.RData")
load("ms2_assembly_collapsed_swM_hcd_pos.RData")

## 7) 
## convert final deconvoluted spectra to Spectra objects from Spectra
## load the Spectra package
library(Spectra)

## neg 
spectra_MK_cid_neg <- create_Spectra(assembly_collapsed_MK_cid_neg)
spectra_MK_hcd_neg <- create_Spectra(assembly_collapsed_MK_hcd_neg)
spectra_ML_cid_neg <- create_Spectra(assembly_collapsed_ML_cid_neg)
spectra_ML_hcd_neg <- create_Spectra(assembly_collapsed_ML_hcd_neg)
spectra_swM_cid_neg <- create_Spectra(assembly_collapsed_swM_cid_neg)
spectra_swM_hcd_neg <- create_Spectra(assembly_collapsed_swM_hcd_neg)

save(spectra_MK_cid_neg, file = "ms2_spectra_MK_cid_neg.RData")
save(spectra_MK_hcd_neg, file = "ms2_spectra_MK_hcd_neg.RData")
save(spectra_ML_cid_neg, file = "ms2_spectra_ML_cid_neg.RData")
save(spectra_ML_hcd_neg, file = "ms2_spectra_ML_hcd_neg.RData")
save(spectra_swM_cid_neg, file = "ms2_spectra_swM_cid_neg.RData")
save(spectra_swM_hcd_neg, file = "ms2_spectra_swM_hcd_neg.RData")
load("ms2_spectra_MK_cid_neg.RData")
load("ms2_spectra_MK_hcd_neg.RData")
load("ms2_spectra_ML_cid_neg.RData")
load("ms2_spectra_ML_hcd_neg.RData")
load("ms2_spectra_swM_cid_neg.RData")
load("ms2_spectra_swM_hcd_neg.RData")

## pos
spectra_MK_cid_pos <- create_Spectra(assembly_collapsed_MK_cid_pos)
spectra_MK_hcd_pos <- create_Spectra(assembly_collapsed_MK_hcd_pos)
spectra_ML_cid_pos <- create_Spectra(assembly_collapsed_ML_cid_pos)
spectra_ML_hcd_pos <- create_Spectra(assembly_collapsed_ML_hcd_pos)
spectra_swM_cid_pos <- create_Spectra(assembly_collapsed_swM_cid_pos)
spectra_swM_hcd_pos <- create_Spectra(assembly_collapsed_swM_hcd_pos)

save(spectra_MK_cid_pos, file = "ms2_spectra_MK_cid_pos.RData")
save(spectra_MK_hcd_pos, file = "ms2_spectra_MK_hcd_pos.RData")
save(spectra_ML_cid_pos, file = "ms2_spectra_ML_cid_pos.RData")
save(spectra_ML_hcd_pos, file = "ms2_spectra_ML_hcd_pos.RData")
save(spectra_swM_cid_pos, file = "ms2_spectra_swM_cid_pos.RData")
save(spectra_swM_hcd_pos, file = "ms2_spectra_swM_hcd_pos.RData")

load("ms2_spectra_MK_cid_pos.RData")
load("ms2_spectra_MK_hcd_pos.RData")
load("ms2_spectra_ML_cid_pos.RData")
load("ms2_spectra_ML_hcd_pos.RData")
load("ms2_spectra_swM_cid_pos.RData")
load("ms2_spectra_swM_hcd_pos.RData")


## 8) 
## calculate pair-wise similarities using the normalizeddotproduct
## neg
##library(MetCirc)
similarityMat_MK_cid_neg <- compareSpectra(spectra_MK_cid_pos, ppm = 20, FUN = ndotproduct())
similarityMat_MK_hcd_neg <- compareSpectra(spectra_MK_hcd_neg, ppm = 20, FUN = ndotproduct())
similarityMat_ML_cid_neg <- compareSpectra(spectra_ML_cid_neg, ppm = 20, FUN = ndotproduct())
similarityMat_ML_hcd_neg <- compareSpectra(spectra_ML_hcd_neg, ppm = 20, FUN = ndotproduct())
similarityMat_swM_cid_neg <- compareSpectra(spectra_swM_cid_neg, ppm = 20, FUN = ndotproduct())
similarityMat_swM_hcd_neg <- compareSpectra(spectra_swM_hcd_neg, ppm = 20, FUN = ndotproduct())

# similarityMat_MK_cid_neg <- compare_Spectra(spectra_MK_cid_neg, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_MK_hcd_neg <- compare_Spectra(spectra_MK_hcd_neg, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_ML_cid_neg <- compare_Spectra(spectra_ML_cid_neg, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_ML_hcd_neg <- compare_Spectra(spectra_ML_hcd_neg, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_swM_cid_neg <- compare_Spectra(spectra_swM_cid_neg, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_swM_hcd_neg <- compare_Spectra(spectra_swM_hcd_neg, fun=normalizeddotproduct, binSize=0.01)

## save
save(similarityMat_MK_cid_neg, file = "ms2_similarityMat_MK_cid_neg.RData")
save(similarityMat_MK_hcd_neg, file = "ms2_similarityMat_MK_hcd_neg.RData")
save(similarityMat_ML_cid_neg, file = "ms2_similarityMat_ML_cid_neg.RData")
save(similarityMat_ML_hcd_neg, file = "ms2_similarityMat_ML_hcd_neg.RData")
save(similarityMat_swM_cid_neg, file = "ms2_similarityMat_swM_cid_neg.RData")
save(similarityMat_swM_hcd_neg, file = "ms2_similarityMat_swM_hcd_neg.RData")

## load
load("ms2_similarityMat_MK_cid_neg.RData")
load("ms2_similarityMat_MK_hcd_neg.RData")
load("ms2_similarityMat_ML_cid_neg.RData")
load("ms2_similarityMat_ML_hcd_neg.RData")
load("ms2_similarityMat_swM_cid_neg.RData")
load("ms2_similarityMat_swM_hcd_neg.RData")


## pos 
similarityMat_MK_cid_pos <- compareSpectra(spectra_MK_cid_pos, ppm = 20, FUN = ndotproduct())
similarityMat_MK_hcd_pos <- compareSpectra(spectra_MK_hcd_pos, ppm = 20, FUN = ndotproduct())
similarityMat_ML_cid_pos <- compareSpectra(spectra_ML_cid_pos, ppm = 20, FUN = ndotproduct())
similarityMat_ML_hcd_pos <- compareSpectra(spectra_ML_hcd_pos, ppm = 20, FUN = ndotproduct())
similarityMat_swM_cid_pos <- compareSpectra(spectra_swM_cid_pos, ppm = 20, FUN = ndotproduct())
similarityMat_swM_hcd_pos <- compareSpectra(spectra_swM_hcd_pos, ppm = 20, FUN = ndotproduct())

# similarityMat_MK_cid_pos <- compare_Spectra(spectra_MK_cid_pos, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_MK_hcd_pos <- compare_Spectra(spectra_MK_hcd_pos, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_ML_cid_pos <- compare_Spectra(spectra_ML_cid_pos, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_ML_hcd_pos <- compare_Spectra(spectra_ML_hcd_pos, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_swM_cid_pos <- compare_Spectra(spectra_swM_cid_pos, fun=normalizeddotproduct, binSize=0.01)
# similarityMat_swM_hcd_pos <- compare_Spectra(spectra_swM_hcd_pos, fun=normalizeddotproduct, binSize=0.01)

## save
save(similarityMat_MK_cid_pos, file = "ms2_similarityMat_MK_cid_pos.RData")
save(similarityMat_MK_hcd_pos, file = "ms2_similarityMat_MK_hcd_pos.RData")
save(similarityMat_ML_cid_pos, file = "ms2_similarityMat_ML_cid_pos.RData")
save(similarityMat_ML_hcd_pos, file = "ms2_similarityMat_ML_hcd_pos.RData")
save(similarityMat_swM_cid_pos, file = "ms2_similarityMat_swM_cid_pos.RData")
save(similarityMat_swM_hcd_pos, file = "ms2_similarityMat_swM_hcd_pos.RData")

## load
load("ms2_similarityMat_MK_cid_pos.RData")
load("ms2_similarityMat_MK_hcd_pos.RData")
load("ms2_similarityMat_ML_cid_pos.RData")
load("ms2_similarityMat_ML_hcd_pos.RData")
load("ms2_similarityMat_swM_cid_pos.RData")
load("ms2_similarityMat_swM_hcd_pos.RData")

## 9) 
## visualize network 
library(igraph)
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

