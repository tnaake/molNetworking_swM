## set working directory
setwd("~/AG-Fernie/Thomas/Data/From Shijuan/maize_pos and neg_ms1/swMaize_K_neg_2019_ms1/mzML/")

## load xcms
library(xcms)
xset_neg <- xcmsSet(file = "./", method="centWave", ppm=20, 
                    snthresh=10, peakwidth=c(5,20), prefilter = c(3, 5000))
xset2_neg <- group(xset_neg, method="density", minfrac=0.5, minsamp=1, bw=5, mzwid=0.025)
xset3_neg <- retcor(xset2_neg, family= "s", plottype= "m", missing=1, extra=1, span=1)
xset4_neg <- group(xset3_neg, method="density", mzwid=0.025, minfrac=0.5, minsamp=1, bw=5)
xset5_neg <- fillPeaks(xset4_neg, method = "chrom")
save("xset_neg", "xset2_neg", "xset3_neg", 
     "xset4_neg", "xset5_neg", file = "./sweetMaize_neg_xcms.RData")

## CAMERA
library(CAMERA)
an_neg <- xsAnnotate(xset5_neg)
anF_neg <- groupFWHM(an_neg, perfwhm = 0.6)
anI_neg <- findIsotopes(anF_neg, mzabs=0.01)
anIC_neg <- groupCorr(anI_neg, cor_eic_th=0.75, graphMethod = "lpc")
anFA_neg <- findAdducts(anIC_neg, polarity="negative")
pl_neg <- getPeaklist(anFA_neg)
save("an_neg", "anF_neg", "anI_neg", "anIC_neg", "anFA_neg", "pl_neg", file = "./sweetMaize_neg_CAMERA.RData")

## set working directory
setwd("~/AG-Fernie/Thomas/Data/From Shijuan/maize_pos and neg_ms1/swMaize_K_pos_2019_ms1/mzML/")

## load xcms
library(xcms)
xset_pos <- xcmsSet(file = "./", method="centWave", ppm = 20, 
                    snthresh = 10, peakwidth = c(5, 20), prefilter = c(3, 5000))
## take a higer bw since some shifts occur
xset2_pos <- group(xset_pos, method = "density", minfrac = 0.5, minsamp = 1, bw = 10, mzwid = 0.025)
xset3_pos <- retcor(xset2_pos, family = "s", plottype = "m", missing = 1, extra = 1, span = 1)
xset4_pos <- group(xset3_pos, method = "density", mzwid = 0.025, minfrac = 0.5, minsamp = 1, bw = 10)
xset5_pos <- fillPeaks(xset4_pos, method = "chrom")
save("xset_pos", "xset2_pos", "xset3_pos", 
     "xset4_pos", "xset5_pos", file = "./sweetMaize_pos_xcms.RData")

## CAMERA
library(CAMERA)
an_pos <- xsAnnotate(xset5_pos)
anF_pos <- groupFWHM(an_pos, perfwhm = 0.6)
anI_pos <- findIsotopes(anF_pos, mzabs=0.01)
anIC_pos <- groupCorr(anI_pos, cor_eic_th=0.75, graphMethod = "lpc")
anFA_pos <- findAdducts(anIC_pos, polarity="positive")
pl_pos <- getPeaklist(anFA_pos)
## with no fill
an_NF_pos <- xsAnnotate(xset4_pos)
anF_NF_pos <- groupFWHM(an_NF_pos, perfwhm = 0.6)
anI_NF_pos <- findIsotopes(anF_NF_pos, mzabs=0.01)
anIC_NF_pos <- groupCorr(anI_NF_pos, cor_eic_th=0.75, graphMethod = "lpc")
anFA_NF_pos <- findAdducts(anIC_NF_pos, polarity="positive")
pl_NF_pos <- getPeaklist(anFA_NF_pos)
save("an_pos", "anF_pos", "anI_pos", "anIC_pos", "anFA_pos", "pl_pos", file = "./sweetMaize_pos_CAMERA.RData")
save("an_NF_pos", "anF_NF_pos", "anI_NF_pos", "anIC_NF_pos", "anFA_NF_pos", "pl_NF_pos", file = "./sweetMaize_pos_CAMERA_noFill.RData")

################################################################################
############################### injection order ################################
################################################################################
##
## remove features from 1.2 min to 13.0 min 
##
cut_rt <- function(peaklist, lower=72, upper=780) {
    peaklist[peaklist[, "rt"] >= lower & peaklist[, "rt"] <= upper,]
} 
## pl_pos_bc <- cut_rt(pl_pos_pos_bc)

pl_neg <- cut_rt(pl_neg)
pl_pos <- cut_rt(pl_pos)
pl_NF_pos <- cut_rt(pl_NF_pos)

## order according to sample order list
setwd("~/AG-Fernie/Thomas/Data/From Shijuan/maize_pos and neg_ms1")
order_inj_neg <- openxlsx::read.xlsx("swMaize_K_neg_2019_ms1/injection_order_neg_byFolder.xlsx", sheet = 1)
order_inj_pos <- openxlsx::read.xlsx("swMaize_K_pos_2019_ms1/injection_order_pos_byFolder.xlsx", sheet = 1)
## remove blanks if any
order_inj_neg <- order_inj_neg[!grepl(order_inj_neg[,"sample_name"], pattern="[B|b]lank"),]
## change all special characters to "."
order_inj_neg[, "sample_name"] <- gsub(pattern = "-", x = order_inj_neg[, "sample_name"], replacement = ".")
order_inj_neg[, "sample_name"] <- gsub(pattern = "_", x = order_inj_neg[, "sample_name"], replacement = ".")
order_inj_neg[, "sample_name"] <- gsub(pattern = "[(]", x = order_inj_neg[, "sample_name"], replacement = ".")
order_inj_neg[, "sample_name"] <- gsub(pattern = "[)]", x = order_inj_neg[, "sample_name"], replacement = ".")
## remove blanks if any
order_inj_pos <- order_inj_pos[!grepl(order_inj_pos[,"sample_name"], pattern="[B|b]lank"),]
## change all special characters to "."
order_inj_pos[, "sample_name"] <- gsub(pattern = "-", x = order_inj_pos[, "sample_name"], replacement = ".")
order_inj_pos[, "sample_name"] <- gsub(pattern = "_", x = order_inj_pos[, "sample_name"], replacement = ".")
order_inj_pos[, "sample_name"] <- gsub(pattern = "[(]", x = order_inj_pos[, "sample_name"], replacement = ".")
order_inj_pos[, "sample_name"] <- gsub(pattern = "[)]", x = order_inj_pos[, "sample_name"], replacement = ".")

## functions to remove columns from matrix or elements from vector
remove_samples_df <- function(peaklist, samples) {
    if (!is.data.frame(peaklist)) stop("peaklist is not a data.frame")
    if (length(samples) > 0) {
        peaklist[, !colnames(peaklist) %in% samples]    
    } else {
        peaklist
    }
    
}

remove_samples_vector <- function(x, samples, names=TRUE) {
    if (length(samples) > 0) {
        if (names) {
            if(is.null(names(x))) stop("x does not contain names")
            x <- x[!names(x) %in% samples]    
        } else {
            x <- x[!x %in% samples]
        }
    }
    x
}

## get batch and peaklist only containing samples
batch_neg <- as.character(sampclass(xset_neg))
batch_pos <- as.character(sampclass(xset_pos))
pl_smp_neg <- pl_neg[, which(colnames(pl_neg) == "C1.2"):which(colnames(pl_neg) == "Y79")]
pl_smp_pos <- pl_pos[, which(colnames(pl_pos) == "C101.1"):which(colnames(pl_pos) == "Y79")]
pl_smp_NF_pos <- pl_NF_pos[, which(colnames(pl_NF_pos) == "C101.1"):which(colnames(pl_NF_pos) == "Y79")]
names(batch_neg) <- colnames(pl_smp_neg)
names(batch_pos) <- colnames(pl_smp_pos)


## check colnames: do they occur in order_inj --> remove if otherwise
colnames(pl_neg)[! colnames(pl_neg) %in% order_inj_neg[, "sample_name"] ]
colnames(pl_pos)[! colnames(pl_pos) %in% order_inj_pos[, "sample_name"] ]
colnames(pl_NF_pos)[! colnames(pl_NF_pos) %in% order_inj_pos[, "sample_name"] ]

## rename "C92.3_1" to "C92.3.1", remove "QC20batch1.pre.[1-7]" and "qc120batch2.pre.[1-8]", "mixsampeNN" from pl_neg, pl_smp_neg 
## and batch_neg
colnames(pl_neg)[which(colnames(pl_neg) == "C92.3_1")] <- "C92.3.1"
colnames(pl_smp_neg)[which(colnames(pl_smp_neg) == "C92.3_1")] <- "C92.3.1"
names(batch_neg)[names(batch_neg) == "C92.3_1"] <- "C92.3.1"

tmp <- c("QC20batch1.pre.1", "QC20batch1.pre.2", "QC20batch1.pre.3", 
    "QC20batch1.pre.4", "QC20batch1.pre.5", "QC20batch1.pre.6", 
    "QC20batch1.pre.7", "qc120batch2.pre.1", "qc120batch2.pre.2", 
    "qc120batch2.pre.3", "qc120batch2.pre.4", "qc120batch2.pre.5",
    "qc120batch2.pre.6", "qc120batch2.pre.7", "qc120batch2.pre.8",
    "mixsample10", "mixsample20", "mixsample21", "mixsample22", "mixsample25",
    "mixsample31", "mixsample33", "mixsample34", "mixsample35", "mixsample36", 
    "mixsample37", "mixsample42", "mixsample43", "mixsample44", "mixsample45")   
pl_neg <- remove_samples_df(peaklist = pl_neg, samples = tmp)
pl_smp_neg <- remove_samples_df(peaklist = pl_smp_neg, samples = tmp)
batch_neg <- remove_samples_vector(x=batch_neg, samples = tmp)

## remove no column from pl_pos, pl_NF_pos and batch_pos
# tmp <-  c("C159.3", "C169.1")
# pl_pos <- remove_samples_df(peaklist = pl_pos, samples = tmp)
# pl_smp_pos <- remove_samples_df(peaklist = pl_smp_pos, samples = tmp)
# pl_NF_pos <- remove_samples_df(peaklist = pl_NF_pos, samples = tmp)
# pl_smp_NF_pos <- remove_samples_df(peaklist = pl_smp_NF_pos, samples = tmp)
# batch_pos <- remove_samples_vector(x = batch_pos, samples = tmp)

order_inj_neg[!order_inj_neg[, "sample_name"] %in% colnames(pl_neg),]
order_inj_pos[!order_inj_pos[, "sample_name"] %in% colnames(pl_pos),]
order_inj_pos[!order_inj_pos[, "sample_name"] %in% colnames(pl_NF_pos),]
## remove nothing from order_inj_neg
## remove nothing from order_inj_pos
order_inj_neg_s <- order_inj_neg[, "sample_name"]
names(order_inj_neg_s) <- order_inj_neg_s
order_inj_pos_s <- order_inj_pos[, "sample_name"]
names(order_inj_pos_s) <- order_inj_pos_s

################################################################################
################# retention time correction for positive mode ##################
################################################################################

## batchCorr
library(devtools)
install_git("https://gitlab.com/CarlBrunius/batchCorr.git")
library(batchCorr)
## create matrix with information on mz and rt
peakIn <- as.matrix(pl_pos[, c("mz", "rt")])
mode(peakIn) <- "numeric"

##create meta data.frame
grp <- order_inj_pos[, "type"]## ifelse(grepl(colnames(pl_smp_pos), pattern = "[q|Q][c|C]"), "QC", "Ref")
inj <- match(colnames(pl_smp_pos), order_inj_pos[, "sample_name"])
meta <- data.frame(batch=batch_pos, grp=grp, inj=inj)

## alignBatches 
alignBat <- alignBatches(peakInfo = peakIn, PeakTabNoFill = t(pl_smp_NF_pos), 
        PeakTabFilled = t(pl_smp_pos), batches = meta$batch, 
        sampleGroups = meta$grp, selectGroup = "QC", rtdiff = 30)
pl_smp_pos_bc <- t(alignBat$PTalign)
dim(pl_smp_pos)
dim(pl_smp_pos_bc) ## does not change dims --> no removed features


## define function pca_plot to create a PCA plot
pca_plot <- function(peaklist, batch, file, text=FALSE) {
    p <- prcomp(t(peaklist))   
    pdf(file)
    plot(p$x[,1], p$x[,2], col=as.numeric(as.factor(batch))+1)
    if (text) text(p$x[,1], p$x[,2], labels = colnames(peaklist), cex = 0.3)
    dev.off()
}

## define function boxplot_plot to create boxplot
boxplot_plot <- function(peaklist, file) {
    pdf(file)
    boxplot(as.matrix(peaklist))
    dev.off()    
}


setwd("~/AG-Fernie/Thomas/Data/From Shijuan/maize_pos and neg_ms1/results_ms1")
pca_plot(pl_smp_pos_bc, batch_pos, file = "pca_peaklist_pos_batchCorr.pdf")
pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_neg_sample_raw.pdf")
pca_plot(pl_smp_pos, batch_pos, "pca_peaklist_pos_sample_raw.pdf")

################################################################################
################################ normalization #################################
################################################################################

## get log2 of intensities and plot PCA
pl_smp_neg <- log2(pl_smp_neg + 1)
pl_smp_pos <- log2(pl_smp_pos + 1)

pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_neg_sample_log.pdf", text=TRUE)
pca_plot(pl_smp_pos, batch_pos, "pca_peaklist_pos_sample_log.pdf", text=TRUE)

## remove QC14, E32, C90.2 from negative since they are outliers
tmp <- c("QC14", "E32", "C90.2")
pl_smp_neg <- pl_smp_neg[, !colnames(pl_smp_neg) %in% tmp]
batch_neg <- batch_neg[!names(batch_neg) %in% tmp]
order_inj_neg <- order_inj_neg[!order_inj_neg[, "sample_name"] %in% tmp, ]
pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_neg_sample_log_withoutOutlier.pdf")

################################################################################
####################### divide by sum of total ion count #######################
################################################################################
sumTIC <- apply(pl_smp_neg, 2, sum)
pl_smp_neg_tic <- sweep(pl_smp_neg, MARGIN=2, STATS=sumTIC, FUN="/")
sumTIC <- apply(pl_smp_pos, 2, sum)
pl_smp_pos_tic <- sweep(pl_smp_pos, MARGIN=2, STATS=sumTIC, FUN="/")

## boxplots
boxplot_plot(pl_smp_neg, "boxplot_peaklist_neg_sample_log_raw.pdf")
boxplot_plot(pl_smp_pos, "boxplot_peaklist_pos_sample_log_raw.pdf")
boxplot_plot(pl_smp_neg_tic, "boxplot_peaklist_neg_sample_log_tic.pdf")
boxplot_plot(pl_smp_pos_tic, "boxplot_peaklist_pos_sample_log_tic.pdf")

## PCA plots
pca_plot(pl_smp_neg_tic, batch_neg, "pca_peaklist_neg_sample_log_tic.pdf")
pca_plot(pl_smp_pos_tic, batch_pos, "pca_peaklist_pos_sample_log_tic.pdf")
#pca_plot(pl_smp_neg_batch_tic, batch_neg, "pca_peaklist_neg_sample_log_batch_tic.pdf")
#pca_plot(pl_smp_pos_batch_tic, batch_pos, "pca_peaklist_pos_sample_log_batch_tic.pdf")

################################################################################
########## MS-Dial (LOWESS normalization tool): lowess normalization ###########
################################################################################

## reorder pl and pl_smp according to order_inj
pl_smp_neg_tic <- pl_smp_neg_tic[, order_inj_neg[, "sample_name"]]
pl_smp_pos_tic <- pl_smp_pos_tic[, order_inj_pos[, "sample_name"]]
batch_neg <- batch_neg[order_inj_neg[, "sample_name"]]
batch_pos <- batch_pos[order_inj_pos[, "sample_name"]]

## convert to a file format that is readable by 
## Name Type Order Feat1 Feat2 Feat3 ...
## char log  num   num   num   num   ...
pl_smp_neg_tic_t <- cbind(Name = order_inj_neg[, "sample_name"], 
    Type = order_inj_neg[, "type"], Order = 1:ncol(pl_smp_neg_tic), t(pl_smp_neg_tic))
pl_smp_neg_tic_t[, "Type"] <- grepl(pl_smp_neg_tic_t[, "Type"], pattern = "QC")
write.table(pl_smp_neg_tic_t, file = "peaklist_neg_log_tic.txt", sep = "\t", 
    dec = ".", row.names = FALSE, quote = FALSE)

pl_smp_pos_tic_t <- cbind(Name = order_inj_pos[, "sample_name"], 
    Type = order_inj_pos[, "type"], Order = 1:ncol(pl_smp_pos_tic), t(pl_smp_pos_tic))
pl_smp_pos_tic_t[, "Type"] <- grepl(pl_smp_pos_tic_t[, "Type"], pattern = "QC")
write.table(pl_smp_pos_tic_t, file = "peaklist_pos_log_tic.txt", sep = "\t", 
    dec = ".", row.names = FALSE, quote = FALSE)
## open in MS-DIAL and normalize


## !! adjust file names of returned files !!
pl_smp_neg_t_lowess <- read.table(file = "peaklist_neg_log_tic_20199271620.txt", 
    sep = "\t", header = TRUE, row.names = 1, quote = "'")
pl_smp_pos_t_lowess <- read.table(file = "peaklist_pos_log_tic_20199271620.txt", 
    sep = "\t", header = TRUE, row.names = 1, quote = "'")
##pl_smp_lowess <- t(pl_smp_t_lowess)
tmp <- c("TYPE", "ORDER")
pl_smp_neg_lowess <- remove_samples_df(pl_smp_neg_t_lowess, tmp)
pl_smp_pos_lowess <- remove_samples_df(pl_smp_pos_t_lowess, tmp)

pl_smp_neg_lowess <- t(pl_smp_neg_lowess)
pl_smp_pos_lowess <- t(pl_smp_pos_lowess)

##pl_smp_lowess <- pl_smp_lowess[!rownames(pl_smp_lowess) %in% c("TYPE", "ORDER"), ]
mode(pl_smp_neg_lowess) <- mode(pl_smp_pos_lowess) <- "numeric"

pca_plot(pl_smp_neg_lowess, batch_neg, "pca_peaklist_neg_sample_log_lowess.pdf", text = TRUE)
pca_plot(pl_smp_pos_lowess, batch_pos, "pca_peaklist_pos_sample_log_lowess.pdf", text = TRUE)
## do not use LOWESS

################################################################################
############################### batch correction ###############################
################################################################################
## batch correction using removeBatchEffect
smp_neg <- data.frame(sample = colnames(pl_smp_neg_tic), colour = batch_neg, day = batch_neg)
smp_pos <- data.frame(sample = colnames(pl_smp_pos_tic), colour = batch_pos, day = batch_pos)
pl_smp_neg_batch <- limma::removeBatchEffect(x = pl_smp_neg_tic, batch = smp_neg[, 3])
pl_smp_pos_batch <- limma::removeBatchEffect(x = pl_smp_pos_tic, batch = smp_pos[, 3])

pca_plot(pl_smp_neg_batch, batch_neg, "pca_peaklist_neg_sample_log_batch.pdf")
pca_plot(pl_smp_pos_batch, batch_pos, "pca_peaklist_pos_sample_log_batch.pdf", text = TRUE)

## outliers in PCA --> samples are all at the end of batch2 --> assign these 
## to batch3 (batchC) and rerun removeBatchEffect
tmp <- c("C261.2", "C28.2", "E14", "C282.1", "C56.2", "C189.2", "C63.1", 
    "C285.2", "C92.3.2", "C170.1", "C110.2", "QC24", "Y119", "C315.1", 
    "C92.3.1", "C126.1", "C143.1", "C73.2", "C289.2")
smp_pos[smp_pos[, "sample"] %in% tmp, "colour"] <- "batchC"
smp_pos[smp_pos[, "sample"] %in% tmp, "day"] <- "batchC"
pl_smp_pos_batch <- limma::removeBatchEffect(x = pl_smp_pos_tic, batch = smp_pos[, 3])
pca_plot(pl_smp_pos_batch, batch_pos, "pca_peaklist_pos_sample_log_batch_rename.pdf", text = TRUE)

boxplot_plot(pl_smp_neg_batch, "boxplot_peaklist_neg_sample_log_batch.pdf")
boxplot_plot(pl_smp_pos_batch, "boxplot_peaklist_pos_sample_log_batch.pdf")


################################################################################
############################## write to a file  ################################
################################################################################

mzRT <- pl_neg[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "isotopes", "adduct", "pcgroup")]
pl_final <- cbind(mzRT, pl_smp_neg_batch)
write.table(pl_final, file = "peaklist_neg_log_tic_batch_final_swM.txt", sep="\t", dec=".", quote=FALSE) 

mzRT <- pl_pos[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "isotopes", "adduct", "pcgroup")]
pl_final <- cbind(mzRT, pl_smp_pos_batch)
write.table(pl_final, file = "peaklist_pos_log_tic_batch_final_swM.txt", sep="\t", dec=".", quote=FALSE) 