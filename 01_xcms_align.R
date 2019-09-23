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
xset_pos <- xcmsSet(file = "./", method="centWave", ppm=20, 
                snthresh=10, peakwidth=c(5,20), prefilter = c(3, 5000))
xset2_pos <- group(xset_pos, method="density", minfrac=0.5, minsamp=1, bw=5, mzwid=0.025)
xset3_pos <- retcor(xset2_pos, family= "s", plottype= "m", missing=1, extra=1, span=1)
xset4_pos <- group(xset3_pos, method="density", mzwid=0.025, minfrac=0.5, minsamp=1, bw=5)
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
save("an_pos", "anF_pos", "anI_pos", "anIC_pos", "anFA_pos", "pl_pos", file = "./sweetMaize_pos_CAMERA.RData")

################################################################################
################################ normalization #################################
################################################################################

##
## remove features from 1.2 min to 13.0 min 
##
cut_rt <- function(peaklist, lower=72, upper=780) {
    peaklist[peaklist[, "rt"] >= lower & peaklist[, "rt"] <= upper,]
} 
pl_neg <- cut_rt(pl_neg)
pl_pos <- cut_rt(pl_pos)

##
## normalization
##

## get batch
batch_neg <- as.character(sampclass(xset_neg))
batch_pos <- as.character(sampclass(xset_pos))
pl_smp_neg <- pl_neg[, which(colnames(pl_neg) == "C1.2"):which(colnames(pl_neg) == "Y79")]
pl_smp_pos <- pl_pos[, which(colnames(pl_pos) == "C1.2"):which(colnames(pl_pos) == "Y79")]
names(batch_neg) <- colnames(pl_smp_neg)
names(batch_pos) <- colnames(pl_smp_pos)

## define function pca_plot to create a PCA plot
pca_plot <- function(peaklist, batch, file, text=FALSE) {
    p <- prcomp(t(peaklist))   
    pdf(file)
    plot(p$x[,1], p$x[,2], col=as.numeric(as.factor(batch))+1)
    if (text) text(p$x[,1], p$x[,2], labels = colnames(peaklist), cex = 0.3)
    dev.off()
}

setwd("~/AG-Fernie/Thomas/Data/From Shijuan/maize_pos and neg_ms1/results_ms1")
pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_neg_sample_raw.pdf")
pca_plot(pl_smp_pos, batch_pos, "pca_peaklist_pos_sample_raw.pdf")


## get log2 of intensities and plot PCA
pl_smp_neg <- log2(pl_smp_neg + 1)
pl_smp_pos <- log2(pl_smp_pos + 1)

pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_neg_sample_log.pdf", text=TRUE)
pca_plot(pl_smp_pos, batch_pos, "pca_peaklist_pos_sample_log.pdf", text=TRUE)

## remove QC14, E32, C90.2 from negative since they are outliers
pl_smp_neg <- pl_smp_neg[, !colnames(pl_smp_neg) %in% c("QC14", "E32", "C90.2")]
batch_neg <- batch_neg[!names(batch_neg) %in% c("QC14", "E32", "C90.2")]
pca_plot(pl_smp_neg, batch_neg, "pca_peaklist_sample_log_withoutOutlier.pdf")

################################################################################
## MS-Dial (LOWESS normalization tool): lowess normalization 
## order according to sample order list
order_inj <- read.table("../swMaize_K_neg_2019_ms1/mzML/Sweet_kernel-run_sequence-20190918.txt", stringsAsFactors = FALSE)
## remove blanks
order_inj <- order_inj[-grep(order_inj[,1], pattern="[B|b]lank"),]
## change all special characters to "."
order_inj <- gsub(pattern = "-", x = order_inj, replacement = ".")
order_inj <- gsub(pattern = "_", x = order_inj, replacement = ".")
order_inj <- gsub(pattern = "[(]", x = order_inj, replacement = ".")
order_inj <- gsub(pattern = "[)]", x = order_inj, replacement = ".")

## functions to remove columns from matrix or elements from vector
remove_samples_matrix <- function(peaklist, samples) {
    if (!is.matrix(peaklist)) stop("peaklist is not a matrix")
    peaklist[, -which(colnames(peaklist) %in% samples)]
}

remove_samples_vector <- function(x, samples, names=TRUE) {
    if (names) {
        if(is.null(names(x))) stop("x does not contain names")
        x <- x[!names(x) %in% samples]    
    } else {
        x <- x[!x %in% samples]
    }
    x
}

## check colnames: do they occur in order_inj --> remove if otherwise
colnames(pl_smp_neg)[! colnames(pl_smp_neg) %in% order_inj ]
colnames(pl_smp_pos)[! colnames(pl_smp_pos) %in% order_inj ]

## rename "C92.3_1" to "C92.3.1", remove "C159.3" from pl_neg, pl_smp_neg 
## and batch_neg
colnames(pl_neg)[which(colnames(pl_neg) == "C92.3_1")] <- "C92.3.1"
colnames(pl_smp_neg)[which(colnames(pl_smp_neg) == "C92.3_1")] <- "C92.3.1"
tmp <- "C159.3"
pl_neg <- remove_samples_matrix(peaklist = as.matrix(pl_neg), samples = tmp)
pl_smp_neg <- remove_samples_matrix(peaklist = as.matrix(pl_smp_neg), samples = tmp)
batch_neg <- remove_samples_vector(x=batch_neg, samples = tmp)

## remove "C159.3", "QC120batch1.pre.[1-6]", "C169.1" from pl_pos, pl_smp_pos,
## and batch_pos
tmp <-  c("C159.3", "QC120batch1.pre.1", "QC120batch1.pre.2", 
          "QC120batch1.pre.3", "QC120batch1.pre.4", "QC120batch1.pre.5", 
          "QC120batch1.pre.6", "C169.1")
pl_pos <- remove_samples_matrix(peaklist = as.matrix(pl_pos), samples = tmp)
pl_smp_pos <- remove_samples_matrix(peaklist = as.matrix(pl_smp_pos), samples = tmp)
batch_pos <- remove_samples_vector(x = batch_pos, samples = tmp)

order_inj[!order_inj %in% colnames(pl_smp_neg)]
order_inj[!order_inj %in% colnames(pl_smp_pos)]
## remove "C119.3", "QC14", "C193.3", "C90.2", "E32" from order_inj and assign
## to order_inj_neg
## remove "qc120batch2.pre.[1-8]", "C119.3", "C261.1" and "C145.3.2" from 
## order_inj and assign to order_inj_pos
names(order_inj) <- order_inj
order_inj_neg <- remove_samples_vector(order_inj, c("C119.3", "QC14", "C193.3", "C90.2", "E32"), names = FALSE)
order_inj_pos <- remove_samples_vector(order_inj, c("qc120batch2.pre.1", 
    "qc120batch2.pre.2", "qc120batch2.pre.3", "qc120batch2.pre.4", 
    "qc120batch2.pre.5", "qc120batch2.pre.6", "qc120batch2.pre.7", 
    "qc120batch2.pre.8", "C119.3", "C261.1", "C145.3.2"), names = FALSE)

## reorder pl and pl_smp according to order_inj
pl_smp_neg <- pl_smp_neg[, order_inj_neg]
pl_smp_pos <- pl_smp_pos[, order_inj_pos]

## convert to a file format that is readable by 
## Name Type Order Feat1 Feat2 Feat3 ...
## char log  num   num   num   num   ...
pl_smp_neg_t <- cbind(Name=colnames(pl_smp_neg), Type=colnames(pl_smp_neg), Order=1:ncol(pl_smp_neg), t(pl_smp_neg))
pl_smp_neg_t[, "Type"] <- grepl(pl_smp_neg_t[, "Type"], pattern="[Q|q][C|c]")
pl_smp_neg_t[grep(pl_smp_neg_t[, "Name"], pattern="batch"), "Type"] <- "FALSE"
write.table(pl_smp_neg_t, file="peaklist_neg_log_withoutOutlier.txt", sep="\t", dec=".", row.names=FALSE, quote=FALSE)

pl_smp_pos_t <- cbind(Name=colnames(pl_smp_pos), Type=colnames(pl_smp_pos), Order=1:ncol(pl_smp_pos), t(pl_smp_pos))
pl_smp_pos_t[, "Type"] <- grepl(pl_smp_pos_t[, "Type"], pattern="[Q|q][C|c]")
pl_smp_pos_t[grep(pl_smp_pos_t[, "Name"], pattern="batch"), "Type"] <- "FALSE"
write.table(pl_smp_pos_t, file="peaklist_pos_log_withoutOutlier.txt", sep="\t", dec=".", row.names=FALSE, quote=FALSE)
## open in MS-DIAL and normalize

## !! adjust file names of returned files !!
pl_smp_neg_t_lowess <- read.table(file="peaklist_neg_log_withoutOutlier_20199231017.txt", sep="\t", header=TRUE, row.names = 1, quote="'")
pl_smp_pos_t_lowess <- read.table(file="peaklist_pos_log_withoutOutlier_20199231028.txt", sep="\t", header=TRUE, row.names = 1, quote="'")
##pl_smp_lowess <- t(pl_smp_t_lowess)
pl_smp_neg_lowess <- remove_samples_matrix(pl_smp_neg_lowess, c("TYPE", "ORDER"))
pl_smp_pos_lowess <- remove_samples_matrix(pl_smp_pos_lowess, c("TYPE", "ORDER"))

pl_smp_neg_lowess <- t(pl_smp_neg_lowess)
pl_smp_pos_lowess <- t(pl_smp_pos_lowess)

##pl_smp_lowess <- pl_smp_lowess[!rownames(pl_smp_lowess) %in% c("TYPE", "ORDER"), ]
mode(pl_smp_neg_lowess) <- mode(pl_smp_pos_lowess) <- "numeric"

pca_plot(pl_smp_neg_lowess, batch_neg, "pca_peaklist_neg_sample_log_lowess.pdf")
pca_plot(pl_smp_pos_lowess, batch_pos, "pca_peaklist_pos_sample_log_lowess.pdf")
## do not use LOWESS

boxplot_plot <- function(peaklist, file) {
    pdf(file)
    boxplot(as.matrix(peaklist))
    dev.off()    
}
boxplot_plot(pl_smp_neg, "boxplot_peaklist_neg_sample_log.pdf")
boxplot_plot(pl_smp_pos, "boxplot_peaklist_pos_sample_log.pdf")


################################################################################
############################### batch correction ###############################
################################################################################
## batch correction using removeBatchEffect
smp_neg <- data.frame(sample = colnames(pl_smp_neg), colour = batch_neg, day = batch_neg)
smp_pos <- data.frame(sample = colnames(pl_smp_pos), colour = batch_pos, day = batch_pos)
pl_smp_neg_batch <- limma::removeBatchEffect(x = pl_smp_neg, batch = smp_neg[, 3])
pl_smp_pos_batch <- limma::removeBatchEffect(x = pl_smp_pos, batch = smp_pos[, 3])

pca_plot(pl_smp_neg_batch, batch_neg, "pca_peaklist_neg_sample_log_batch.pdf")
pca_plot(pl_smp_pos_batch, batch_pos, "pca_peaklist_pos_sample_log_batch.pdf")

boxplot_plot(pl_smp_neg_batch, "boxplot_peaklist_neg_sample_log_batch.pdf")
boxplot_plot(pl_smp_pos_batch, "boxplot_peaklist_pos_sample_log_batch.pdf")

################################################################################
####################### divide by sum of total ion count #######################
################################################################################
sumTIC <- apply(pl_smp_neg, 2, sum)
pl_smp_tic <- sweep(pl_smp_neg, MARGIN=2, STATS=sumTIC, FUN="/")
sumTIC <- apply(pl_smp_pos, 2, sum)
pl_smp_tic <- sweep(pl_smp_pos, MARGIN=2, STATS=sumTIC, FUN="/")

sumTIC <- apply(pl_smp_neg_batch, 2, sum)
pl_smp_neg_batch_tic <- sweep(pl_smp_neg_batch, MARGIN=2, STATS=sumTIC, FUN="/")
sumTIC <- apply(pl_smp_pos_batch, 2, sum)
pl_smp_pos_batch_tic <- sweep(pl_smp_pos_batch, MARGIN=2, STATS=sumTIC, FUN="/")

## boxplots
boxplot_plot(pl_smp_neg_tic, "boxplot_peaklist_neg_sample_log_tic.pdf")
boxplot_plot(pl_smp_pos_tic, "boxplot_peaklist_pos_sample_log_tic.pdf")
boxplot_plot(pl_smp_neg_batch_tic, "boxplot_peaklist_neg_sample_log_batch_tic.pdf")
boxplot_plot(pl_smp_pos_batch_tic, "boxplot_peaklist_pos_sample_log_batch_tic.pdf")

## PCA plots
pca_plot(pl_smp_neg_tic, batch_neg, "pca_peaklist_neg_sample_log_tic.pdf")
pca_plot(pl_smp_pos_tic, batch_pos, "pca_peaklist_pos_sample_log_tic.pdf")
pca_plot(pl_smp_neg_batch_tic, batch_neg, "pca_peaklist_neg_sample_log_batch_tic.pdf")
pca_plot(pl_smp_pos_batch_tic, batch_pos, "pca_peaklist_pos_sample_log_batch_tic.pdf")

################################################################################
################## choose pl_smp_batch_tic and write to file ###################
################################################################################

mzRT <- pl_neg[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "isotopes", "adduct", "pcgroup")]
pl_final <- cbind(mzRT, pl_smp_neg_batch_tic)
write.table(pl_final, file = "peaklist_neg_log_batch_tic_final.txt", sep="\t", dec=".", quote=FALSE) 

mzRT <- pl_pos[, c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "npeaks", "isotopes", "adduct", "pcgroup")]
pl_final <- cbind(mzRT, pl_smp_pos_batch_tic)
write.table(pl_final, file = "peaklist_pos_log_batch_tic_final.txt", sep="\t", dec=".", quote=FALSE) 