# For running on the remote computing system

module load R

setwd("/oasis/tscc/scratch/mioverto/mutAccum/ambiRef/bam/depth_metrics")

fnFile <- "N_A_filenames.txt"
fnHead <- scan(fnFile, character(), quote = "", sep = "\n")
smplsHead <- substr(fnHead, 1, 5)
depthCols <- c("CHROM", "POS", smplsHead)

tableFile <- "N_A_coverage.tsv"
depth_df <- read.table(tableFile, header=F, col.names = depthCols)

cvrgFiles <- list.files(pattern="coverage")
nmFiles <- list.files(pattern="filenames")
for (f in 1:length(cvrgFiles)) {
    fnFile <- nmFiles[f]
    fnHead <- scan(fnFile, character(), quote = "", sep = "\n")
    smplsHead <- substr(fnHead, 1, 5)
    depthCols <- c("CHROM", "POS", smplsHead)

    tableFile <- cvrgFiles[f]
    depth_df <- read.table(tableFile, header=F, col.names = depthCols)
    for (c in 3:ncol(depth_df)) {
    # c=2
        colN <- c
        nPOS <- nrow(depth_df[depth_df[,colN] > 7,])
        sNm <- colnames(depth_df)[colN]
        sDf <- data.frame(ID=sNm, nDP8=nPOS)
        if (c == 1) {
            tblOut <- sDf
        } else {
            tblOut <- rbind(tblOut, sDf)
        } 
    }
    if (f == 1) {
        fnlTbl <- tblOut
    } else {
        fnlTbl <- rbind(fnlTbl, tblOut)
    }

}

write.table(fnlTbl, "coverage_DP8.tsv")


gLen <- 12071326