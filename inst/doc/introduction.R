## ------------------------------------------------------------------------
library("oncosnputils")
qcFile <- system.file("extdata", "HCC1395.qc", package = "oncosnputils")
qcDf <- LoadOncosnpQcFile(qcFile)
knitr::kable(qcDf)

## ------------------------------------------------------------------------
# Loading the OncoSNP CNV file
cnvFile <- system.file("extdata", "HCC1395.cnvs", package = "oncosnputils")
cnvDf <- LoadOncosnpCnvFile(cnvFile)

# Show only first 10 rows and 10 columns for vignette purposes
knitr::kable(cnvDf[1:10, 1:10])

## ------------------------------------------------------------------------
# Loading the PennCNV probe file
probeFile <- system.file("extdata", "logR_BAF.snp_probes.txt", package = "oncosnputils")
probeDf <- LoadPennCNVProbeData(probeFile)

# Show only first 10 rows for vignette purposes
knitr::kable(probeDf[1:10, ])

## ------------------------------------------------------------------------
head(cnvDf)

## ---- message = FALSE----------------------------------------------------
# only add for the first 10 segments for vignette purposes
cnvDf.LRR.BAF <- AddLRRBAF2OncosnpCNV(cnvDf[1:10, ], qcDf, probeDf)

## ---- message = FALSE----------------------------------------------------
# only add for the first 5000 probes for vignette purposes
probeDf.oncosnp <- AddOncosnp2PennCNVProbe(cnvDf, qcDf, probeDf[1:5000, ])
head(probeDf.oncosnp)

