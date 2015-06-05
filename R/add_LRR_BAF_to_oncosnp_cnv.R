#' Add the LRR and BAF to each segment reported in a OncoSNP CNV file.
#'
#' @param cnvDt The OncoSNP standard file (.cnvs) data.table loaded using the 
#'				loadOncosnpCnvFile() function.
#' @param qcDt The OncoSNP quality control (.qc) data.table loaded using the 
#'				LoadOncosnpQcFile() function.
#' @param probeDt The PennCNV probe data file data.table loaded using the 
#'        LoadPennCNVProbeData() function.
#' @return A modified data.table of the cnvDt that contains extra columns for
#'				 LRR, LRRShifted, BAF, etc.
AddLRRBAF2OncosnpCNV <- function(cnvDt, qcDt, probeDt) {
	message("Getting LRR shift values for each ploidy")
	LRRShift.ploidyConfig1 <- qcDt[ploidyNo == 1L, LRRShift]
	LRRShift.ploidyConfig2 <- qcDt[ploidyNo == 2L, LRRShift]
	cnvDt.LRR.BAF <- dplyr::mutate(cnvDt, 
                                 LRR = as.double(NA), 
                                 LRRShifted = as.double(NA),
                                 BAF = as.double(NA), 
                                 numProbes = as.numeric(NA), 
                                 numSnpProbes = as.numeric(NA))
	
	message("Adding LRR and BAF to each CNV segment")
	for (i in seq(nrow(cnvDt.LRR.BAF))) {
		# display "." for each 10 iterations for progress bar
		if (i %% 10 == 0) {
			cat(".")
		}
		cnv.chr <- cnvDt.LRR.BAF[i, chr]
		cnv.start <- cnvDt.LRR.BAF[i, start]
		cnv.end <- cnvDt.LRR.BAF[i, end]
		cnv.ploidyConfig <- cnvDt.LRR.BAF[i, ploidyNum]

		probeDt.sub <- probeDt[chr == cnv.chr & pos >= cnv.start & pos <= cnv.end, 
													 c("probeID", "LRR", "BAF"), with = FALSE]

		probeDt.sub.LRRShifted <- dplyr::mutate(probeDt.sub, LRRShifted = NA)
		if (cnv.ploidyConfig == 1) {
			probeDt.sub.LRRShifted[, LRRShifted := LRR + LRRShift.ploidyConfig1]
		} else {
			probeDt.sub.LRRShifted[, LRRShifted := LRR + LRRShift.ploidyConfig2]
		}

		# calculate BAF of segment using only SNP probes
		probeDt.sub.LRRShifted.snpID <- probeDt.sub.LRRShifted[grep("^SNP", probeID), ]
		cnvDt.LRR.BAF[i, BAF := round(mean(probeDt.sub.LRRShifted.snpID[, BAF]), 3)]

		# add additional segment information
		cnvDt.LRR.BAF[i, LRR := round(mean(probeDt.sub.LRRShifted[, LRR]), 3)]
		cnvDt.LRR.BAF[i, LRRShifted := round(mean(probeDt.sub.LRRShifted[, LRRShifted]), 3)]
		cnvDt.LRR.BAF[i, numProbes := nrow(probeDt.sub)]
		cnvDt.LRR.BAF[i, numSnpProbes := nrow(probeDt.sub.LRRShifted.snpID)]
	}
	return(cnvDt.LRR.BAF)
}
