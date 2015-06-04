#' Add the OncoSNP information to the PennCNV probe data. 
#'
#' Each probe will have it's LRR value shifted and assigned a rank 
#' (or multiple) according to the output of OncoSNP
#'
#' @param cnvDt The OncoSNP standard file (.cnvs) data.table loaded using the 
#'    loadOncosnpCnvFile() function.
#' @param qcDt The OncoSNP quality control (.qc) data.table loaded using the 
#'    LoadOncosnpQcFile() function.
#' @param probeDt The PennCNV probe data file data.table loaded using the 
#'    paramDescriptionLoadPennCNVProbeData() function.
#' @param ploidyConfig The ploidy configuration from OncoSNP that you want to 
#'    to map the information to.
#' @return A modified data.table of the probeDt that contains extra columns for
#'    shifted LRR and ranks from OncoSNP
#' @examples
#' AddOncosnp2PennCNVProbe(..., ploidyConfig = 1)
#' AddOncosnp2PennCNVProbe(..., ploidyConfig = 2)
AddOncosnp2PennCNVProbe <- function(cnvDt, qcDt, probeDt, ploidyConfig = 1){
	LRR.shift <- qcDt[ploidyNo == ploidyConfig, LRRShift]
	probeDt.modified <- dplyr::mutate(probeDt, 
													 LRRShifted = round(LRR + LRR.shift, 4),
													 rankState1 = as.numeric(NA),
													 rankState2 = as.numeric(NA),
													 rankState3 = as.numeric(NA),
													 rankState4 = as.numeric(NA),
													 rankState5 = as.numeric(NA))
	data.table::setkey(probeDt.modified, probeID)

	# Finding overlapping probes with each segment and adding shifted LRR
	# and ranks of the segment states
	cnvDf.sub <- cnvDt[ploidyNum == ploidyConfig, ]
	for (i in seq(nrow(cnvDt))) {
		# display "." for each 10 iterations for progress bar
		if (i %% 10 == 0) {
			cat('.')
		}
		segChr <- cnvDt[i, chr]
		segStart <- cnvDt[i, start]
		segEnd <- cnvDt[i, end]
		segRank <- cnvDt[i, rank]
		segState <- cnvDt[i, state]
		overlappingProbes <- probeDt.modified[chr == segChr & pos >= segStart & pos <= segEnd, 
																 probeID]
		probeDt.modified[overlappingProbes, 
						paste("rankState", segRank, sep = "") := segState, 
						with = FALSE]
	}
	return(probeDt.modified)
}
