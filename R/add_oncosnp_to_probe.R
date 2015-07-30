#' Add the OncoSNP information to the PennCNV probe data. 
#'
#' Each probe will have its LRR value shifted and assigned a rank 
#' (or multiple) according to the output of OncoSNP
#'
#' @param cnvDf The OncoSNP standard file (.cnvs) data.frame loaded using the 
#'    loadOncosnpCnvFile() function.
#' @param qcDf The OncoSNP quality control (.qc) data.frame loaded using the 
#'    LoadOncosnpQcFile() function.
#' @param probeDf The PennCNV probe data file data.frame loaded using the 
#'    LoadPennCNVProbeData() function.
#' @param ploidyConfig The ploidy configuration from OncoSNP that you want to 
#'    to map the information to.
#' @return A modified data.frame of the probeDf that contains extra columns for
#'    shifted LRR and ranks from OncoSNP
#' @examples
#' AddOncosnp2PennCNVProbe(..., ploidyConfig = 1)
#' AddOncosnp2PennCNVProbe(..., ploidyConfig = 2)
AddOncosnp2PennCNVProbe <- function(cnvDf, qcDf, probeDf, ploidyConfig = 1L){
  LRR.shift <- dplyr::filter_(qcDf, ~ploidyNo == ploidyConfig)[, "LRRShift"]
  LRR.shift <- as.double(LRR.shift)

  probeDf[, "LRRShifted"] <- round(probeDf$LRR + LRR.shift, 4L)
  probeDf[, paste("rankState", 1:5, sep = "")] <- NA_integer_

  # Finding overlapping probes with each segment and adding shifted LRR
  # and ranks of the segment states
  message("Adding shifted LRR and ranks for each probe")
  cnvDf.sub <- dplyr::filter_(cnvDf, ~ploidyNum == ploidyConfig)
  for (i in seq(nrow(cnvDf.sub))) {
    segChr <- as.character(cnvDf.sub[i, "chr"])
    segStart <- as.integer(cnvDf.sub[i, "start"])
    segEnd <- as.integer(cnvDf.sub[i, "end"])
    segRank <- as.integer(cnvDf.sub[i, "rank"])
    segState <- as.integer(cnvDf.sub[i, "state"])
    
    vars <- ~chr == segChr & pos >= segStart & pos <= segEnd
    overlappingProbes <- unlist(dplyr::filter_(probeDf, vars)[, "probeID"])
    probeDf.overlappingProbes.flag <- probeDf$probeID %in% overlappingProbes
    probeDf[probeDf.overlappingProbes.flag, 
            paste("rankState", segRank, sep = "")] <- segState
  }
  probeDf
}
