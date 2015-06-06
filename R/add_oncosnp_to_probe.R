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
AddOncosnp2PennCNVProbe <- function(cnvDf, qcDf, probeDf, ploidyConfig = 1){
  LRR.shift <- dplyr::filter_(qcDf, ~ploidyNo == ploidyConfig)[, "LRRShift"]
  LRR.shift <- as.double(LRR.shift)

  varval <- lazyeval::interp(~round(LRR + LRR.shift, 4))
  probeDf.modified <- dplyr::mutate_(probeDf, 
                                     .dots = setNames(list(varval), 
                                                      "LRRShifted"))
  probeDf.modified[, paste("rankState", 1:5, sep = "")] <- as.numeric(NA)

  # Finding overlapping probes with each segment and adding shifted LRR
  # and ranks of the segment states
  message("Adding shifted LRR and ranks for each probe")
  cnvDf.sub <- dplyr::filter_(cnvDf, ~ploidyNum == ploidyConfig)
  for (i in seq(nrow(cnvDf))) {
    # display "." for each 10 iterations for progress bar
    if (i %% 10 == 0) {
      cat(".")
    }

    segChr <- as.character(cnvDf[i, "chr"])
    segStart <- as.numeric(cnvDf[i, "start"])
    segEnd <- as.numeric(cnvDf[i, "end"])
    segRank <- as.numeric(cnvDf[i, "rank"])
    segState <- as.numeric(cnvDf[i, "state"])
    
    vars <- ~chr == segChr & pos >= segStart & pos <= segEnd
    overlappingProbes <- unlist(dplyr::filter_(probeDf.modified, vars)[, "probeID"])
    varval <- lazyeval::interp(~ifelse(probeID %in% overlappingProbes, 
                                       segState, 
                                       NA))
    probeDf.modified <- dplyr::mutate_(probeDf.modified,
                                       .dots = setNames(list(varval), 
                                                        paste("rankState", 
                                                              segRank, 
                                                              sep = "")))
  }
  return(probeDf.modified)
}
