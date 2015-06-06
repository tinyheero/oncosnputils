#' Add the LRR and BAF to each segment reported in a OncoSNP CNV file.
#'
#' @param cnvDt The OncoSNP standard file (.cnvs) data.table loaded using the 
#'        loadOncosnpCnvFile() function.
#' @param qcDt The OncoSNP quality control (.qc) data.table loaded using the 
#'        LoadOncosnpQcFile() function.
#' @param probeDt The PennCNV probe data file data.table loaded using the 
#'        LoadPennCNVProbeData() function.
#' @return A modified data.frame of the cnvDt that contains extra columns for
#'         LRR, LRRShifted, BAF, etc.
AddLRRBAF2OncosnpCNV <- function(cnvDt, qcDt, probeDt) {
  message("Getting LRR shift values for each ploidy")

  LRRShift.ploidyConfig1 <- dplyr::filter_(qcDt, ~ploidyNo == 1L)[, "LRRShift"]
  LRRShift.ploidyConfig1 <- as.double(LRRShift.ploidyConfig1)
  LRRShift.ploidyConfig2 <- dplyr::filter_(qcDt, ~ploidyNo == 2L)[, "LRRShift"]
  LRRShift.ploidyConfig2 <- as.double(LRRShift.ploidyConfig2)

  cnvDt.LRR.BAF <- cnvDt
  cnvDt.LRR.BAF[, "LRR"] <- NA_real_
  cnvDt.LRR.BAF[, "LRRShifted"] <- NA_real_
  cnvDt.LRR.BAF[, "BAF"] <- NA_real_
  cnvDt.LRR.BAF[, "numProbes"] <- NA_integer_
  cnvDt.LRR.BAF[, "numSnpProbes"] <- NA_integer_
  
  message("Adding LRR and BAF to each CNV segment")
  for (i in seq(nrow(cnvDt.LRR.BAF))) {
    cnv.chr <- as.character(cnvDt.LRR.BAF[i, "chr"])
    cnv.start <- as.integer(cnvDt.LRR.BAF[i, "start"])
    cnv.end <- as.integer(cnvDt.LRR.BAF[i, "end"])
    cnv.ploidyConfig <- as.integer(cnvDt.LRR.BAF[i, "ploidyNum"])

    vars <- ~chr == cnv.chr & pos >= cnv.start & pos <= cnv.end
    probeDt.sub <- dplyr::filter_(probeDt, vars)[, c("probeID", "LRR", "BAF")]

    probeDt.sub.LRRShifted <- dplyr::mutate(probeDt.sub, LRRShifted = NA_real_)
    varname <- "LRRShifted"
    if (cnv.ploidyConfig == 1L) {
      varval <- lazyeval::interp(~LRR + shift, shift = LRRShift.ploidyConfig1)
      probeDt.sub.LRRShifted <- dplyr::mutate_(probeDt.sub, 
                                               .dots = setNames(list(varval), 
                                                                varname))
    } else {
      varval <- lazyeval::interp(~LRR + shift, shift = LRRShift.ploidyConfig2)
      probeDt.sub.LRRShifted <- dplyr::mutate_(probeDt.sub, 
                                               .dots = setNames(list(varval), 
                                                                varname))
    }

    # calculate BAF of segment using only SNP probes
    probeDt.sub.LRRShifted.snpID <- dplyr::filter_(probeDt.sub.LRRShifted, 
                                                   ~grepl("^SNP", probeID))
    cnvDt.LRR.BAF[i, "BAF"] <- dplyr::summarise_(probeDt.sub.LRRShifted.snpID,
                                                 ~round(mean(BAF), 3L))

    # add additional segment information
    cnvDt.LRR.BAF[i, "LRR"] <- dplyr::summarise_(probeDt.sub.LRRShifted, 
                                          ~round(mean(LRR), 3L))
    cnvDt.LRR.BAF[i, "LRRShifted"] <- dplyr::summarise_(probeDt.sub.LRRShifted,
                                                 ~round(mean(LRRShifted), 3L))
    cnvDt.LRR.BAF[i, "numProbes"] <- nrow(probeDt.sub.LRRShifted)
    cnvDt.LRR.BAF[i, "numSnpProbes"] <- nrow(probeDt.sub.LRRShifted.snpID)
  }
  cnvDt.LRR.BAF
}
