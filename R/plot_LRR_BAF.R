#' Plot a LRR and BAF plot 
#'
#' @param probeDf.modified The PennCNV probe data file data.table that has been modified
#'     using the AddOncosnp2PennCNVProbe() function.
#' @return Plots the LRR and BAF 
#' @export
plot_LRR_BAF <- function(probeDf.modified, chr) {
  stateRankDf <- probeDf.modified[, paste("rankState", 1:5, sep = "")]

  probeDf.modified <- cbind(probeDf.modified, 
                            finalState = summarize_probe_state(stateRankDf, maxRank = 5))

  probeDf.modified.melt <- reshape2::melt(probeDf.modified, 
                                id.vars = c("probeID", "chr", 
                                            "pos", "finalState"),
                                measure.vars = c("LRRShifted", "BAF"), 
                                variable.name = "cnMeasure")

    ggplot2::ggplot(probeDf.modified.melt, ggplot2::aes(pos, value, color = factor(finalState))) + 
      ggplot2::geom_point(shape = 1) + 
      ggplot2::facet_grid(cnMeasure ~ ., scales = "free") + 
      ggplot2::xlab("Position") + 
      ggplot2::ylab("")
}

#' Summarize the OncoSNP state of a probe
#'
#' As OncoSNP may classify a segments belonging to several possible states,
#' each probe can belong to several states. This function resolve this and
#' identify the state that has the highest (non NA) rank from OncoSNP. 
#' Specifically, given a stateRankDf (from rankState1 to rankState5) and a 
#' maxrank, it will return the state with the highest rank. If the highest 
#' rank is NA, then it will go down the list of ranks until it finds a 
#' state that is not NA.
#' 
#' @param maxRank Specifies the highest OncoSNP rank to retrieve state values 
#'     from.
summarize_probe_state <- function(stateRankDf, maxRank = 5) {
	states.summarized <- stateRankDf[[maxRank]]
	if (maxRank != 1) {
		for (i in (maxRank-1):1) {
			states.summarized[is.na(states.summarized)] <- 
        stateRankDf[is.na(states.summarized), i][[1]]
		}
	}
	states.summarized
}
