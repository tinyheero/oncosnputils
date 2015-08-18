#' Get Oncosnp State Information
#'
#' Returns information about the Oncosnp states in a list structure. This is
#' dependent on the input state type list that was used in the Oncosnp run.
#'
#' @param state.type The "state" level.
#' @return A list with names corresponding to various details "categories" of
#'   state information
#' @export 
get_state_info <- function(state.type = "CLL") {
  if (state.type == "CLL"){
    states <- list( 
      "states" = factor(1:11), 
      "states.del" = c(1:2), 
      "states.del.loh" = c(1,2,6), # counts deletion and the somatic 2N deletion state
      "states.neu" = c(3,9:11), # counts only germline neutral states (excluding the somatic 2N LOH state)
      "states.neu.loh" = c(3,6,9:11), # counts all neutral states (including the somatic 2N LOH state)
      "states.gain" = c(4,5,7,8), # only somatic gains
      "states.gain.germline" = c(4,5,7,8,10,11), # all gain states including somatic and germline
      "states.somloh" = c(6:8), 
      "states.germloh" = c(9:11), 
      "states.som" = c(1,2,4,5,7,8),
      "states.HOMD" = c(1), 
      "states.HETD" = c(2), 
      "states.cn.neu" = 3, # copy-number neutral state 2
      "states.cn.neu.loh" = 6, # copy-number neutral LOH state 2
      "states.HLAMP" = c(5,8), # anything somatic 4N 
      "states.names" = setNames(c("HOMD", "HETD", "NEU", "3N_DUP", 
                         "4N_MONODUP", "2N_SOMLOH", "3N_SOMLOH", 
                         "4N_SOMLOH", "2N_GERMLOH", "3N_GERMLOH", 
                         "4N_GERMLOH"), 1:11),
      "states.full_names" = setNames(c("Homozygous Deletion", 
                              "Hemizygous Deletion", "2N Copy Neutral", 
                              "3N Monoallelic Duplication", 
                              "4N Monoallelic Duplication", "2N Somatic LOH", 
                              "3N Somatic LOH", "4N Somatic LOH", 
                              "2N Germline LOH", "3N Germline LOH", 
                              "4N Germline LOH"), 1:11),
      "states.cols" = setNames(c("#1F78B4", "#A6CEE3", "lightgrey", "#FB9A99", 
                        "#E31A1C", "#FDAE6B", "#F16913", "#A63603","#B2DF8A", 
                        "#33A02C", "#006837"), 1:11),
      "states.pointSymbol" = c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3),
      "copynum.cols" = setNames(c("#619CFF", "lightgrey", "F8766D"),
                        c("del", "neu", "gain"))
    )
  }
  states
}
