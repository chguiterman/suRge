#' Identify defoliation events in host trees
#'
#' @param host_tree a data.frame rwl object containing the tree-level growth
#'   series for all host trees to be compared to the non-host chronology
#' @param nonhost_chron a data.frame rwl object comtaining a single non-host
#'   chronology
#' @param duration_years the mimimum number of years in which to consider a
#'   defolation event
#' @param max_surge the minimum level of tree growth to be considered in
#'   defoliation
#'@param bridge_events Binary, defaults to \code{FALSE}. This option allows for
#'   two successive events separated by a single year to be bridged and called one
#'   event. It should be used cautiously and closely evaluated to ensure the
#'   likelihood that the two events are actually one long event.
#' @param series_end_event Binary, defaults to \code{FALSE}. This option allows
#'   the user to identify an event ocuring at the time of sampling as a
#'   defoliation event, regardless of duration. Including it will help to
#'   quantify periodicity and extent of an outbreak. This should only be used if
#'   the user has direct knowledge of an ongoing defoliation event when the
#'   trees were sampled.
#' @param list_output defaults to \code{FALSE}. This option is to output a long
#'   list object containing a separate data.frame for each series in
#'   \code{host_tree} that includes the input series and the
#'   \code{nonhost_chron}, the corrected series, and the character string
#'   identifying the defoliation events.
#'
#' @return By default this returns a long-form data frame of tree-level growth
#'   suppression indices and identified defoliation events. If \code{list_output
#'   = TRUE}, it returns a list object with each element containing a data.frame
#'   rwl object of the host and non-host series, plus the outputs from
#'   \code{gsi}. The list object is useful for assessing the effects of running
#'   \code{gsi} on the host and nonhost data.
#'
#' @note Other functions in \code{dfoliatR}, like \code{outbreak} and
#'   \code{plot_defol}, require a long-form data frame identifiable as a
#'   \code{defol} object. Selecting \code{list_output = TRUE} will trigger
#'   errors in running other functions.
#'
#' @export
defoliate_trees <- function(host_tree, nonhost_chron, duration_years = 8,
                            max_surge = 1.28, bridge_events = FALSE,
                            series_end_event = FALSE, list_output = FALSE) {
  if(ncol(nonhost_chron) > 1) stop("nonhost_chron can only contain 1 series")
  # if(max_reduction > 0) max_reduction <- max_reduction * -1
  # To DO: Add provision that if only host series are given, no correction is made, but series are scanned for defol_status
  host_tree <- data.frame(host_tree)
  nonhost_chron <- data.frame(nonhost_chron)
  nseries <- ncol(host_tree)
  tree_list <- lapply(seq_len(nseries), function(i){
    input_series <- stats::na.omit(dplR::combine.rwl(host_tree[, i, drop=FALSE],
                                                     nonhost_chron))
    corrected_series <- dfoliatR::gsi(input_series)
    defoliated_series <- suRge::id_defoliation(corrected_series,
                                        duration_years = duration_years,
                                        bridge_events = bridge_events,
                                        max_surge = max_surge,
                                        series_end_event = series_end_event)
    return(defoliated_series)
  }
  )
  if (list_output) return(tree_list)
  else return(dfoliatR::stack_defoliation(tree_list))
}
