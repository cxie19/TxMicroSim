# The helper function is to plot the multistate model's diagram

circle_layout_custom <- function(labels, radius = 1, start_angle = 0, clockwise = TRUE) {
  n <- length(labels)
  if (n == 0) return(data.frame())
  
  if (clockwise) {
    angles <- seq(start_angle, start_angle - 2*pi, length.out = n + 1)[- (n+1)]
  } else {
    angles <- seq(start_angle, start_angle + 2*pi, length.out = n + 1)[- (n+1)]
  }
  
  data.frame(
    id    = labels,
    label = labels,
    x     = radius * cos(angles),
    y     = -radius * sin(angles)   # flip Y for visNetwork
  )
}

# The helper function is to count the number of patients of each transition from the multistate model 

updateTransitionCounts <- function(msdata, edges = NULL, tm = NULL) {
  
  mapping <- data.frame(
    trans = seq(nrow(edges)),
    label = paste0(edges$from, " → ",edges$to))
  
  ms_sub <- msdata[msdata$status == 1 & !is.na(msdata$trans), , drop = FALSE]
  
  if (nrow(ms_sub) == 0) {
    # No events at all -> return all transitions with 0
    mapping$n_patients <- 0L
    return(mapping)
  }
  
  trans_counts <- aggregate(
    id ~ trans,
    data = ms_sub,
    FUN = function(x) length(unique(x))
  )
  names(trans_counts)[2] <- "n_patients"
  
  out <- merge(mapping, trans_counts, by = "trans", all.x = TRUE, sort = FALSE)
  out$n_patients[is.na(out$n_patients)] <- 0L
  
  out <- out[match(mapping$trans, out$trans), , drop = FALSE]
  
  return(out)
}
