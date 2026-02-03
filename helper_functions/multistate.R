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

updateTransitionCounts <- function(msdata, edges, tm) {
  msdata_sub <- subset(msdata, status == 1)
  trans_counts <- aggregate(msdata_sub$id, by = list(transition = msdata_sub$trans),
                            FUN = function(x) length(unique(x)))
  colnames(trans_counts)[2] <- "n_patients"
  mapping <- data.frame(
    trans = seq(nrow(edges)),
    label = paste0(edges$from, " → ", edges$to),
    stringsAsFactors = FALSE
  )
  trans_counts <- merge(trans_counts, mapping, by.x = "transition", by.y = "trans", all.x = TRUE)
  trans_counts$n_patients[is.na(trans_counts$n_patients)] <- 0
  trans_counts <- trans_counts[match(mapping$label, trans_counts$label), ]
  trans_counts
}