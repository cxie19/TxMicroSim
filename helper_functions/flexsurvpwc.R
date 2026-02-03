# The helper functions are related to the flexible parametric proportional-hazards function
# with the piecewise constant (exponential) function as the baseline hazard function

# split the number of events into time intervals 
make_event_intervals <- function(data, events_per_interval, max_intervals, max_time){
  # times: vector of event/censoring times
  # status: 1 = event, 0 = censored
  # events_per_interval: target number of events in each interval
  times = data$time
  status = data$status
  
  # Extract event times
  event_times <- times[status == 1]
  n_events <- length(event_times)
  
  if (n_events == 0) return(numeric(0))  # no events
  
  # How many intervals can we make?
  n_intervals <- max(1, ceiling(n_events / events_per_interval))
  n_intervals <- min(n_intervals,max_intervals)
  
  if (n_intervals == 1) {
    return(numeric(0))  # no cutpoints needed, constant hazard
  }
  
  # Quantiles of event times
  probs <- seq(0, 1, length.out = n_intervals+1 )
  cuts <- quantile(event_times, probs = probs)
  cuts[1] <- 0
  cuts[length(cuts)] <- max_time
  cuts <- unique(sort(cuts))
  
  return(as.numeric(cuts))
}

# plot the baseline hazard function and the baseline cumulative hazard function
plot_piecewise_hazard <- function(cuts, haz, tr, xmax_val,ylim_haz = NULL, ylim_cumhaz = NULL) {
  
  # Initialize with full cuts and cumulative hazard
  H_vals <- numeric(length(cuts))
  for (i in seq_len(length(cuts) - 1)) {
    left <- cuts[i]
    right <- cuts[i + 1]
    H_vals[i + 1] <- H_vals[i] + (right - left) * haz[i]
  }
  
  # Default cumulative hazard up to all cuts
  cuts_ext <- cuts
  H_ext <- H_vals
  
  # Default limits
  if (is.null(ylim_haz)) ylim_haz <- max(haz)
  if (is.null(ylim_cumhaz)) ylim_cumhaz <- max(H_ext)
  
  par(mfrow = c(1, 2))
  # Plot hazard
  plot(NULL, xlim = c(0, xmax_val), ylim = c(0, ylim_haz),
       xlab = "Time", ylab = "Baseline hazard function",
       main = paste("Transition", tr, ": Baseline hazard"))
  for (i in seq_along(haz)) {
    lines(c(cuts[i], cuts[i+1]), rep(haz[i], 2), lwd = 2, col = "darkblue")
    if (i < length(haz)) {lines(rep(cuts[i+1], 2), c(haz[i], haz[i+1]), lty = 2, col = "gray")}
  }
  
  plot(cuts_ext, H_ext, type = "l", lwd = 2, col = "darkred",
       xlab = "Time", ylab = "Baseline cumulative hazard function",
       main = paste("Transition", tr, ": Baseline cumulative hazard"),
       xlim = c(0,xmax_val),ylim = c(0, ylim_cumhaz))
  
  par(mfrow = c(1, 1))
  
}

# Obtain the survival time 
invert_pwexp_ph <- function(u, cuts, h0s, lp = NULL) {
  stopifnot(is.numeric(u), all(u > 0 & u < 1),
            is.numeric(cuts), is.numeric(h0s),
            length(h0s) == length(cuts) )
  
  # PH multiplier
  w <- exp(lp)
  # recycle w if scalar
  if (length(w) == 1L) w <- rep(w, length(u))
  
  # Target *baseline* cumulative hazard to hit
  Hstar <- -log(u) / w
  
  # Precompute baseline cumulative hazard at each cut (finite part only)
  widths <- diff(c(0,cuts))              # K-1 widths (finite intervals)
  H0_at_cuts <- c(0, cumsum(h0s * widths))  # length K
  
  # For each Hstar, find the interval index k such that:
  # H0_at_cuts[k] <= Hstar < H0_at_cuts[k+1] for k = 1..K-1,
  # or k = K when it falls beyond last cut (into the last open interval)
  k <- findInterval(Hstar, H0_at_cuts, rightmost.closed = TRUE) 
  # clamp to [1, K]
  K <- length(h0s)
  k[k < 1L] <- 1L
  k[k > K]  <- K
  
  # Start time of interval k and its upper time bound
  t_start <- c(0,cuts)[k]                 # length(u)
  # cumulative baseline hazard accumulated up to t_start
  H0_start <- H0_at_cuts[k]
  
  # time increment needed inside interval k to accumulate remaining hazard
  # deltaH = Hstar - H0_start; since h0 is constant = h0s[k], delta_t = deltaH / h0s[k]
  dt <- (Hstar - H0_start) / h0s[k]
  # dt <- (Hstar[8] - H0_start[8]) / h0s[k][8]
  
  # Final time; for k == K this naturally goes beyond last cut
  t <- t_start + dt
  
  # Guard tiny negative numerical noise (e.g., u ~ 1 gives ~0 time)
  return(pmax(t, 0))
  
}
