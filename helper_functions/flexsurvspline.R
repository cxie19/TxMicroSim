# For the flexible parametric proportional hazards function with the natural cubic spline function 
# for the logarithm of baseline cumulative hazard function
# The helper function is to find the best number of knots. The criterion is BIC.

fit_flexsurvspline <- function(df, covariates){
  if (length(covariates) == 0) return(NULL)
  formula_str <- paste("Surv(time, status) ~", paste(covariates, collapse = " + "))
  f <- as.formula(formula_str)
  
  BIC_check <- numeric()
  for (k in 1:10) {
    current <- try(flexsurvspline(f, data = df, k = k, scale = "hazard"), silent = TRUE)
    if (inherits(current, "try-error")) next
    BIC_check[k] <- current$BIC
  }
  if (length(BIC_check) == 0) return(NULL)
  best_k <- which.min(BIC_check)
  best_fit <- flexsurvspline(f, data = df, k = best_k, scale = "hazard")
  return(best_fit)
}


