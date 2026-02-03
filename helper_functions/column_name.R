# The helper function is to extract the column names of the dataset in the case-insensitive way, 
# where the column names are used for the data transformation using mstate::msprep
# given that there is any dataset uploaded

resolve_cols_ci <- function(df_names, requested) {
  out <- character(length(requested))
  for (i in seq_along(requested)) {
    want <- requested[i]
    exact_idx <- which(df_names == want)
    if (length(exact_idx) >= 1) {
      out[i] <- df_names[exact_idx[1]]
    } else {
      ci_idx <- which(tolower(df_names) == tolower(want))
      if (length(ci_idx) >= 1) {
        out[i] <- df_names[ci_idx[1]]
      } else {
        stop(sprintf("Column not found (case-insensitive): '%s'", want))
      }
    }
  }
  out
}
resolve_single_ci <- function(df_names, requested_one) resolve_cols_ci(df_names, requested_one)[1]
