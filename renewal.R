## The renewal equation

## function to create cases from a time-series of R values
get_I_timeseries <- function(w, R, I_init) {

  I <- rep(NA, length(R))
  I[1:length(I_init)] <- I_init

  for(t in (length(I_init) + 1):length(I)) {
    tau_end = min(length(w), t - 1)
    c_mat <- w[1:tau_end] %*% I[t - 1:tau_end]
    I[t] <- R[t] * c_mat
  }
  return(I)
}

##
get_Rt <- function(m, w) {

  R <- rep(NA, length(m))

  for(t in 2:length(m)) {
    tau_end = min(length(w), t - 1)
    c_mat <- w[1:tau_end] %*% m[t - 1:tau_end]
    R[t] <- m[t] / c_mat
  }

  return(R)
}
