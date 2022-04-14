clustered_multitask_multiple_kernel_survival_test <- function(Km, state) {
  T <- length(Km)
  Z <- state$Z
    Keta <- vector("list", T)
  for (t in 1:T) {
    Keta[[t]] <- calculate_Keta(Km[[t]], Z[t,] %*% state$eta) 
  }
  
  y <- vector("list", T)
  for (t in 1:T) {
    y[[t]] <- Keta[[t]] %*% state$alpha[[t]] + state$b[[t]] 
  }
  
  prediction <- list(y = y)
}

