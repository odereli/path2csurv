clustered_multitask_multiple_kernel_survival_train <- function(Km, y, delta, parameters) {
  T <- length(Km)
  K <- parameters$cluster_count
  P <- dim(Km[[1]])[3]
  eta <- matrix(1 / P, K, P)
  
  Z <- matrix(runif(T*K), nrow = T)
  Z <- t(apply(Z, 1, function(x){ as.numeric(x == max(x))}))
  csums <- colSums(Z)
  update_list <- which(csums == 0)
  for (k in update_list){
    from_list <- which(csums >=2)
    if(length(from_list) == 1){
      fromCol <- from_list
    }else{
      fromCol <- sample(from_list, 1)
    }
    fromRow <- sample(which(Z[ ,fromCol] == 1), 1)
    Z[fromRow, fromCol] <- 0
    Z[fromRow, k] <- 1
    csums <- colSums(Z)
  }
  
  Z_initial <- Z
  
  Keta <- vector("list", T)
  for (t in 1:T) {
    Keta[[t]] <- calculate_Keta(Km[[t]], Z[t,] %*% eta)
  }
  
  models <- vector("list", T)
  for (t in 1:T) {
    models[[t]] <- solve_survival_svm(Keta[[t]], y[[t]], delta[[t]], parameters$C, parameters$tube, parameters$epsilon)
  }
  objectives <- sum(sapply(models, function(model) {model$objective}))
  print(c(length(objectives), sum(sapply(models, function(model) {model$objective}))))
  
  nrep <- 10
  objective <- sum(sapply(models, function(model) {model$objective}))
  temperature <- objective
  initial_temp <- temperature
  final_temp <- 1
  iter_count <- 1
  max_iter_count <- parameters$iteration_count
  task_set <- c(1:T)
  cluster_set <- c(1:K)
  
  accept_prob <- NULL
  Z_all <- list()
  
  while (1) {
    for (i in 1:nrep) {
      cvec <- matrix(0, T, K)
      for (t in 1:T) {
        for (c in 1:K) {
          divisor <- 1 / eta[c, ]
          divisor[(is.infinite(divisor))] <- 0 
          cvec[t, c] <- (t(models[[t]]$alpha) %*% calculate_Keta(Km[[t]], ((Z[t, ]*eta)*divisor)) %*% models[[t]]$alpha)
        }
      }
      Z_old <- Z
      models_old <- models
      
      t1 <- sample(task_set, 1)
      k1 <- which(Z[t1, ] == 1)
      
      new_objectives <- rep(0, K)
      obj <- Inf
      for (k2 in cluster_set) {
        Z_new <- Z_old
        Z_new[t1, ] <- 0
        Z_new[t1, k2] <- 1
        if (colSums(Z_new)[k1] == 0) {
          rownames(Z_new) <- NULL
          candidate_clusters_2 <- which(colSums(Z_new) > 1)
          if(length(candidate_clusters_2) == 1) {
            k <- candidate_clusters_2
          }else {
            k <- sample(which(colSums(Z_new) > 1), 1)
          }
          candidate_tasks <- which(Z_new[ , k] == 1)
          if(t1 %in% candidate_tasks){
            candidate_tasks <- candidate_tasks[-which(candidate_tasks == t1)]
          }
          t2 <- candidate_tasks[order(cvec[candidate_tasks, k], decreasing = T)][1]
          Z_new[t2, k] <- 0
          Z_new[t2, k1] <- 1
        }

        if(sum(colSums(Z_new) == 0) >= 1 | sum(rowSums(Z_new) != 1) >= 1) {
          print("Znew is infeasible")
          stop()
        }
        
        eta_new <- matrix(0, K, P)
        for (k in 1:K) {
          for (m in 1:P) {
            norm_tm <- sqrt(sum(sapply(1:T, function(t) {Z_new[t, k] * (t(models[[t]]$alpha) %*% Km[[t]][,,m] %*% models[[t]]$alpha)})))
            eta_new[k, m] <- eta[k, m] * norm_tm
          }
        }
        
        
        eta_new <- eta_new/ matrix(rowSums(eta_new), K, P)
        eta_new[eta_new< parameters$epsilon] <- 0
        eta_new <- eta_new/ matrix(rowSums(eta_new), K, P)
        
        for (t in 1:T) {
          Keta[[t]] <- calculate_Keta(Km[[t]], Z_new[t,] %*% eta_new)
        }
        
        for (t in 1:T) {
          models[[t]] <- solve_survival_svm(Keta[[t]], y[[t]], delta[[t]], parameters$C, parameters$tube, parameters$epsilon)
        }
        
        new_objectives[k2] <- sum(sapply(models, function(model) {model$objective}))
        
        if(new_objectives[k2] < obj){
          obj <- new_objectives[k2]
          Z_candidate <- Z_new
          eta_candidate <- eta_new
          models_candidate <- models
        }
      }
      
      Z_new <- Z_candidate
      eta_new <- eta_candidate
      models <- models_candidate
      obj_new <- obj

      if (obj_new < objective) {
        Z <- Z_new
        objective <- obj_new
        eta <- eta_new
        objectives <- c(objectives, objective)
      }else {
        n <- runif(1)
        accept_prob <- exp(-(obj_new - objective)/temperature)
        
        if(n < accept_prob) {
          Z <- Z_new
          objective <- obj_new
          eta <- eta_new
          objectives <- c(objectives, objective)
          
        }else {
          objectives <- c(objectives, objective)
          models <- models_old
        }
      }
    }
    
    temperature <- initial_temp * ((final_temp/initial_temp)^(iter_count / max_iter_count)) 

    print(iter_count)
    if (temperature <= final_temp | iter_count == max_iter_count) {
      break
    }
    iter_count <- iter_count + 1
 
  }
  
  state <- list(alpha = sapply(models, function(model) model$alpha, simplify = FALSE), b = sapply(models, function(model) model$b, simplify = FALSE), eta = eta, Z = Z, Z_initial = Z_initial, objectives = objectives, parameters = parameters)
}

