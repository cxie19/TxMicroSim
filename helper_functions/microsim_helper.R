# function to generate the candidate times for states
gen_candidate_time <- function(paramsurv_list,num_sample,tr,model,lp_ind_all){ # num_sample is the number of patients for generation, tr is transition
  if(paramsurv_list[[tr]]$dist=="fixed"){
    return(rep(paramsurv_list[[tr]]$coefs$est,num_sample))
  }else{
    
    if(model=="piecewise"){
      gen_t <- invert_pwexp_ph(u=runif(n=num_sample,min=0,max=1), 
                               cuts=paramsurv_list[[tr]]$aux$time, 
                               h0s=sapply(1:length(paramsurv_list[[tr]]$coefs),function(x) paramsurv_list[[tr]]$coefs[[x]][1]), 
                               lp = lp_ind_all[[as.character(tr)]])
    }else if (model=="spline"){
      gamma_est <- sapply(1:length(paramsurv_list[[tr]]$coefs), # gamma estimates
                          function(x) paramsurv_list[[tr]]$coefs[[x]][1])
      gen_t <- rsurvspline(num_sample,
                           gamma=gamma_est+c(lp_ind_all[[as.character(tr)]],rep(0,length(paramsurv_list[[tr]]$coefs)-1)), 
                           knots = paramsurv_list[[tr]]$aux$knots, 
                           scale = "hazard", timescale = "log", spline = "rp")
    }
    return(gen_t)
  }
}

# The helper function is to generate the transition times and states
gen_trans_time <- function(paramsurv_list, trans_mat, sim_data, start_state, prior_intervention_state, intervention_state,
                           thor, n_sample, model){
  
  na_rows <- apply(trans_mat, 1, function(x) all(is.na(x)))
  if(sum(na_rows)>0){absorbing_state <- which(na_rows)}
  
  #risk scores in transitions not affected by treatment
  all_tr <- as.vector(trans_mat)
  all_tr <- sort(all_tr[!is.na(all_tr)])
  tx_tr <- trans_mat[,intervention_state]
  tx_tr <- tx_tr[!is.na(tx_tr)]
  non_tx_tr <- all_tr[-tx_tr]
  lp_ind_all <- setNames(vector("list", length(non_tx_tr)), as.character(non_tx_tr))
  for (tr in non_tx_tr){
    coef <- paramsurv_list[[tr]]$coefs[[1]][-1]
    lp_ind <-as.matrix(sim_data[,colnames(paramsurv_list[[tr]]$coefs[[1]])[-1]])%*%matrix(coef,byrow=F)
    lp_ind_all[[as.character(tr)]] <-lp_ind[1,1]
  }
  
  # build the microsimulation data
  microsim_data <- data.frame(id=1:n_sample,
                              from=start_state,
                              total_time=0)
  
  combined_gen_data <- vector("list", 0L)  
  i <- 0L
  
  gen_start_state <- start_state
  while(length(gen_start_state)!=0){
    
    for (gss in gen_start_state){

      next_to_state <- trans_mat[gss,] # possible transitions
      if (gss!=prior_intervention_state){next_to_state <- next_to_state[-intervention_state]}
      next_to_state <- next_to_state[!is.na(next_to_state)]
      n_gss <- sum(microsim_data$from == gss)
      
      next_to_state_time_list <- lapply(
        next_to_state,
        function(j) gen_candidate_time(paramsurv_list, num_sample = n_gss, tr = j, model = model, lp_ind_all = lp_ind_all)
      )
      
      next_to_state_time <- do.call(cbind, next_to_state_time_list)
      colnames(next_to_state_time)<- next_to_state
      min_index <- max.col(-next_to_state_time, ties.method = "first")
      microsim_data$transition[microsim_data$from==gss] <- next_to_state[min_index]
      microsim_data$min_value[microsim_data$from==gss] <- apply(next_to_state_time, 1, min)
      microsim_data$to[microsim_data$from==gss] <- sapply(microsim_data$transition[microsim_data$from==gss], function(x)
        which(trans_mat == x, arr.ind = TRUE)[2])
      microsim_data$total_time[microsim_data$from==gss] <- microsim_data$min_value[microsim_data$from==gss]+ microsim_data$total_time[microsim_data$from==gss]
    }
    microsim_data$to[microsim_data$total_time>thor] <-  microsim_data$from[microsim_data$total_time>thor]
    microsim_data$min_value[microsim_data$total_time>thor] <- microsim_data$min_value[microsim_data$total_time>thor]-(microsim_data$total_time[microsim_data$total_time>thor] - thor)
    microsim_data$total_time[microsim_data$total_time>thor] <- thor
    i <- i + 1L
    combined_gen_data[[i]] <- microsim_data
    
    microsim_data <- microsim_data[microsim_data$total_time < thor&microsim_data$to!=absorbing_state,]
    if (nrow(microsim_data)>0){
      microsim_data <- microsim_data[ , !names(microsim_data) %in% c("from","transition","min_value"), drop = FALSE]
      names(microsim_data)[names(microsim_data)=="to"] <- "from"
      microsim_data$to <- NA
      gen_start_state <- unique(microsim_data$from)
    }else{
      gen_start_state <- NULL
    }
  }
  
  result <- do.call(rbind, combined_gen_data)
  
  return(result)
  
}
