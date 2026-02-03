# The helper function is to generate the transition times and states
gen_trans_time <- function(paramsurv_list, trans_mat, sim_data, start_state, prior_intervention_state, intervention_state,
                           thor, n_sample, model){
  
  na_rows <- apply(trans_mat, 1, function(x) all(is.na(x)))
  if(sum(na_rows)>0){absorbing_state <- which(na_rows)}
  
  next_to_state <- trans_mat[start_state,]
  if (start_state!=prior_intervention_state){next_to_state <- next_to_state[-intervention_state]}
  next_to_state <- next_to_state[!is.na(next_to_state)]
  if(!is.null(next_to_state)){
    next_to_state_time <- foreach(j= next_to_state,.combine = "cbind")%do%{
      
      if(paramsurv_list[[j]]$dist=="fixed"){
        return(rep(paramsurv_list[[j]]$coefs$est,n_sample))
      }
      else{
        coef=paramsurv_list[[j]]$coefs[[1]][-1]
        lp_ind <-as.matrix(sim_data[,colnames(paramsurv_list[[j]]$coefs[[1]])[-1]])%*%matrix(coef,byrow=F)
        lp_ind <-lp_ind[1,1]
        
        if(model=="piecewise"){
          gen_t <- invert_pwexp_ph(u=runif(n=n_sample,min=0,max=1), 
                                   cuts=paramsurv_list[[j]]$aux$time, 
                                   h0s=sapply(1:length(paramsurv_list[[j]]$coefs),function(x) paramsurv_list[[j]]$coefs[[x]][1]), 
                                   lp = lp_ind)
        }else if (model=="spline"){
          # gamma estimates
          gamma_est <- sapply(1:length(paramsurv_list[[j]]$coefs),
                              function(x) paramsurv_list[[j]]$coefs[[x]][1])
          
          gen_t <- rsurvspline(n_sample,
                               gamma=gamma_est+c(lp_ind,rep(0,length(paramsurv_list[[j]]$coefs)-1)), 
                               knots = paramsurv_list[[j]]$aux$knots, 
                               scale = "hazard", timescale = "log", spline = "rp")
        }
        
        return(c(gen_t))
      }
    }
    next_to_state_time <- as.matrix(next_to_state_time,byrow=F,nrow=length(next_to_state))
    colnames(next_to_state_time)<- next_to_state
    next_to_state_time <- as.data.frame(next_to_state_time)
    min_index <- max.col(-next_to_state_time[, 1:length(next_to_state)], ties.method = "first")
    microsim_data <- data.frame(id=1:n_sample,
                                from=start_state,
                                transition = next_to_state[min_index],
                                min_value = apply(next_to_state_time, 1, min))
    # microsim_data <- cbind(microsim_data,next_to_state_time)
    microsim_data$to <- sapply(microsim_data$transition, function(x)
      which(trans_mat == x, arr.ind = TRUE)[2])
    
  }
  
  # patients still at risk 
  go_to_next_trans <- microsim_data[microsim_data$min_value < thor,]
  next_to_state <- unique(go_to_next_trans$to)
  if(!is.null(absorbing_state)){  #absorbing state
    next_to_state <- next_to_state[!next_to_state%in%absorbing_state]
    go_to_next_trans <- go_to_next_trans[go_to_next_trans$to!=absorbing_state,]
  }
  
  while(length(next_to_state)!=0){
    
    microsim_data$total_time <- NULL
    
    next_to_state_time <- foreach(j=next_to_state)%do%{
      
      potential_next_tr <-  trans_mat[j,]
      potential_next_tr <- potential_next_tr[!is.na(potential_next_tr)]
      num_next_trans_j <- sum(go_to_next_trans$to==j)
      
      potential_next_tr_time <- foreach(k=potential_next_tr,.combine = "cbind")%do%{
        if(paramsurv_list[[k]]$dist=="fixed"){
          return(rep(paramsurv_list[[k]]$coefs$est,num_next_trans_j))
        }else{
          coef=paramsurv_list[[k]]$coefs[[1]][-1]
          lp_ind <-as.matrix(sim_data[,colnames(paramsurv_list[[k]]$coefs[[1]])[-1]])%*%matrix(coef,byrow=F)
          lp_ind <-lp_ind[1,1]
          
          if(model=="piecewise"){
            gen_t <- invert_pwexp_ph(u=runif(n=num_next_trans_j,min=0,max=1), 
                                     cuts=paramsurv_list[[k]]$aux$time, 
                                     h0s=sapply(1:length(paramsurv_list[[k]]$coefs),function(x) paramsurv_list[[k]]$coefs[[x]][1]), 
                                     lp = lp_ind)
          }else if (model=="spline"){
            # gamma estimates
            gamma_est <- sapply(1:length(paramsurv_list[[k]]$coefs),
                                function(x) paramsurv_list[[k]]$coefs[[x]][1])
            
            gen_t <- rsurvspline(num_next_trans_j,
                                 gamma=gamma_est+c(lp_ind,rep(0,length(paramsurv_list[[k]]$coefs)-1)), 
                                 knots = paramsurv_list[[k]]$aux$knots, 
                                 scale = "hazard", timescale = "log", spline = "rp")
          }
          
          return(c(gen_t))
        }
      }
      
      potential_next_tr_time <- as.matrix(potential_next_tr_time,byrow=F,nrow=length(potential_next_tr))
      colnames(potential_next_tr_time)<- potential_next_tr
      potential_next_tr_time <- as.data.frame(potential_next_tr_time)
      min_index <- max.col(-potential_next_tr_time[, 1:length(potential_next_tr)], ties.method = "first")
      add_microsim_data <- data.frame(id=go_to_next_trans$id[go_to_next_trans$to==j],
                                      from=j,
                                      transition = potential_next_tr[min_index],
                                      min_value = apply(potential_next_tr_time, 1, min))
      # microsim_data <- cbind(microsim_data,next_to_state_time)
      add_microsim_data$to <- sapply(add_microsim_data$transition, function(x)
        which(trans_mat == x, arr.ind = TRUE)[2])
      
      microsim_data <- rbind(microsim_data,add_microsim_data)
      return(add_microsim_data)
    }
    
    microsim_data <- microsim_data %>% 
      arrange(id) %>% 
      group_by(id) %>% 
      mutate(total_time=cumsum(min_value))
    
    # patients still at risk 
    go_to_next_trans <- microsim_data %>%
      group_by(id) %>%
      filter(row_number() == n() & total_time < thor)
    next_to_state <- unique(go_to_next_trans$to)
    if(!is.null(absorbing_state)){  #absorbing state
      next_to_state <- next_to_state[!next_to_state%in%absorbing_state]
      go_to_next_trans <- go_to_next_trans[go_to_next_trans$to!=absorbing_state,]
    }
  }
  
  if(is.null(microsim_data$total_time)){
    microsim_data$total_time <- microsim_data$min_value
  }
  microsim_data$to[microsim_data$total_time>thor] <- microsim_data$from[microsim_data$total_time>thor]
  microsim_data$total_time[microsim_data$total_time>thor] <- thor
  return(microsim_data)
  
}