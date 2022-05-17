compute_elpd_ps <- function(jags_data, samples, iteration){
  
  # Unpack Simulated Data
  nf        <- jags_data$nf - 20
  nm        <- jags_data$nm - 20
  k         <- jags_data$k 
  sf        <- jags_data$sf[1:nf, 1:k]
  sm        <- jags_data$sm[1:nm, 1:k]  
  recap_f   <- jags_data$recap_f[1:nf, 1:k]
  recap_m   <- jags_data$recap_m[1:nm, 1:k]  
  recruit_f <- jags_data$recruit_f[1:nf, 1:k]
  recruit_m <- jags_data$recruit_m[1:nm, 1:k] 
  mating_f  <- jags_data$mating_f[1:nf,1:k]
  mating_m  <- jags_data$mating_m[1:nm,1:k]
  pairs_f   <- jags_data$pairs_f[1:nf, 1:k]
  
  # Unpack parameters
  PF    <- unname(samples[iteration,"PF"])
  PM    <- unname(samples[iteration,"PM"])
  rho   <- unname(samples[iteration,"rho"])
  
  # Recapture 
  sig.PF <- sqrt(PF*(1-PF))
  sig.PM <- sqrt(PM*(1-PM))
  
  Pfm <- rho * sig.PF * sig.PM + PF * PM
  P00 <- 1 - PF - PM + Pfm
  Pf0 <- PF - Pfm
  Pm0 <- PM - Pfm

  # Male recapture rates
  p_cond_m      <- PM * sm[1:nm,1:k] * recruit_m[1:nm,1:k]
  single_female <- matrix(NA, nrow = nf, ncol = k)
  p_cond_f      <- matrix(NA, nrow = nf, ncol = k)
  
  # Female recapture rates 
  for(t in 1:k){
    
    single_female[1:nf,t] <- vectorMatch(pairs_f[1:nf,t], nm + 1) 
    
    p_cond_f[1:nf, t] <- compute_prob_condF(single_female[1:nf,t],
                                            recap_m[1:nm,t],
                                            pairs_f[1:nf,t],
                                            PF,
                                            PM,
                                            Pfm,
                                            Pf0,
                                            nf,
                                            nm) * sf[1:nf,t] * recruit_f[1:nf,t]
    
  }
  
  p_cond_mf  <- rbind(p_cond_m,p_cond_f)
  recap      <- rbind(recap_m, recap_f)
  elpd_hat        <- 0

  for(i in 1:(nm + nf)){
    for(t in 1:k){
      if(is.na( p_cond_mf[i,t])|p_cond_mf[i,t]==0) next
      elpd_hat <- elpd_hat + dbinom(recap[1:(nm+nf),t], 1, p_cond_mf[1:(nm+nf),t], log =T)
    }
  }
  
  return(nll)
}


generate_pp_stat <- function(jags_data, samples){
  
  reps <- nrow(samples)
  
  nf <- jags_data$nf - 20
  nm <- jags_data$nm - 20
  n <- nf + nm 
  k <- jags_data$k
  rep_data <- vector("list", reps)
  
  for(i in 1:reps){
    print(i)
    rep_data[[i]] <- sim_dat(parameter_list = list(
      n = n,
      k = k,
      lf = 0,
      lm = 0,
      prop.female = nf/n,
      delta = rep(unname(samples[i,"delta"]), k),
      phi.f = rep(unname(samples[i,"PhiF"]), k),
      phi.m = rep(unname(samples[i,"PhiM"]), k),
      gam = rep(unname(samples[i,"gamma"]), k),
      p.f = rep(unname(samples[i,"PF"]), k),
      p.m = rep(unname(samples[i,"PM"]), k),
      rho = rep(unname(samples[i,"rho"]), k),
      betas = list(beta0 = unname(samples[i,"beta0"]), beta1 =unname(samples[i,"beta1"])),
      rand_init = F,
      init = sample(1, n, TRUE),
      show_unmated = T,
      data_aug = F
    ))
    
  }
  
}