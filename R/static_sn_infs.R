#' Simulation of outbreak in static network
#' @export
sim_outbreak_static_sn <- function(N=1e4, 
                                   gamma=1.8, k0=0.5, phi=1e4, 
                                   n_infs = 1, 
                                   beta = 1, psi = 1,  sigma = 1, 
                                   nu = 1, eta = 1, mu = 1, rho = 0,
                                   t=1, max.iter=1e6,sn = 0) {
  
  start <- Sys.time()
  
  inputs <- list(N = N, sn = sn, 
                 gamma = gamma, k0 = k0, phi = phi,
                 n_infs = n_infs,
                 beta = beta, psi = psi, sigma = sigma,
                 nu = nu, eta = eta, mu =  mu, rho = rho, 
                 t = t, max.iter = max.iter)
  
  time.window <- c(0, t)
  time <- time.window[1]
  tau <- 0 # initialise
  times <- rep(0, max.iter) # storing times at which events happen
  
  if(!is.list(sn)) {
    sn <- sim_static_sn(gamma = gamma, k0 = k0, N = N, phi = phi)
    sn$rec_rels <- which(sn$res == 1, arr.ind = T)
    colnames(sn$rec_rels) <- c("p1", "p2")
  }
  
  lambdas <- sn$lambdas
  
  rels <- sn$rec_rels[, c("p1", "p2")]
  dd <- sn$dd
  prop0 <- sn$prop0
  const <- sn$const
  
  NRel <- nrow(rels)
  N_vec <- 1:N # used to speed up sampling in loop
  
  n_stage <- 4 # (U, E, A, T)
  
  infs <- matrix(F, nrow = N,  ncol = n_stage, dimnames = list(NULL, c("U", "E", "A", "T")))
  repeat_infs <- rep(0, N) # vector used to calculate total number of infections for each infectee
  N_vec <- 1:N # used to speed up sampling in loop
  w <- order(lambdas, decreasing = T)[1:n_infs] # choose initial infections
  infs[w, "A"] <- T
  
  strains <- rep(0, N)
  strains[w] <- 1 # add strain to each infectee
  
  Ninf <- rep(0, n_stage)
  names(Ninf) <- c("U", "E", "A", "T")
  Ninf["A"] <- n_infs # introduce infections
  Ninf0 <- Ninf
  
  symptoms_vec <- sample(c("E", "A"), size = max.iter, replace = T, prob = c(psi, 1 - psi)) # pre-generate probs of developing symptoms
  
  # create objects to record number of replacement infections vs new infections
  new_inf_event <- 0  # record which strains are causing new infections
  log_infs <- matrix(data = NA, nrow = max.iter, ncol = 8, dimnames = list(NULL, c("time", "c1", "c2",  "p", "infector", "infector_ninf","infector_stage", "p_ninf"))) 
  
   
  rinf <- NRel * beta  # R^_i(t_b) # set initial rate of infs = 0 (as implemented after burn-in)
  rsymp <- 0 # total rate of developing symptoms
  rtreat <- 0 # set initial rate of seeking treatment = 0 (as implemented after burn-in)
  rscreen <-  eta * n_infs # increase screening rate by eta
  rrec <- nu * n_infs # increase recovery rate by nu 
  rcure <- 0 # set initial rate of cure (or failure) = 0 (as implemented after burn-in)
  
  step <- 1
  
  while(time < time.window[2] & step <= max.iter){
    if (sum(Ninf)==0) break

    prob <- c(rinf, rsymp, rtreat, rrec, rscreen, rcure) # now a vector of length 6
   
    e <- sample(c('rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure'), size = 1, prob = prob)
    

    tau <- rexp(n = 1, rate = sum(prob)) # generate event time
    time <- time + tau # update time
    times[step] <- time    #store time
    
    if (e=='rinf') { #infection event
      
      w <- sample.int(NRel, 1) # choose one couple
      w_inf_status <- strains[rels[w, c("p1", "p2")]] # find their infection statuses
      
      if( !any(w_inf_status == 1)) next # neither are infected 
      if( any(w_inf_status == 99)) next # one is vaccinated
      if( all(w_inf_status > 0 )) next # both are infected
      if( any(infs[rels[w, c("p1", "p2")], "T"])) next # in treatment so not infectious
      
      infectee <- rels[w, which(w_inf_status != 1)]
      infector <- rels[w, which(w_inf_status == 1)]
      
      rsymp <- rsymp + sigma # increase overall rate of leaving incubation stage
      new_inf_event <- new_inf_event + 1 # increase ticker
      infs[infectee, "U"] <- T # record infected status for infectee
      strains[infectee] <- 1
      Ninf[ "U"] <- Ninf[ "U"] + 1 # increase total infs 1
      repeat_infs[infectee] <- repeat_infs[infectee] + 1 # increase total number of infections for infectee
      log_infs[step, c("time", "c1", "c2", "p", "infector", "infector_ninf","infector_stage", "p_ninf")] <- c(time, 1, 2, infectee, infector, repeat_infs[infector], which(infs[infector, ])+1, repeat_infs[infectee])
    }      
    
    else if (e == 'rsymp') { #infection leaves incubation stage
      
      if(Ninf["U"] == 1) w <- N_vec[infs[,"U"]] # choose one incubating individual
      else w <- sample(x = N_vec[infs[,"U"]], size = 1) # equal probability of leaving incubation stage for all new infections
      
      infs[w, c("U", symptoms_vec[step])] <- c(F, T) # add infection to either symptomatic or asymptomatic compartment
      Ninf[ c("U", symptoms_vec[step])] <- Ninf[ c("U", symptoms_vec[step])] + c(-1, 1) # remove infection from U 
      
      rsymp <- rsymp - sigma # decrease overall rate of leaving incubation stage
      
      if (symptoms_vec[step] == "E") { # if symptoms develop
        rtreat <- rtreat + mu # increase treatment rate by mu
        # record event
        log_infs[step, c("time", "c1", "c2",  "p", "p_ninf")] <- c(time, 2, 3,  w, repeat_infs[w])
      }
      else if (symptoms_vec[step] == "A") { # if infection is asymptomatic
        rscreen <- rscreen + eta # increase screening rate by eta
        rrec <- rrec + nu # increase recovery rate by nu 
        # record event
        log_infs[step, c("time", "c1", "c2",  "p", "p_ninf")] <- c(time, 2, 4,  w, repeat_infs[w])
      }
    }
    
    else if (e == 'rtreat') { # symptomatic individual seeks treatment
      
      if(Ninf["E"] == 1) w <- N_vec[infs[,"E"]]  # choose one symptomatic individual
      else w <- sample(x = N_vec[infs[,"E"]], size = 1) # equal probability of treatment seeking for all symptomatic carriers
      
      infs[w, c("E", "T")] <- c(F, T) # remove infection from symptomatic stage, add infection to treatment stage
      Ninf[c("E", "T")] <- Ninf[c("E", "T")] + c(-1, 1) # remove infection from symptomatic stage, add infection to treatment stage
      
      rtreat <- rtreat - mu # decrease overall rate of seeking treatment by mu
      rcure <- rcure + rho # increase the overall rate of cure by rho
      # record event
      log_infs[step, c("time", "c1", "c2",  "p", "p_ninf") ] <- c(time, 3, 5,  w, repeat_infs[w])
    }
    
    else if (e == 'rscreen') { # asymptomatic individual is screened
      
      if(Ninf["A"] == 1) w <- N_vec[infs[, "A"]] # choose one asymptomatic individual
      else w <- sample(x = N_vec[infs[, "A"]], size = 1) # equal probability of being screened for all asymptomatic carriers
      
      infs[w, c("A", "T")] <- c(F, T) # remove infection from asymptomatic stage, add infection to treatment stage
      Ninf[c("A", "T")] <- Ninf[ c("A", "T")] + c(-1, 1) # remove infection from asymptomatic stage, add infection to treatment stage
      
      rrec <- rrec - nu # decrease overall rate of recovery by nu 
      rscreen <- rscreen - eta # decrease screening rate by eta
      rcure <- rcure + rho # increase the overall rate of cure by rho
      # record event
      log_infs[step, c("time", "c1", "c2", "p", "p_ninf")] <- c(time, 4, 5,  w, repeat_infs[w])
    }
    
    else if (e == 'rrec') { # asymptomatic individual recovers
      
      if(Ninf["A"] == 1) w <- N_vec[strains == 1 & infs[, "A"]] # choose one asymptomatic individual
      else w <- sample(N_vec[strains == 1 & infs[, "A"]], 1) # equal probability of recovery for all asymptomatic carriers with that strain
      
      Ninf[ "A"] <- Ninf[ "A"] - 1 # decrease total infs by 1 
      infs[w, "A"] <- F # remove infection from asymptomatic stage and set strain to 0
      strains[w] <- 0
      
      rrec <- rrec - nu # decrease overall rate of recovery by nu
      rscreen <- rscreen - eta # also decrease screening rate by eta
      # record event
      log_infs[step, c("time", "c1", "c2", "p", "p_ninf")] <- c(time, 4, 1, w, repeat_infs[w])
    }
    
    else if (e == 'rcure') { # treated individual is cured
      
      if(Ninf["T"] == 1) w <- N_vec[infs[, "T"]] # choose one individual in treatment 
      else w <- sample(x = N_vec[infs[, "T"]], size = 1) # equal probability of cure for all treated infections
      
      Ninf["T"] <- Ninf[ "T"] - 1 # decrease total infs by 1
      infs[w, "T"] <- F # remove infection from treatment stage and set strain to 0
      rcure <- rcure - rho # decrease the overall rate of cure by rho
      strains[w] <- 0
      # record event
      log_infs[step, c("time", "c1", "c2", "p", "p_ninf")] <- c(time, 5, 1, w, repeat_infs[w])
    }
    
    step <- step + 1
  }
  
  #print(Sys.time() - start)
  
  # remove unused storage
  log_infs <- log_infs[!is.na(log_infs[, 1]), , drop=F]
  
  
  log_infs  <- as.data.frame(log_infs)
  if( nrow(log_infs) > 0 ) {
    log_infs$c1 <- factor(x = log_infs$c1, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$c2 <- factor(x = log_infs$c2, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$infector_stage <- factor(x = log_infs$infector_stage, levels = 2:4, labels = c( "U", "E", "A"), exclude = NA)
    log_infs$infector_full <- paste0(log_infs$infector, "_", log_infs$infector_ninf)
    log_infs$p_full <- paste0(log_infs$p, "_", log_infs$p_ninf)
    }
  
  sn <- list(inputs = inputs, kmax = sn$kmax, const = sn$const, 
             log_infs = log_infs,
             dd = sn$dd, prop0 = sn$prop0,
             comp_time = Sys.time () - start)
  
  return(sn)
}
