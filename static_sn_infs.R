source("sn/src/static_sn.R")
source("sn/src/extract_offspring_dist.R")
source("general/mean_ci.R")
sim_outbreak_static_sn <- function(N, sn = 0, 
                                   gamma = 1.5, k0 = 0.3, phi = 1e4, kmax = N-1, 
                                   lambdas = 0, opt_phi = TRUE,
                                   n_strain = 1,  n_infs = 0, prop_pres = 1, efficacy = 1,
                                   beta = 1, psi = 1,  sigma = 1, alpha = 1,
                                   nu = 1, eta = 1, mu = 1, rho = 0,
                                   t, max.iter, burn.in = 1/365) {
  
  start <- Sys.time()
  
  inputs <- list(N = N, sn = sn, 
                 gamma = gamma, k0 = k0, kmax = kmax, phi = phi, lambdas = lambdas,
                 n_strain = n_strain,  n_infs = n_infs, prop_pres = prop_pres, efficacy = efficacy,
                 beta = beta, psi = psi, sigma = sigma, alpha = alpha,
                 nu = nu, eta = eta, mu =  mu, rho = rho, 
                 t = t, max.iter = max.iter, burn.in = burn.in)
  
  if (n_strain != length(beta)) stop("number of strains and supplied betas do not match")
  if (n_strain != length(alpha)) stop("number of strains and supplied alphas do not match")
  if (n_strain != length(n_infs)) stop("number of strains and supplied infs do not match")
  if(sum(prop_pres) != 1) stop('proportion of prescriptions does not sum to 1')
  if(sum(dim(efficacy) - n_strain) != 0) stop('efficacy matrix does not match number of strains')
  if(any(efficacy > 1) | any(efficacy < 0)) stop('efficacy must be between 0 and 1')
  
  
  time.window <- c(burn.in, t)
  time <- time.window[1]
  tau <- 0 # initialise
  times <- rep(0, max.iter) # storing times at which events happen
  
  efficacy <- as.matrix(efficacy)
  nualpha <- nu * alpha
  n_pres <- length(prop_pres)
  
  if(!is.list(sn)) {
    sn <- sim_static_sn(gamma = gamma, k0 = k0, N = N, plot = F, kmax = kmax, mu = 1, phi = phi)
    sn$rec_rels <- which(sn$res == 1, arr.ind = T)
    colnames(sn$rec_rels) <- c("p1", "p2")
    lambdas <- sn$inputs$lambdas
  } else {
    lambdas <- sn$inputs$lambdas
  }
  
  rels <- sn$rec_rels[, c("p1", "p2")]
  dd <- sn$dd
  prop0 <- sn$prop0
  const <- sn$const
  
  NRel <- nrow(rels)
  N_vec <- 1:N # used to speed up sampling in loop
  
  n_stage <- 4 # (U, E, A, Tc, To)
  
  infs <- matrix(F, nrow = N,  ncol = n_stage, dimnames = list(NULL, c("U", "E", "A", "T")))
  w <- order(lambdas, decreasing = T)[1:n_infs] # choose initial infections
  infs[w, "A"] <- T
  repeat_infs <- rep(0,N) # vector used to calculate total number of infections for each infectee
  
  strains <- rep(0,N)
  strains[w] <- rep(1:n_strain, n_infs) # determine which strain infects each infectee
  
  Ninf <- matrix(0, nrow = n_strain, ncol = n_stage,  dimnames = list(NULL, c("U", "E", "A", "T"))) # number of infs at start
  Ninf[ ,"A"] <- n_infs # introduce infections
  Ninf0 <- Ninf
  
  strain_vec <- sample.int(n_strain, size = max.iter, replace = T, prob = beta) # pre-generate strains resulting from infection
  symptoms_vec <- sample(c("E", "A"), size = max.iter, replace = T, prob = c(psi, 1 - psi)) # pre-generate probs of developing symptoms
  pres_vec <- sample.int(n_pres, size = max.iter, replace = T, prob = prop_pres) # pregenerate prescripton proportions
  treat_success_vec <- runif(n = max.iter) # pre-generate probs of developing symptoms
  # # can reduce number generated here if need be for comp time
  
  # create objects to record number of replacement infections vs new infections
  replacement_event <- matrix(0, nrow = n_strain, ncol = n_strain, 
                              dimnames = list(paste0("from.", 1:n_strain), paste0("to.",  1:n_strain)))# record which strains replace each other
  new_inf_event <- rep(0, n_strain)  # record which strains are causing new infections
  log_infs <- matrix(data = NA, nrow = max.iter, ncol = 9, dimnames = list(NULL, c("time", "c1", "c2", "strain", "p", "infector", "infector_ninf","infector_stage", "p_ninf"))) 
  
  # infs.t <- array(data = 0, dim = c(n_strain, n_stage, max.iter), dimnames = list(NULL, c("U", "E", "A", "T"), NULL)) # storing number of infs with each strain over time
  
  
  
  
  
  rinf <- NRel * sum(beta)  # R^_i(t_b) # set initial rate of infs = 0 (as implemented after burn-in)
  rsymp <- 0 # total rate of developing symptoms
  rtreat <- 0 # set initial rate of seeking treatment = 0 (as implemented after burn-in)
  rscreen <-  eta*sum(n_infs) # increase screening rate by eta
  rrec <- sum(nualpha[strains[w]]) # increase recovery rate by nu * alpha
  rcure <- 0 # set initial rate of cure (or failure) = 0 (as implemented after burn-in)
  
  step <- 1
  ext <- FALSE
  ext_day <- t*365
  
  while(time < time.window[2] & step <= max.iter){
    if (sum(Ninf)==0) break
    # print(c(radd, rdel, rinf, rrec))
    prob <- c(rinf, rsymp, rtreat, rrec, rscreen, rcure) # now a vector of length 6
    # names(prob) <- c('rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure')
    # if(any(prob < -1e-10)) browser()
     if(any(prob < 0)) prob <- pmax(0, prob)
    # if(sum(infs) != sum(Ninf)) browser()
    # if(length(prob) != 6) browser()
    # names(prob) <- c('rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure')
    # if(abs(sum(prob[c('rsymp', 'rtreat', 'rscreen', 'rcure')] - colSums(Ninf)*c(sigma, mu, eta, rho))) > 0.1)  browser()
    # if(abs(prob['rrec'] - sum(nualpha * Ninf[, "A"])) > 0.1) browser()
    # # prob[c('rsymp', 'rtreat', 'rscreen', 'rcure')]
    # colSums(Ninf)*c(sigma, mu, eta, rho)
    # infs.t[,,step -4:0]
    # if((time > burn.in) & (abs(sum(colSums(infs[,-1]) - Ninf )) > 1e-10)) browser()
    # colSums(infs)
    # Ninf
    e <- sample(c('rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure'), size = 1, prob = prob)
    
    # generate event time
    tau <- rexp(n = 1, rate = sum(prob))
    # update time
    time <- time + tau
    #store time
    times[step] <- time
    
    
    if (e=='rinf') {
      strain <- strain_vec[step]
      if(strain == 0) browser()
      # infection with strain 
      # print(paste0("rinf = ", rinf ))
      w <- sample.int(NRel, 1) # choose one couple
      
      w_inf_status <- strains[rels[w, c("p1", "p2")]] 
      
      if( !any(w_inf_status == strain)) next # neither are infected with strain
      if( all(w_inf_status > 0 )) next # both are infected
      if( any(infs[rels[w, c("p1", "p2")],"T"])) next # in treatment so not infectious
      
      infectee <- rels[w, which(w_inf_status != strain)]
      infector <- rels[w, which(w_inf_status == strain)]
      if(length(strains[infectee]) ==0) browser()
      
      # if(strains[infectee] > 0) { # if infectee partner is currently infected with other strain
      #   old_strain <- strains[infectee] # calculate which strain prev infection was
      #   Ninf[old_strain, ] <- Ninf[old_strain, ] - infs[infectee, ] # decrease total number of infections with prev strain
      #   infs[infectee, ] <- F # remove previous infection
      #   replacement_event[old_strain, strain] <-  replacement_event[old_strain, strain] + 1 # increase ticker
      # } else { # if infectee was not previously infected with another
      #   rsymp <- rsymp + sigma # increase overall rate of leaving incubation stage
      #   new_inf_event[strain] <- new_inf_event[strain] + 1 # increase ticker
      # }
      
      rsymp <- rsymp + sigma # increase overall rate of leaving incubation stage
      new_inf_event[strain] <- new_inf_event[strain] + 1 # increase ticker
      infs[infectee, "U"] <- T # record infected status for infectee
      strains[infectee] <- strain
      Ninf[strain, "U"] <- Ninf[strain, "U"] + 1 # increase total infs  with strain by 1
      repeat_infs[infectee] <- repeat_infs[infectee] + 1
      log_infs[step, c("time", "c1", "c2", "strain", "p", "infector", "infector_ninf","infector_stage", "p_ninf")] <- c(time, 1, 2, strain, infectee, infector, repeat_infs[infector], which(infs[infector, ])+1, repeat_infs[infectee])
    }      
    
    else if (e == 'rsymp') {
      # print(paste0("rrec = ", rrec , " vs. ", sum(infs)*rho))
      if(sum(Ninf[,"U"]) == 1) w <- N_vec[infs[,"U"]]
      else w <- sample(x = N_vec[infs[,"U"]], size = 1) # equal probability of leaving incubation stage for all new infections
      
      # if the above line is slow can think re-parametrise as per Rels
      strain <- strains[w]
      if(strain == 0) browser()
      
      infs[w, c("U", symptoms_vec[step])] <- c(F, T) # add infection to either symptomatic or asymptomatic compartment
      Ninf[strain, c("U", symptoms_vec[step])] <- Ninf[strain, c("U", symptoms_vec[step])] + c(-1, 1) # remove infection from U in correct strain
      
      rsymp <- rsymp - sigma # decrease overall rate of leaving incubation stage
      
      if (symptoms_vec[step] == "E") {
        rtreat <- rtreat + mu # if develop symptoms, increase treatment rate by mu
        log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 2, 3, strain, w, repeat_infs[w])
      }
      else if (symptoms_vec[step] == "A") { # if remains asymptomatic
        rscreen <- rscreen + eta # increase screening rate by eta
        rrec <- rrec + nualpha[strain] # increase recovery rate by nu * alpha
        log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 2, 4, strain, w, repeat_infs[w])
      }
    }
    
    else if (e == 'rtreat') {
      # print(paste0("rrec = ", rrec , " vs. ", sum(infs)*rho))
      if(length(N_vec[infs[,"E"]]) ==0) browser()
      if(sum(Ninf[,"E"]) == 1) w <- N_vec[infs[,"E"]]
      else w <- sample(x = N_vec[infs[,"E"]], size = 1) # equal probability of treatment seeking for all symptomatic carriers
      
      strain <- strains[w]
      if(strain == 0) browser()
      # if the above line is slow can think re-parametrise as per Rels
      # print(paste0("infs = ", colSums(infs[,1,])))
      # print(paste0("NInf = ", Ninf))
      
      infs[w, c("E", "T")] <- c(F, T) # remove infection from symptomatic stage, add infection to treatment stage
      Ninf[strain, c("E", "T")] <- Ninf[strain, c("E", "T")] + c(-1, 1) # remove infection from symptomatic stage, add infection to treatment stage
      
      rtreat <- rtreat - mu # decrease overall rate of seeking treatment by mu
      rcure <- rcure + rho # increase the overall rate of cure by rho
      log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf") ] <- c(time, 3, 5, strain, w, repeat_infs[w])
    }
    
    else if (e == 'rscreen') {
      # print(paste0("rrec = ", rrec , " vs. ", sum(infs)*rho))
      if(sum(Ninf[,"A"]) == 1) w <- N_vec[infs[,"A"]]
      else w <- sample(x = N_vec[infs[,"A"]], size = 1) # equal probability of being screened for all asymptomatic carriers
      strain <- strains[w]
      if(strain == 0) browser()
      infs[w, c("A", "T")] <- c(F, T) # remove infection from asymptomatic stage, add infection to treatment stage
      Ninf[strain, c("A", "T")] <- Ninf[strain, c("A", "T")] + c(-1, 1) # remove infection from asymptomatic stage, add infection to treatment stage
      
      rrec <- rrec - nualpha[strain] # decrease overall rate of recovery by nu * alpha
      rscreen <- rscreen - eta # decrease screening rate by eta
      rcure <- rcure + rho # increase the overall rate of cure by rho
      log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 4, 5, strain, w, repeat_infs[w])
    }
    
    else if (e == 'rrec') {
      # print(paste0("rrec = ", rrec , " vs. ", sum(infs)*rho))
      strain <- sample.int(n = n_strain, size = 1, prob = Ninf[, "A"] * nualpha) # sample which strain with prob A_s * nu * alpha
      if(Ninf[strain, "A"] == 1) w <- N_vec[strains == strain & infs[, "A"]]
      else w <- sample(N_vec[strains == strain & infs[, "A"]], 1) # equal probability of recovery for all asymptomatic carriers with that strain
      
      if(strain == 0) browser()
      Ninf[strain, "A"] <- Ninf[strain, "A"] - 1 # decrease total infs by 1 in correct strain
      infs[w, "A"] <- F # remove infection from asymptomatic stage and set strain to 0
      strains[w] <- 0
      
      # if the above line is slow can think re-parametrise as per Rels
      
      rrec <- rrec - nualpha[strain] # decrease overall rate of recovery
      rscreen <- rscreen - eta # also decrease screening rate by eta
      log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 4, 1, strain, w, repeat_infs[w])
    }
    
    else if (e == 'rcure') {
      # print(paste0("rrec = ", rrec , " vs. ", sum(infs)*rho))
      # probw <- rowSums(infs[, , "T", drop = F])
      # if(sum(probw) < 1) browser()
      # print(sum(probw))
      if(sum(Ninf[,"T"]) == 1) w <- N_vec[infs[,"T"]]
      else w <- sample(x = N_vec[infs[,"T"]], size = 1) # equal probability of cure for all treated infections
      strain <- strains[w]
      if(strain == 0) browser()
      
      
      Ninf[strain, "T"] <- Ninf[strain, "T"] - 1 # decrease total infs by 1 in correct strain
      infs[w, "T"] <- F # remove infection from treatment stage and set strain to 0
      rcure <- rcure - rho # decrease the overall rate of cure by rho
      
      pres <- pres_vec[step]
      if(strain == 0) browser()
      if(treat_success_vec[step] < efficacy[strain, pres]) { # treatment success
        strains[w] <- 0
        log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 5, 1, strain, w, repeat_infs[w])
      } else { # treatment failure
        Ninf[strain, "A"] <- Ninf[strain, "A"] + 1 # increase asymptomatics by 1 in correct strain
        infs[w,  "A"] <- T #  asymptomatic = 1 in correct strain
        rscreen <- rscreen + eta # increase screening rate by eta
        rrec <- rrec + nualpha[strain] # increase recovery rate by nu * alpha
        log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 5, 4, strain, w, repeat_infs[w])
      } 
      
    }

    # infs.t[, , step] <- Ninf # record number of infs of each strain
    if(any(Ninf < 0)) browser() # check
    step <- step + 1
  }
  
  print(Sys.time() - start)
  
  # remove unused storage
  Nstep <- step - 1
  times <- c(0,times[1:Nstep])
  log_infs <- log_infs[!is.na(log_infs[,1]), , drop=F]

  # infs.t <- infs.t[, , 1:Nstep, drop=F]
  
  
  log_infs  <- as.data.frame(log_infs)
  if( nrow(log_infs) > 0 ) {
    log_infs$c1 <- factor(x = log_infs$c1, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$c2 <- factor(x = log_infs$c2, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$infector_stage <- factor(x = log_infs$infector_stage, levels = 2:4, labels = c( "U", "E", "A"), exclude = NA)
    log_infs$infector_full <- paste0(log_infs$infector, "_", log_infs$infector_ninf)
    log_infs$p_full <- paste0(log_infs$p, "_", log_infs$p_ninf)
    }
  sn <- list(inputs = inputs, const = const, 
             # rels = rels, 
             times = times, 
             # infs.t = infs.t, 
             # inf_event = rbind(replacement_event, "new" = new_inf_event), 
             log_infs = log_infs,
             degree = rowSums(sn$res) + colSums(sn$res),
             dd = dd, prop0 = prop0,
             comp_time = Sys.time () - start)
  # dists <- extract_dists(sn)
  sn$dists <- NULL
  sn$R0 <- rep(0,3)
  if (sum(Ninf) == 0 ) {
    ext_day <- floor(max(log_infs$time*365))
    if (ext_day < 365) ext <- TRUE
  }
  sn$ext = ext
  sn$ext_day = ext_day
  sn$total_treated <-  sum(((log_infs$c1==5) | (log_infs$c1 == "T")) & log_infs$time <1 )

  return(sn)
}
