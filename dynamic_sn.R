source("static_sn.R")
source("calc_dd.R")


sim_dynamic_sn <- function (N, 
                            gamma, k0,  phi, 
                            n_strain = 1,  n_infs = 0,
                            beta = 1, psi = 1,  sigma = 1, alpha = 1,
                            nu = 1, eta = 1, mu = 1,rho = 0, n_vax = 0, vax_strat=NA,
                            t, max.iter, burn.in = t,
                            record=FALSE, record_term = 1, record_lengths = FALSE) {
  
  start <- Sys.time ()
  
  inputs <- list(N = N, gamma = gamma, k0 = k0, phi = phi, 
                 n_strain = n_strain,  n_infs = n_infs,
                 beta = beta, psi = psi, sigma = sigma, alpha = alpha,
                 nu = nu, eta = eta, mu =  mu, rho = rho, n_vax = n_vax, vax_strat = vax_strat,
                 t = t, max.iter = max.iter, burn.in = burn.in,
                 record_term = record_term)
  
  lambdas <- rexp(n = N, rate = 1) # sample N lambda ~ Exp(1)
  
  
  kmax <- optimise(f = function(z) {
    abs(calculate_g(x = Inf, gamma = gamma, k0 = k0, N = N, kmax = z)$g - 1) 
  }, 
  interval = c(1, N))
  kmax <- floor(kmax$minimum)

 
  # n_infs is a vector of length n_strain
  # burn.in is now in days not iterations
  # infs is now a N x n_strain matrix of 0s and 1s
  
  
  #input checks
  if (n_strain != length(beta)) stop("number of strains and supplied betas do not match")
  if (n_strain != length(alpha)) stop("number of strains and supplied alphas do not match")
  if (n_strain != length(n_infs)) stop("number of strains and supplied infs do not match")
 
  
  nualpha <- nu * alpha

  
  time.window <- c(0, t)
  time <- time.window[1]
  tau <- 0 # initialise
  times <- rep(0, max.iter) # storing times at which events happen
  
  
  # generate partnerships at equilibrium M0
  NS <- N * (N - 1) / 2
  q0 <- 1 / (1 + phi)
  NRel <- rbinom(n=1, size = N ^ 2 / 2, prob = q0) # note the size is  N^2/2 because will do rejection sampling
  
  rels0 <- matrix(sample.int(n = N, size = 2 * NRel, replace = T), ncol = 2) # sample partnerships (with equal probability)
  rels0 <- rels0[rels0[, 1 ] != rels0[, 2], ] # remove duplicates (rejection sampling)
  rels0 <- t(apply(X = rels0, MARGIN = 1, FUN = sort)) # put in right order (also transposes)
  NRel <- nrow(rels0) # update NRel
  NRel0 <- NRel
  
  gs <- calculate_g(x = lambdas, gamma = gamma, k0 = k0, kmax = kmax, N = N) # calculate g(lamda)
  g_fac <- 1 + 1 / phi

  gl <- gs$g * sqrt(g_fac) # doing calcs on vector rather than on large matrix
  const <- gs$c
  
  
  radd_mat <- -log(g_fac - gl %o% gl) # work on this to speed up?
  radd_mat <- pmax(radd_mat, 1/50/phi) # if not optimising phi force min a(,) to be 0
  diag(radd_mat) <- 0 # so if in relationship it is for the rest of their life
  
  radd <- sum(radd_mat)/2 #total rate of adding relationships

  degree_vec <- rep(0, N) # storing degree distribution if 'record' is enabled
  
  length_mat <- NULL # storing length of each relationship if enabled
  if (record_lengths) length_mat <- matrix(NA, nrow = max.iter/2, ncol = 4, dimnames = list(NULL, c('p1', 'p2', 'start', 'end')))
  
  
  # refine matrix size to be within 2sds of mean NS * q0
  rels <- matrix(NA, nrow = NS * q0 * 2, ncol = 4, dimnames = list(NULL, c('p1', 'p2', 'rdel', 't')))
  rels[1:NRel, c("p1", "p2") ] <- rels0  # which relationships currently exist?
  rels[1:NRel, "rdel"] <- radd_mat[rels0] * phi # add column of rdel for current relationships
  rels[1:NRel, "t"]  <- -rexp(n = NRel, rate = rels[1:NRel, "rdel"]) # add col of time at which each relationship formed
  
  rdel <- sum(rels[, "rdel"], na.rm = T) #total rate of deleting relationships             
    
  rels.t <- rep(0, max.iter) # storing number of relationships over time
  total_rel <- 0 # counter of relationships ever in network, count as they break up
  
  
  # pre-generate matrix of proposed relationships (rather than doing in loop)
  p <- sample.int(N^2, size = max.iter, replace = T, prob = radd_mat )
  W <- cbind((p - 1) %% N + 1 , (p - 1) %/% N + 1) # populate matrix based on p
  
  n_stage <- 4 # (U, E, A, Tc, To)
   
  infs <- matrix(F, nrow = N,  ncol = n_stage, dimnames = list(NULL, c("U", "E", "A", "T")))
  repeat_infs <- rep(0,N) # vector used to calculate total number of infections for each infectee
  N_vec <- 1:N # used to speed up sampling in loop
  
  strains <- rep(0,N)
  
  Ninf <- matrix(0, nrow = n_strain, ncol = n_stage,  dimnames = list(NULL, c("U", "E", "A", "T"))) # number of infs at start
  Ninf0 <- Ninf
  
  
  strain_vec <- sample.int(n_strain, size = max.iter, replace = T, prob = beta) # pre-generate strains resulting from infection
  symptoms_vec <- sample(c("E", "A"), size = max.iter, replace = T, prob = c(psi, 1 - psi)) # pre-generate probs of developing symptoms
   treat_success_vec <- runif(n = max.iter) # pre-generate probs of developing symptoms
  # # can reduce number generated here if need be for comp time
  
  # create objects to record number of replacement infections vs new infections
  
  replacement_event <- matrix(0, nrow = n_strain, ncol = n_strain, 
                               dimnames = list(paste0("from.", 1:n_strain), paste0("to.",  1:n_strain)))# record which strains replace each other
  new_inf_event <- rep(0, n_strain)  # record which strains are causing new infections
  
   # infs.t <- array(data = 0, dim = c(n_strain, n_stage, max.iter), dimnames = list(NULL, c("U", "E", "A", "T"), NULL)) # storing number of infs with each strain over time
  log_infs <- matrix(data = NA, nrow = max.iter, ncol = 9, dimnames = list(NULL, c("time", "c1", "c2", "strain", "p", "infector", "infector_ninf","infector_stage", "p_ninf"))) 
    # record c as (S=1, U=2, E=3, A=4, T=5)
  rinf <- 0 # set initial rate of infs = 0 (as implemented after burn-in)
  rsymp <- 0 # set initial rate of leaving incubation stage = 0 (as implemented after burn-in)
  rtreat <- 0 # set initial rate of seeking treatment = 0 (as implemented after burn-in)
  rrec <-  0  # set initial rate of recovery = 0 (as implemented after burn-in)
  rscreen <- 0 
  rcure <- 0
  step <- 1
  ext <- FALSE 
  ext_day <- t*365
  v <- NA
  
  breaktime <- ifelse(record, record_term, burn.in)
  
  while(time < time.window[2] & step <= max.iter){
    if (time < burn.in) rinf <- 0 # total rate of infections = 0 (as implemented after burn-in)
    if ((time > burn.in) & (time - tau < burn.in)) {      # start infections/recoveries
      ranked_lambdas <- order(lambdas, decreasing = T)
      w <- ranked_lambdas[1:n_infs] # choose initial infections
      infs[w, "A"] <- T
      Ninf[ ,"A"] <- n_infs # introduce infections
      if (n_infs>0) strains[w] <- rep(1:n_strain, n_infs) # determine which strain infects each infectee
      if (n_vax > 0) {
        if(vax_strat =="target") v <- ranked_lambdas[n_infs + 1:n_vax]
        else if (vax_strat == "random") v <- sample(x = ranked_lambdas[-(1:n_infs)], size = n_vax)
         strains[v] <- 99 ### vaccinate
      }
      if (record) {
        tab <- table(rels[, 1:2])
        degree_vec[as.numeric(names(tab))] <- tab
      }
      rscreen <- eta*sum(n_infs) # increase screening rate by eta
      rrec <- sum(nualpha[strains[w]]) # increase recovery rate by nu * alpha
      rinf <- NRel * sum(beta)  # R^_i(t_b)
    }

    if (time > breaktime & sum(infs)==0) break
    # print(c(radd, rdel, rinf, rrec))
    prob <- c(radd, rdel, rinf, rsymp, rtreat, rrec, rscreen, rcure) # now a vector of length 7 
    # names(prob) <- c('radd', 'rdel', 'rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure')
    if(any(prob < -1e-10)) browser()
    if(any(prob < 0)) prob <- pmax(0, prob)
    if(length(prob) != 8) browser()
    # if(abs(sum(prob[-(1:3)] - colSums(Ninf)[c("U", "E", "A", "A", "T")]*c(sigma, mu, nu, eta, rho))) > 0.1)  browser()
    # prob[-(1:3)]
    # colSums(Ninf)[c("U", "E", "A", "A", "T")]*c(sigma, mu, nu, eta, rho)
    # infs.t[,,step -2:0]
      # if((time > burn.in) & (abs(sum(colSums(infs) - colSums(Ninf) )) > 1e-10)) browser()
    # colSums(infs)
    # Ninf
    e <- sample(c('radd', 'rdel', 'rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure'), size = 1, prob = prob)
    
    # generate event time
    tau <- rexp(n = 1, rate = sum(prob))
    # update time
    time <- time + tau
    #store time
    times[step] <- time
    
    if (e =='radd') {
      # print(paste0("radd = ", radd ))
      #add relation
      w <- W[step,]
      # if (w[1]==w[2]) next # cannot have relationship with self (no longer needed)
      # if (!sex_poss[type[w[1]], type[w[2]]]) next # relation not possible (no longer needed)
      if (w[1] > w[2]) w <- w[c(2, 1)]
      if (NRel > 0 & {any({rels[1:NRel, 1] == w[1]} & {rels[1:NRel, 2] == w[2]})}) next #relation already exists
        
      r <- radd_mat[w[1], w[2]] * phi

      rels[NRel + 1, ] <- c(w, r, time)
      NRel <- NRel + 1
      rdel <- rdel + r
      rinf <- rinf + sum(beta) # increase overall rate of infection
      if(record & time <= record_term) degree_vec[w] <- degree_vec[w] + 1
    } 
    
    else if (e=='rdel') {
      #remove relation
      # print(paste0("rdel = ", rdel, " vs. ", sum(rels[,3]) ))
      w <- sample.int(NRel, 1, prob = rels[1:NRel, 'rdel'])
      
      if(record_lengths) { # update this on breakup for simplicity
        total_rel <- total_rel + 1 # increase total number of completed relationships by 1
        length_mat[total_rel, c('p1', 'p2', 'start')] <- rels[w, c('p1', 'p2', 't')] # store partners and start time
        length_mat[total_rel, 'end'] <- time  # record current time as relationship ends
      }
      
      rdel <- rdel - rels[w, 'rdel']
      if(rdel < 0 ) browser()
      rinf <- rinf - sum(beta) # reduce overall rate of infection
      rels[w, ] <- rels[NRel, ] # move last rel recorded to replace deleted rel
      rels[NRel, ] <- NA
      NRel <- NRel - 1     # decrease total rels by 1
      
    } 
    
    
    if (e=='rinf') {
      strain <- strain_vec[step]
      if(strain == 0) browser()
      # infection with strain 
      # print(paste0("rinf = ", rinf ))
      w <- sample.int(NRel, 1) # choose one couple
      
      w_inf_status <- strains[rels[w, c("p1", "p2")]] 
      
      if( !any(w_inf_status == strain)) next # neither are infected with strain
      if( any(w_inf_status == 99)) next # vaccinated
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
      log_infs[step, c("time", "c1", "c2", "strain", "p","p_ninf")] <- c(time, 4, 5, strain, w, repeat_infs[w])
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
      
      if(strain == 0) browser()
        strains[w] <- 0
        log_infs[step, c("time", "c1", "c2", "strain", "p", "p_ninf")] <- c(time, 5, 1, strain, w, repeat_infs[w])
      
      
    }
    rels.t[step] <- NRel # record total number of rels
     # infs.t[, , step] <- Ninf # record number of infs of each strain
    if(any(Ninf < 0)) browser() # check
    
    # speed-up to remove extinct strains
    # if (time > burn.in) { # after infections are introduced
    #   if (any(colSums(infs.t[step - c(1,0),, drop=F]) == 1)) { # if a strain has just gone extinct
    #     beta[rowSums(Ninf) == 0] <- 0 # remove strain from circulation to speed up algorithm
    #     rinf <- NRel * sum(beta)  # recalculate rinf
    #     if (rinf > 0) strain_vec[(step+1):max.iter] <- sample.int(n_strain, size = max.iter - step, replace = T, prob = beta) # remove strain from future calcs
    #     # as if all strains have died out can just continue without infection
    #   }
    # }
    # 
    # records the degree distribution of the network over specified term
    #   if (times[step] > t - record_term) {
    #     if (times[step - 1] < t - record_term) {
    #       tab <- table(rels[, 1:2])
    #       degree_vec[as.numeric(names(tab))] <- tab
    #       rec_rels <- matrix(data = NA, nrow = max.iter, ncol = 2, dimnames = list(NULL, c("p1", "p2")))
    #       rec_rels[1:NRel, ] <- rels[1:NRel, c("p1", "p2")] # partnerships in existence
    #       NRel_rec <- NRel
    #     } 
    #     else if (times[step - 1] > t - record_term) {
    #       if( rels.t[step] > rels.t[step - 1]) { # only counts rels that have actually been added
    #         degree_vec[w] <- degree_vec[w] + 1
    #         NRel_rec <- NRel_rec + 1 # only ever increase at this point
    #         rec_rels[NRel_rec, ] <- w # add partners
    #     }
    #   }
    # }
     
    step <- step + 1
  }
  
  print(Sys.time () - start)
  
  # remove unused storage
  Nstep <- step - 1
  times <- c(0,times[1:Nstep])
  log_infs <- log_infs[!is.na(log_infs[,1]), ,drop=F]
  # times <- c(0, times)
  rels.t <- c(NRel0, rels.t[1:Nstep])
   # infs.t <- infs.t[, , 1:Nstep, drop=F]
  
  log_infs  <- as.data.frame(log_infs)
  # print(log_infs)
  if( nrow(log_infs) > 0) {
    log_infs$c1 <- factor(x = log_infs$c1, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$c2 <- factor(x = log_infs$c2, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$infector_stage <- factor(x = log_infs$infector_stage, levels = 2:4, labels = c( "U", "E", "A"), exclude = NA)
    log_infs$infector_full <- paste0(log_infs$infector, "_", log_infs$infector_ninf)
    log_infs$p_full <- paste0(log_infs$p, "_", log_infs$p_ninf)
    
    }
  
  
  

  dd <- calc_dd(degree_vec, N)
  const <- (1-gamma) / ((kmax * (1 - dd$prop0)) ^ (1-gamma) - k0 ^ (1-gamma))
  

  
  deg_curr <- rep(0, N)
  deg_curr_table <- table(rels[, c("p1", "p2")])
  deg_curr[as.numeric(names(deg_curr_table))] <- as.numeric(deg_curr_table)
  prop0_curr <- sum(deg_curr == 0)/N
  
  
  if(record_lengths) {
    length_mat <- length_mat[1:total_rel, ] # remove unused storage
    length_mat <- cbind(length_mat, "length" = length_mat[,"end"] - length_mat[,"start"])
  }
  
  
  sn <- list(inputs = inputs,  kmax = kmax, const = const, 
             log_infs = log_infs,
             vax_check = intersect(v, log_infs$p),
             # infs.t = infs.t, 
             dd = dd$dd, prop0 = dd$prop0,
             inf_event = rbind(replacement_event, "new" = new_inf_event), 
             comp_time = Sys.time () - start)
  
  
  
  if(record) {
    net_out <- list( 
      # rels=rels[1:NRel,c("p1", "p2") ], 
                     # rec_rels = rec_rels[1:NRel_rec, ],
                     # rel_lengths = length_mat[,"length"], 
                     # rels.t = rels.t, 
                     degree = degree_vec, 
                     degree_curr = deg_curr, prop0_curr = prop0_curr)
    sn <- c(sn, net_out)
  }
  
  sn$dists <- NA
  sn$R0 <- rep(0,3)
  
  if(sum(new_inf_event)>0) {
    # dists <- extract_dists(sn)
    # sn$dists <- dists
    sn$R0 <- rep(0,3)
  }

  if (sum(Ninf) == 0 ) {
    ext_day <- floor(max(log_infs$time*365))
    if (ext_day < 365) ext <- TRUE
  }
  
  sn$ext = ext
  sn$ext_day = ext_day
  
  sn$total_treated <- sum(((log_infs$c1==5) | (log_infs$c1 == "T")) & log_infs$time <1 )
  
  
  return(sn)
  
}

