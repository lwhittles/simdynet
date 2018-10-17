#' Simulation of dynamic network
#'
#' @param N network size
#' @param gamma exponent of power law distribution
#' @param k0 parameter controlling proportion of unpartnered
#' @param phi parameter controlling relationship turnaround rate
#' @param n_infs number of initially infected individuals
#' @param beta rate of transmission per partnership
#' @param psi proportion of infections becoming symptomatic
#' @param sigma reciprocal of duration of incubation period
#' @param nu reciprocal of duration of carriage
#' @param eta rate of screening when asymptomatic
#' @param mu reciprocal of time to seeking treatment
#' @param rho reciprocal of time to recovery following treatment
#' @param n_vax number of vaccinated individuals
#' @param vax_strat vaccination strategy, can be either 'target' or 'random'
#' @param t duration of outbreak simulation
#' @param max.iter maximum number of Gillespie iterations
#' @param burn.in burn-in period
#' @param record whether or not to record degree distribution
#' @param record_term observation window of the degree distribution
#' @param record_lengths whether or not to record partnership durations
#'
#' @export
sim_dynamic_sn <- function (N = 1e4, 
                            gamma = 1.8, k0 = 0.5,  phi = N, 
                            n_infs = 0,
                            beta = 1, psi = 1,  sigma = 1, 
                            nu = 1, eta = 1, mu = 1,rho = 0, n_vax = 0, vax_strat=NA,
                            t = 1, max.iter = 1e6, burn.in = 0,
                            record=TRUE, record_term = 1, record_lengths = TRUE) {
  
  start <- Sys.time ()
  
  inputs <- list(N = N, gamma = gamma, k0 = k0, phi = phi, 
                 n_infs = n_infs,
                 beta = beta, psi = psi, sigma = sigma, 
                 nu = nu, eta = eta, mu =  mu, rho = rho, n_vax = n_vax, vax_strat = vax_strat,
                 t = t, max.iter = max.iter, burn.in = burn.in,
                 record_term = record_term)
  
  lambdas <- rexp(n = N, rate = 1) # sample N lambda ~ Exp(1)
  
  
  kmax <- optimise(f = function(z) { #  optimise kmax
    abs(calculate_g(x = Inf, gamma = gamma, k0 = k0, N = N, kmax = z)$g - 1) 
  }, 
  interval = c(1, N))
  kmax <- floor(kmax$minimum)

 

  # burn.in is in days not iterations

  time <- 0
  tau <- 0 # initialise
  times <- rep(0, max.iter) # storing times at which events happen
  
  
  # generate partnerships at equilibrium M0
  NS <- N * (N - 1) / 2 # total possible partnerships
  q0 <- 1 / (1 + phi) # equilibrium probability of partnership
  NRel <- rbinom(n=1, size = N ^ 2 / 2, prob = q0) # note the size is  N^2/2 because will do rejection sampling on partnerhips with self
  
  rels0 <- matrix(sample.int(n = N, size = 2 * NRel, replace = T), ncol = 2) # sample partnerships (with equal probability)
  rels0 <- rels0[rels0[, 1 ] != rels0[, 2], ] # remove duplicates (rejection sampling)
  rels0 <- t(apply(X = rels0, MARGIN = 1, FUN = sort)) # put in right order (also transposes)
  NRel <- nrow(rels0) # update NRel
  NRel0 <- NRel # record initial number of rels for output
  
  gs <- calculate_g(x = lambdas, gamma = gamma, k0 = k0, kmax = kmax, N = N) # calculate g(lamda)
  g_fac <- 1 + 1 / phi

  gl <- gs$g * sqrt(g_fac) # calculations are vectorised for speed
  const <- gs$c
  
  
  radd_mat <- -log(g_fac - gl %o% gl) 
  radd_mat <- pmax(radd_mat, 1/50/phi) # set maximum partnership length
  diag(radd_mat) <- 0 # no partnerships with self
  
  radd <- sum(radd_mat)/2 #total rate of adding partnerships

  degree_vec <- rep(0, N) # storing degree distribution if 'record' is enabled
  
  length_mat <- NULL # storing length of each partnership if enabled
  if (record_lengths) length_mat <- matrix(NA, nrow = max.iter/2, ncol = 4, dimnames = list(NULL, c('p1', 'p2', 'start', 'end')))
  
  
  # matrix set to be large enough to store double expected rels (to be safe)
  rels <- matrix(NA, nrow = NS * q0 * 2, ncol = 4, dimnames = list(NULL, c('p1', 'p2', 'rdel', 't')))
  rels[1:NRel, c("p1", "p2") ] <- rels0  # which partnerships currently exist?
  rels[1:NRel, "rdel"] <- radd_mat[rels0] * phi # add column of rdel for current partnerships
  rels[1:NRel, "t"]  <- -rexp(n = NRel, rate = rels[1:NRel, "rdel"]) # add col of time at which each partnership formed
  
  rdel <- sum(rels[, "rdel"], na.rm = T) #total rate of deleting partnerships             
  total_rel <- 0 # counter of partnerships ever in network, count as they break up
  
  
  # pre-generate matrix of proposed partnerships (rather than doing in loop)
  p <- sample.int(N^2, size = max.iter, replace = T, prob = radd_mat )
  W <- cbind((p - 1) %% N + 1 , (p - 1) %/% N + 1) # populate matrix based on p
  
  n_stage <- 4 # (U, E, A, Tc, To)
   
  infs <- matrix(F, nrow = N,  ncol = n_stage, dimnames = list(NULL, c("U", "E", "A", "T")))
  repeat_infs <- rep(0,N) # vector used to calculate total number of infections for each infectee
  N_vec <- 1:N # used to speed up sampling in loop
  
  strains <- rep(0,N) # storage vector for wildtype or vaccine infection
  
  Ninf <- rep(0, n_stage)
  names(Ninf) <- c("U", "E", "A", "T")
  
  
  symptoms_vec <- sample(c("E", "A"), size = max.iter, replace = T, prob = c(psi, 1 - psi)) # pre-generate probs of developing symptoms
  
  new_inf_event <- 0  # record total infs
  
  log_infs <- matrix(data = NA, nrow = max.iter, ncol = 8, dimnames = list(NULL, c("time", "c1", "c2", "p", "infector", "infector_ninf","infector_stage", "p_ninf"))) 
    # record c as (S=1, U=2, E=3, A=4, T=5)
  rinf <- 0 # set initial rate of infs = 0 (as implemented after burn-in)
  rsymp <- 0 # set initial rate of leaving incubation stage = 0 (as implemented after burn-in)
  rtreat <- 0 # set initial rate of seeking treatment = 0 (as implemented after burn-in)
  rrec <-  0  # set initial rate of recovery = 0 (as implemented after burn-in)
  rscreen <- 0 # set initial rate of screening = 0 (as implemented after burn-in)
  rcure <- 0 # set initial rate of cure = 0 (as implemented after burn-in)
  step <- 1 # initialize
  
  breaktime <- ifelse(record, record_term, burn.in)
  
  while(time < t & step <= max.iter){
    
    if (time < burn.in) rinf <- 0 # total rate of infections = 0 (as implemented after burn-in)
    
    if ((time >= burn.in) & (time - tau <= burn.in)) {      # start infections/recoveries
      
      ranked_lambdas <- order(lambdas, decreasing = T)
      w <- ranked_lambdas[1:n_infs] # choose initial infections
      infs[w, "A"] <- T # make them asymptomatic
      Ninf["A"] <- n_infs # record number
      if (n_infs > 0) strains[w] <- 1 
      
      if (n_vax > 0) { # if vaccinating
        
        if(vax_strat =="target") { # with targeted strategy
          v <- ranked_lambdas[n_infs + 1:n_vax]
        } else if (vax_strat == "random") { # with randomised strategy
          v <- sample(x = ranked_lambdas[-(1:n_infs)], size = n_vax)
        }
         strains[v] <- 99 ### vaccinate
      }
      
      if (record) { # if recording relationship history
        tab <- table(rels[, c("p1", "p2")])
        degree_vec[as.numeric(names(tab))] <- tab
      }
      
      rscreen <- eta*sum(n_infs) # increase screening rate by eta
      rrec <- nu*n_infs # increase recovery rate by nu 
      rinf <- NRel * beta  # R^_i(t_b)
    }

    if (time > breaktime & sum(infs) == 0) break # stop if outbreak dies out

    prob <- c(radd, rdel, rinf, rsymp, rtreat, rrec, rscreen, rcure) # now a vector of length 7 
    
    e <- sample(c('radd', 'rdel', 'rinf', 'rsymp', 'rtreat', 'rrec', 'rscreen', 'rcure'), size = 1, prob = prob)
    
    tau <- rexp(n = 1, rate = sum(prob)) # generate event time
    time <- time + tau # update time
    times[step] <- time    #store time
    
    if (e =='radd') {  #add partnership
      
      w <- W[step,] # choose couple
      if (w[1] > w[2]) w <- w[c(2, 1)] # order sequentially
      if (NRel > 0 & {any({rels[1:NRel, 1] == w[1]} & {rels[1:NRel, 2] == w[2]})}) next # partnership already exists
        
      r <- radd_mat[w[1], w[2]] * phi # calculate rate of breakup

      rels[NRel + 1, c('p1', 'p2', 'rdel', 't')] <- c(w, r, time) # store details of new partnership
      NRel <- NRel + 1  # increase overall number of partnerships
      rdel <- rdel + r  # increase overall rate of breakup by r
      rinf <- rinf + beta # increase overall rate of infection by beta
      if(record & time <= record_term) degree_vec[w] <- degree_vec[w] + 1 # keep track of degree dist
    } 
    
    else if (e=='rdel') { #remove partnership

      w <- sample.int(NRel, 1, prob = rels[1:NRel, 'rdel']) # choose couple
      
      if(record_lengths) { # if we are recording partnership lengths
        total_rel <- total_rel + 1 # increase total number of completed partnerships by 1
        length_mat[total_rel, c('p1', 'p2', 'start')] <- rels[w, c('p1', 'p2', 't')] # store partners and start time
        length_mat[total_rel, 'end'] <- time  # record current time as partnership ends
      }
      
      rdel <- rdel - rels[w, 'rdel']  # reduce overall rate of partnership breakup
      rinf <- rinf - sum(beta) # reduce overall rate of infection
      rels[w, ] <- rels[NRel, ] # move last rel recorded to replace deleted rel
      rels[NRel, ] <- NA 
      NRel <- NRel - 1     # decrease total rels by 1
      
    } 
    
    else if (e=='rinf') { #infection event

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
  
  #print(Sys.time () - start)
  
  # remove unused storage
  log_infs <- log_infs[!is.na(log_infs[,1]), ,drop=F]
  log_infs  <- as.data.frame(log_infs)

  if( nrow(log_infs) > 0) {
    log_infs$c1 <- factor(x = log_infs$c1, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$c2 <- factor(x = log_infs$c2, levels = 1:5, labels = c("S", "U", "E", "A", "T"))
    log_infs$infector_stage <- factor(x = log_infs$infector_stage, levels = 2:4, labels = c( "U", "E", "A"), exclude = NA)
    log_infs$infector_full <- paste0(log_infs$infector, "_", log_infs$infector_ninf)
    log_infs$p_full <- paste0(log_infs$p, "_", log_infs$p_ninf)
    }
  


  dd <- calc_dd(degree_vec, N)
  const <- (1-gamma) / ((kmax * (1 - dd$prop0)) ^ (1-gamma) - k0 ^ (1-gamma))
  
  sn <- list(inputs = inputs,  kmax = kmax, const = const, 
             log_infs = log_infs,
             dd = dd$dd, prop0 = dd$prop0,
             comp_time = Sys.time () - start)
  
  if(record_lengths) {
    length_mat <- length_mat[1:total_rel, ] # remove unused storage
    # length_mat <- cbind(length_mat, "length" = length_mat[,"end"] - length_mat[,"start"])
    sn <- c(sn, "log_rels" = length_mat)
  }
  
  if(record) {
    deg_curr <- rep(0, N)
    deg_curr_table <- table(rels[, c("p1", "p2")])
    deg_curr[as.numeric(names(deg_curr_table))] <- as.numeric(deg_curr_table)
    prop0_curr <- sum(deg_curr == 0)/N
    net_out <- list(degree = degree_vec, 
                     degree_curr = deg_curr, prop0_curr = prop0_curr)
    sn <- c(sn, net_out)
  }
  
  return(sn)
  
}

