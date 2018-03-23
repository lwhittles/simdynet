source("calc_dd.R")
source("calc_g.R")

sim_static_sn <- function(gamma, k0,  N=1e4, phi = 1e4) {
  inputs <- list(gamma = gamma, k0 = k0, N = N,  phi = phi)
  l  <- rexp(n = N, rate = 1) # draw N lambdas from Exp(1) dist
  
  kmax <- optimise(f = function(z) { #optimise kmax
    abs(calculate_g(x = Inf, gamma = gamma, k0 = k0, N = N, kmax = z)$g - 1) 
    }, 
    interval = c(1, N))
  kmax <- floor(kmax$minimum) 
  
  gs <- calculate_g(x = l, gamma = gamma, k0 = k0, N = N, kmax = kmax, mu = 1) # calculate g(lamda)
  gl <- gs$g 
  const <- gs$c
  Q <- outer(X = gl, Y = gl, FUN = '*') # calculate matrix Q
  min_q <- 1 - exp(-1/50/phi) / (1 + 1/phi)
  Q <- pmax(Q, min_q)   # set maximum partnership length
  
  U <- matrix(1, nrow=N, ncol=N) # matrix of U[0,1] draws in lower triangle
  U[row(U) < col(U)] <- runif(N * (N - 1) / 2)
  
  M <- matrix(0, nrow = N, ncol = N) # matrix of partnerships (if U < Q)
  M[U < Q] <- 1

  
  degree_vec <- colSums(M) + rowSums(M)
  dd <- calc_dd(degree_vec = degree_vec, N = N)
  

  const <- (1-gamma) / ((kmax * (1 - dd$prop0)) ^ (1-gamma) - k0 ^ (1-gamma))
  
  sn <- list(inputs = inputs, res = M, dd = dd$dd, prop0 = dd$prop0, const = const, kmax = kmax, lambdas = l)
  
  return(sn)
}
