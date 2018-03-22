# case where lambda ~ Exp(1)


source("general/project_colours.R")
source("sn/figs/plot_dd_function.R")
require(igraph)


calculate_g <- function(x, gamma, k0, N, kmax = N - 1, mu = 1) {

  const <- (1-gamma) / (kmax ^ (1-gamma) - k0 ^ (1-gamma))
  
  k <- function(x) { # k in terms of F(lambda)
    ((1-gamma)/const*x + k0^(1-gamma)) ^ (1/(1-gamma)) 
  }
  kdash <- integrate(f = k, lower = 0, upper = 1)
  # check_kdash <- ((const*k0^(1-gamma) + 1-gamma) * (k0^(1-gamma) + (1-gamma)/const) ^ (1/(1-gamma)) - const*k0^(2-gamma)) / (2-gamma)
  # print(kdash$value); print(check_kdash)
  # working fine, so go with simpler option in terms of code
  
  gbar <- sqrt(kdash$value / (N-1))
  Fl <- pexp(x, rate = mu) # since lambdas are dist Exp(mu)
  g <- k(x = Fl) / (N-1) / gbar
  
  return(list(g=g, c=const))
}

# # # try with gamma = -1
 calculate_g_gamma_1 <- function(x, k0, N, kmax = N - 1) {
  gamma <- 1
  bmin <- k0 / (N - 1)
  bmax <- kmax / (N - 1)
  const <- 1 / log(kmax / k0)
  R <- 1 - exp(-x)
  # Ngbar2 <- ((const * k0 ^ beta + beta) * (beta / const + k0 ^ beta) ^ (1 / beta) - const * k0 ^ (beta + 1)) / (beta + 1)
  # gbar <- sqrt(k0 * const * (exp(1/const) - 1) / N)
  g <- bmin * sqrt((log(bmax) - log(bmin)) / (bmax - bmin)) * exp(R*log(bmax / bmin))
  return(g)
}
# calculate_q <- function(x,y, gamma, k0, N) {
#   gx <- calculate_g(x = x, gamma = gamma, k0 = k0, N = N)
#   gy <- calculate_g(x = y, gamma = gamma, k0 = k0, N = N)
#   f <- gx * gy
# }



sim_static_sn <- function(gamma, k0,  N=1e4, plot=T, kmax = N-1, mu = 1, phi = 1e4) {
  l  <- rexp(n = N, rate = mu) # draw lambdas from Exp(mu) dist
  
  kmax_opt <- optimise(f = function(z) {
    abs(calculate_g(x = Inf, gamma = gamma, k0 = k0, N = N, kmax = z)$g - 1) 
    }, 
    interval = c(1, N))
  kmax_opt <- floor(kmax_opt$minimum)
  kmax <- min(kmax, kmax_opt)
  
  gs <- calculate_g(x = l, gamma = gamma, k0 = k0, N = N, kmax = kmax, mu = mu) # calculate g(lamda)
  gl <- gs$g
  const <- gs$c
  Q <- outer(X = gl, Y = gl, FUN = '*')
  min_q <- 1 - exp(-1/50/phi) / (1 + 1/phi)
  Q <- pmax(Q, min_q)
  U <- matrix(1, nrow=N, ncol=N)
  U[row(U) < col(U)] <- runif(N * (N - 1) / 2)
  M <- matrix(0, nrow = N, ncol = N)
  M[U < Q] <- 1
  g <- graph_from_adjacency_matrix(M, mode = "undirected")
  dd <- degree_distribution(g)
  prop0 <- dd[1]
  dd <- dd[-1]/sum(dd[-1])
  const <- (1-gamma) / ((kmax * (1 - prop0)) ^ (1-gamma) - k0 ^ (1-gamma))
  
  inputs <- list(lambdas = l, gamma=gamma, k0 = k0, kmax = kmax, mu = mu)
  sn <- list(inputs = inputs, res=M, g=g, dd=dd, prop0 = prop0, const=const)
  if(plot) plot_dd(sn)
  
  return(sn)
}
