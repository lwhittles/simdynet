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
