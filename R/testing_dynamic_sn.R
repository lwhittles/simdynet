rm(list=ls())
source("dynamic_sn.R")
source("calc_dd.R")
set.seed(2)
test_og <- sim_dynamic_sn(N = 1e4, gamma = 1.8, k0 = 0.5, phi = 1e4,
                          beta = 1, n_strain = 1, rho = 1, n_infs = 1,
                          t = 1+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365)
source("dynamic_sn.R")
set.seed(4)
test_new <- sim_dynamic_sn(N = 1e4, gamma = 1.8, k0 = 0.5, phi = 1e4,
                          beta = 1, rho = 1, n_infs = 1,
                          t = 1+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365)
test_new$log_infs
test_new$prop0_curr==test_og$prop0_curr
all(test_og$degree_curr == test_new$degree_curr)



debug(sim_dynamic_sn)
