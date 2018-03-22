rm(list=ls())
source("sn/src/dynamic_sn.R")
source("sn/src/calc_dd.R")

# test <- sim_dynamic_sn(N = 1e4, gamma = 1.6, k0 = 0.5, phi = 1e4, 
#                        beta = 1, n_strain = 1, rho = 1, n_infs = 1,
#                        t = 3+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
# 



# phis <- c(2e4, 3e4, 5e4)
# tests <- lapply(X = phis, FUN = function(x) {
#   sim_dynamic_sn(N = 1e4, gamma = 1.6, k0 = 0.5, phi = x, 
#                 beta = 1, n_strain = 1, rho = 1, n_infs = 0,
#                 t = 1+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
# })
# par(mfrow = c(2,2), mar=c(3,3,1,1), mgp =c(1.5,0.5,0), bty="n")
# plot_dd(test)
# lapply(X = tests, FUN = plot_dd)




plot_dd <-  function(sn) {
  with(sn, {
  n <- length(dd)
  xx <- 1:n
  with(inputs, {
    plot(x= xx, y= dd, pch=20, log="xy", col = col_line,
         xlim = c(1, 5e2), ylim = c(1e-4, 5e-1), 
         xlab = "number of partners over past year", ylab = "density",
         main= substitute(
           gamma~"="~g*";"~
             k[0]~"="~k0*";"~
             # k[max]~"="~kmax*";"~
             "p(k=0)="~pk0,
           list(g = gamma, 
                k0 = k0, 
                # kmax = kmax, 
                pk0 = round(prop0, 2))
         ))
    yy <- const * xx ^ (-gamma)
    lines(x=xx, y= yy, col=spectral_colours[3], lwd=2)
  })
})
}


gammas <- c(1.6, 1.7, 1.8, 1.9)
k0s <- c(0.4, 0.5, 0.6)


tests <- mapply(k0 = rep(k0s, 4), gamma = rep(gammas, each=3), FUN = function(k0, gamma) {
  sim_dynamic_sn(N = 1e4, gamma =gamma, k0 = k0, phi = 1e4,
                 beta = 1, n_strain = 1, rho = 1, n_infs = 0,
                 t = 1+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
}, SIMPLIFY = F)

pdf("sn/figs/dynamic_network_example.pdf", width = 10, height = 2.5)
par(mfcol = c(1,3), mar=c(3,3,2,0.5), mgp =c(1.7,0.5,0), bty="n")
lapply(X = tests[c(4,8,12)], FUN = plot_dd)
dev.off()
system('open "sn/figs/dynamic_network_example.pdf"')


lapply(X = tests[c(7,8,9)], FUN = plot_dd)
par(mfcol = c(1,1), mar=c(3,3,1,1), mgp =c(1.5,0.5,0), bty="n")
plot_dd(tests[[8]])

# 
# 
# phis <- rep(c(1e4, 1.25e4, 1.5e4), 3)
# k0s <- rep(c(0.3,0.4, 0.5), each = 3)
# 
# tests3 <- mapply(k0 = k0s, phi = phis, SIMPLIFY = F,FUN= function(k0, phi) {
#   sim_dynamic_sn(N = 1e4, gamma = 1.6, k0 = k0, phi = phi, 
#                  beta = 1, n_strain = 1, rho = 1, n_infs = 0,
#                  t = 1+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
# })
# lapply(X = tests3, FUN = plot_dd)

# load("sn/results/dynamic_network_tests.RData")


# 
# #### TESTING INFS
# 
# test_infs <- sim_dynamic_sn(N = 1e4, gamma = 1.6, k0 = 0.5, phi = 1e4, 
#                           n_strain = 1,n_infs = 1,
#                           beta = 70, psi = 0.95, sigma = 77, 
#                           nu = 2.3, eta = 2.3, mu = 150, rho = 53,
#                           t = 3+1/365, max.iter = 1e6, record = F, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
# 
# undebug(sim_dynamic_sn)
#  
# test_infs$ext_day
# 
# sort(sapply(test_infs, object.size))
# 
# head(test_infs$log_infs)
# 
# plot(test_infs$dd, log="xy")
# 
# 
# mean_ci(od)
# subset(test_infs$log_infs, p_full %in% infectors & c2=="S")
# subset(test_infs$log_infs, p==7175 & c1=="S", c("time", "p_ninf", "infector_full"))
# 
# od
# test_infs$log_infs
# test_infs$total_treated
# test_infs$ext
# which.max(test_infs$inputs$lambdas)
# test_infs$inputs$lambdas[test_infs$log_infs$p]
# order(test_infs$inputs$lambdas, decreasing = T)[1:10]
# 
# dim(unique(test$rec_rels))/dim(test$rec_rels) # around 94% of relationships are with new partners
# 
# pdf("sn/figs/example1_infs_plot.pdf", width = 12, height = 5)
# par(mfrow=c(2,2), mar=c(3,3,1,1), mgp =c(1.5,0.5,0), bty="n")
# plot_all(test_infs)
# cols <- spectral_colours[c(5,3, 11,10)]
# with(test_infs, matplot(x= times[-1], y= t(infs.t[1,,]),type = "l", lty=1, col = cols, 
#                         xlab = "time (years)", ylab = "number in infection stage"))
# legend("topleft", legend = c("U", "E", "A", "T"), fill = cols , bty="n")
# dev.off()
# system('open "sn/figs/example1_infs_plot.pdf"')
# 
# 
# test_infs$log_infs[!is.na(test_infs$log_infs$infector_ninf),]
# test_infs$log_infs$infector_stage <- factor(x = log_infs$infector_stage, levels = 2:4, labels = c( "U", "E", "A"), exclude = NA)
# 
# test_infs$log_infs$r0_id <- with(test_infs$log_infs, paste0(infector, "_", infector_ninf))
# test_infs$log_infs$r0_id[test_infs$log_infs$r0_id=="NA_NA"] <- NA
# hist(table(test_infs$log_infs$r0_id), breaks = 0.5 + 0:10)
# mean(table(test_infs$log_infs$r0_id))
# hist(table(unique_infs))
# 
# par(mfrow=c(1,1))
# table(test_infs$log_infs$infector)[-1]
# 
# 
# 
# test_infs2 <- sim_dynamic_sn(N = 1e4, gamma = 1.6, k0 = 0.4, phi = 1e4, 
#                         n_strain = 2, n_infs = c(10, 10), prop_pres = c(0.75,0.25), efficacy = matrix(data = c(1,1,
#                                                                                                                1,0.17), nrow = 2, byrow = T),
#                         beta = c(10,10), psi = 0.6, sigma = 77, alpha = c(1, 1.8),
#                         nu = c(1.8, 1.8), eta = 2, mu = 136, rho = 53,
#                         t = 5+1/365, max.iter = 1e6, record = T, record_term = 1, record_lengths = T, burn.in = 1/365, opt_phi = F)
# 
# # tapply(X = pmax(diff(test_infs2$infs.t[1,"U",]), 0), INDEX = ceiling(test2$times[-(1:2)]), FUN = sum)
# # tapply(X = pmax(diff(test2$infs.t[2,"U",]), 0), INDEX = ceiling(test2$times[-(1:2)]), FUN = sum)
# 
# 
# pdf("sn/figs/example_strains_plot.pdf", width = 12, height = 6)
# par(mfrow=c(2,3), bty="n", mar = c(3,3,1,1), mgp = c(1.5,0.5, 0))
# cols <- spectral_colours[c(5,3, 11,10)]
# with(test_infs2, matplot(x= times[-1], y= t(infs.t[1,,]),type = "l", lty=1, col = cols, ylim = c(0, 200), xlab = "time (years)", ylab = "strain 1 prevalence"))
# legend("topleft", legend = c("U", "E", "A", "T"), fill = cols , bty="n")
# with(test_infs2, matplot(x= times[-1], y= t(infs.t[2,,]),type = "l", lty=1, col = cols, ylim = c(0,200), xlab = "time (years)", ylab = "strain 2 prevalence"))
# legend("topleft", legend = c("U", "E", "A", "T"), fill = cols , bty="n")
# cols <- spectral_colours[c(1, 9)]
# with(test_infs2, matplot(x = times[-1], y = apply(X = infs.t, MARGIN = 1, FUN = colSums),type = "l", 
#                     lty=1, col = cols, ylim = c(0,200),
#                     xlab = "time (years)", ylab = "prevalence"))
# legend("topleft", legend = paste0("strain ", 1:2), fill = cols , bty="n")
# 
# plot_all(test_infs2)
# dev.off()
# system('open "sn/figs/example_strains_plot.pdf"')
# 
# save.image(file = "sn/results/dynamic_network_tests.RData")
