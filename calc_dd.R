calc_dd <- function(degree_vec, N) {
  dd <- rep(0, max(degree_vec))
  degree_table <- table(degree_vec)
  dd[as.numeric(names(degree_table)[-1])] <- degree_table[-1]
  prop0 <- degree_table[1]/N
  dd <- dd/sum(dd)
  return(list(prop0 = prop0, dd=dd))
}