calc_dd <- function(degree_vec, N) {
  
  dd <- rep(0, max(degree_vec))
  
  degree_table <- table(degree_vec)
  
  prop0 <- degree_table[1]/N
  
  degree_table <- degree_table[-1]
  dd[as.numeric(names(degree_table))] <- degree_table

  dd <- dd/sum(dd)
  
  return(list(prop0 = prop0, dd=dd))
}