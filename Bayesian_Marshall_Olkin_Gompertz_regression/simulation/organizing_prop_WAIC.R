setwd("~/Code-for-papers/Bayesian_Marshall_Olkin_Gompertz_regression/simulation/simulation_results/comparing_models")



sim_WAIC_50 = read.csv2("1sim_WAIC_50.csv")[, -1]
sim_WAIC_100 = read.csv2("1sim_WAIC_100.csv")[, -1]
sim_WAIC_300 = read.csv2("1sim_WAIC_300.csv")[, -1]

sim_WAIC_500_1 = read.csv2("11sim_WAIC_500.csv")[, -1] #500
sim_WAIC_500_2 = read.csv2("21sim_WAIC_500.csv" )[, -1] #500

sim_WAIC_500 = rbind(sim_WAIC_500_1, sim_WAIC_500_2)

sim_WAIC_1000 = read.csv2("1sim_WAIC_1000.csv")[, -1]


suma_waic = function(table_waic){
  
  sum = table(factor(
    apply(table_waic, 1, function(x) {
      if (all(is.na(x))) return(NA_integer_)
      which.min(replace(x, is.na(x), Inf))
    }),
    levels = seq_len(ncol(table_waic))
  ))
  
  names(sum) = colnames(table_waic)
  
  return(sum)
}

suma_waic(sim_WAIC_50)/1000
suma_waic(sim_WAIC_100)/1000
suma_waic(sim_WAIC_300)/1000
suma_waic(sim_WAIC_500)/1000
suma_waic(sim_WAIC_1000)/1000



















































