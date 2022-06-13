# Functions to analyse MCMC
# Moisès Bernabeu
# May 2022

# Some of these functions are adapted from ggmcmc and coda to better performance
# sacrificing some useless arguments that also made the functions to be slower.

require(ggmcmc)
require(tidyr)
require(dplyr)
require(cumstats)


get_ac_df <- function(mcmc_df, nLags = 50) {
  # Retrieve autocorrelation per parameter and chain
  wc.ac <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    do(ac(.$value, nLags))

  return(wc.ac)
}


get_full_means <- function(mcmc_df) {
  # Estimates of the mean by iteration
  dm.m <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    summarize(m = mean(value))
  
  return(dm.m)
}

get_running <- function(mcmc_df) {
  # Calculate the running mean
  # Force the object to be sorted by Parameter, and hence avoid 'rm' calculation
  # to be wrong
  dm.rm <- mcmc_df %>%
    arrange(Parameter, Iteration) %>%
    group_by(Parameter, Chain) %>%
    mutate(rm = cumsum(value) / Iteration,
           sd = sqrt(cumvar(value)),
           up = rm + sd,
           do = rm - sd)
  
  return(dm.rm)
}
