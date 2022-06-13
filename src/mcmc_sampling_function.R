# MCMC sampling function using JAGS
# Moisès Bernabeu
# Barcelona, May 2022

require(rjags)

mcmcfun <- function(spto, dat, nchains = 3, niter = 1000,
                    thin = 3, unifmax = 100, oname) {
  # Model definition
  model_text <- sprintf('model{
    # Likelihood
    for (i in 1:n) {
      y[i] ~ dgamma(a, b)
    }
    
    m <- a / b
    v <- a / b^2
    mo <- ifelse(a <= 1, 0, (a - 1) / b)
    
    # Prior distributions
    a ~ dunif(0, %d)
    b ~ dunif(0, %d)
  }', unifmax, unifmax)
  
  # Data generation
  y <- dat[which(dat$sp_to == spto & dat$ndist_A != 0), 'ndist_A']
  N <- length(y)
  
  datlist <- list(y = y, n = N)
  
  # Preparing the model
  model_jags <- jags.model(textConnection(model_text), 
                           data = datlist,
                           n.chains = nchains,
                           n.adapt = round(niter * 0.1, digits = 0))
  
  # Running the mcmc sampling
  post <- coda.samples(model_jags, 
                       variable.names = c('a', 'b', 'm', 'v', 'mo'), 
                       n.iter = niter,
                       thin = thin)
  
  # Assigning species name to the posterior sample
  assign(sprintf('%s_mcmc', spto), post)
  
  # Saving the posterior sample and its information into RData file
  save(list = c('N', 'nchains', 'niter', 'thin', 'unifmax', sprintf('%s_mcmc', spto)),
       file = sprintf('../outputs/jags_%s_%s_mcmc.RData', oname, spto))
  
  return(post)
}
