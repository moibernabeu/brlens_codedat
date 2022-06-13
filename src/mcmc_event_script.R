# MCMC sampling for data
# Moisès Bernabeu
# Barcelona, May 2022

# Loading sampling function
source('mcmc_event_sampling_function.R')
source('mcmc_analysis_funs.R')

require(ggplot2)
require(ggpubr)
require(coda)

theme_set(theme_bw())

# Loading data
load('../test/event_dist.RData')

events <- c('vert', 'met')

for (event in events) {
  # Executing sampling
  system.time(
    mcmcout <- mcmcfun(y = evdat[, paste(event, 'ndist', sep = '_')],
                       nchains = 3, niter = 100000,
                       thin = 3, unifmax = 100, oname = event)
  )
  
  # Plotting posterior distributions
  pdf(sprintf('../outputs/jags_%s_plots.pdf', event),
      width = 6.3, height = 7)
  plot(mcmcout)
  dev.off()
  
  post_list <- mcmc.list(mcmcout)
  post_df <- as.data.frame(ggs(post_list))
  post_df$Chain <- as.factor(post_df$Chain)
  
  hists <- ggplot(post_df, aes(x = value, colour = Chain)) +
    geom_density() +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = paste('yedat', event))
  
  traces <- ggplot(post_df, aes(x = Iteration, y = value, colour = Chain)) +
    geom_line(size = 0.2, aes(lty = Chain)) +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = paste('yedat', event))
  
  pdf(sprintf('../outputs/jags_%s_plots_gg.pdf', event),
      width = 5.27, height = 8.39)
  ggarrange(hists, traces, align = 'hv', common.legend = TRUE, legend = 'bottom')
  dev.off()
  
  diagdf <- as.data.frame(ggs_diagnostics(post_df))
  summ <- summary(mcmcout)
  summdf <- cbind(summ$statistics[, 1:2],
                  'R' = diagdf[which(diagdf$Diagnostic == 'Rhat'), 4],
                  'Effective size' = effectiveSize(mcmcout))
  write.csv(summdf, file = sprintf('../outputs/jags_%s_R_Effsz.csv', event))
  
  post_ac <- get_ac_df(post_df, 100)
  pdf(sprintf('../outputs/jags_%s_ac.pdf', event), width = 7, height = 7)
  ggplot(post_ac, aes(x = Lag, y = Autocorrelation, colour = Chain)) +
    geom_line() +
    facet_grid(Parameter ~ Chain) +
    labs(title = event)
  dev.off()
}


load('../outputs/jags_vert_mcmc.RData')
