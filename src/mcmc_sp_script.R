# MCMC sampling for data
# Moisès Bernabeu
# Barcelona, May 2022

# To run it in every species, it has to be duplicated and change the index of
# the vector sps[1] to sps[n] where n is the species number. This is done to
# allow the paralelisation per species using greasy. Note the dataset has also
# to be changed from yedat to hudat for the human phylome. All this can be done
# with a bash script with the command sed.

# Loading sampling function
source('mcmc_sampling_function.R')
source('../scripts/mcmc_analysis_funs.R')

# Loading data
load('../data/seed2sp_dist.Rdata')

# Getting species list
sps <- yedat$sp_to[!duplicated(yedat$sp_to)]
length(sps)

# Executing sampling
system.time(
  mcmcout <- mcmcfun(sps[1], yedat, nchains = 3, niter = 100000,
                     thin = 3, unifmax = 100, oname = 'yedat')
)

# Plotting posterior distributions
pdf(sprintf('../outputs/jags_yedat_%s_plots.pdf', sps[1]),
    width = 6.3, height = 7)
plot(mcmcout)
dev.off()

post_list <- mcmc.list(mcmcout)
post_df <- as.data.frame(ggs(post_list))
post_df$Chain <- as.factor(post_df$Chain)

hists <- ggplot(post_df, aes(x = value, colour = Chain)) +
  geom_density() +
  facet_wrap(~Parameter, scales = 'free', ncol = 1) +
  labs(title = paste('yedat', sps[1]))

traces <- ggplot(post_df, aes(x = Iteration, y = value, colour = Chain)) +
  geom_line(size = 0.2, aes(lty = Chain)) +
  facet_wrap(~Parameter, scales = 'free', ncol = 1) +
  labs(title = paste('yedat', sps[1]))

pdf(sprintf('../outputs/jags_yedat_%s_plots_gg.pdf', sps[1]),
    width = 5.27, height = 8.39)
ggarrange(hists, traces, align = 'hv', common.legend = TRUE, legend = 'bottom')
dev.off()

diagdf <- as.data.frame(ggs_diagnostics(post_df))
summ <- summary(mcmcout)
summdf <- cbind(summ$statistics[, 1:2],
                'R' = diagdf[which(diagdf$Diagnostic == 'Rhat'), 4],
                'Effective size' = effectiveSize(mcmcout))
write.csv(summdf, file = sprintf('../outputs/jags_yedat_%s_R_Effsz.csv', sps[1]))

post_ac <- get_ac_df(post_df, 100)
pdf(sprintf('../outputs/jags_yedat_%s_ac.pdf', sps[1]), width = 7, height = 7)
ggplot(post_ac, aes(x = Lag, y = Autocorrelation, colour = Chain)) +
  geom_line() +
  facet_grid(Parameter ~ Chain) +
  labs(title = paste('yedat', sp[1]))
dev.off()
