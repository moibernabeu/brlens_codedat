# Parsing seed to event data
# Moisès Bernabeu
# Barcelona, April 2022

evdat <- read.csv('../outputs/0076_dist.csv')
evdat <- evdat[-which(evdat$vert_dist == evdat$met_dist), ]
evdat$wnwidthratio <- evdat$whole_width / evdat$norm_width
evdat <- evdat[-which(evdat$wnwidthratio >=
                        quantile(evdat$wnwidthratio, 0.95)), ]
evdat <- evdat[-which(evdat$vert_dist == 0 | evdat$met_dist == 0), ]

save(evdat, file = '../event/event_dist.RData')
