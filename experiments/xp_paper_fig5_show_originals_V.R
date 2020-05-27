devtools::load_all()

data("unvotes100")
p <- plot_V(unvotes100, xlab="vote", ylab = "country", aspect.ratio=1/3)
dataset <- "unvotes100"
ggsave(p, filename = paste0('fig_V_dmkd_', dataset, '.eps'), 
       height=5, width=12, units='cm')

data("unvotes100coldwar")
p <- plot_V(unvotes100coldwar, xlab="vote", ylab = "country", aspect.ratio=1/3)
dataset <- "unvotes100coldwar"
ggsave(p, filename = paste0('fig_V_dmkd_', dataset, '.eps'), 
       height=5, width=12, units='cm')

data("unvotes100_abs")
p <- plot_V(unvotes100_abs, xlab="vote", ylab = "country", aspect.ratio=1/3)
dataset <- "unvotes100_abs_chronological"
ggsave(p, filename = paste0('fig_V_dmkd_', dataset, '.eps'), 
       height=5, width=12, units='cm')

data("lastfm")
p <- plot_V(lastfm, xlab="band", ylab = "user", aspect.ratio=4.2/1)
dataset <- "lastfm"
ggsave(p, filename = paste0('fig_V_dmkd_', dataset, '.eps'), 
       height=20, width=6, units='cm')

data("paleo")
p <- plot_V(paleo, xlab="location", ylab = "genus", aspect.ratio=0.4)
dataset <- "paleo"
ggsave(p, filename = paste0('fig_V_dmkd_', dataset, '.eps'), 
       height=5, width=10, units='cm')

data("catalanparliament")
p <- plot_V(catalanparliament, xlab="MP", ylab = "MP", aspect.ratio=1)
dataset <- "parliament"
ggsave(p, filename = paste0('Fig5c_V_dmkd_', dataset, '.eps'), 
       height=4.2, width=4.2, units='cm')

data("animals")
p <- plot_V(t(animals), xlab="animal", ylab = "attribute", aspect.ratio=1.15)
dataset <- "animals"
ggsave(p, filename = paste0('Fig5b_V_dmkd_', dataset, '.eps'), 
       height=3.93, width=3.5, units='cm')
