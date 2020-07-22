devtools::load_all()

plots_file <- 'Fig5'

data("unvotes_original")
p <- plot_V(unvotes_original, xlab="vote", ylab = "country", aspect.ratio=1/3)
dataset <- "unvotes100_bw"
ggsave(p, filename = paste0(plots_file, '_', dataset, '.eps'), 
       height=5, width=12, units='cm')

data("lastfm")
p <- plot_V(lastfm, xlab="band", ylab = "user", aspect.ratio=4.2/1)
dataset <- "lastfm"
ggsave(p, filename = paste0(plots_file, '_', dataset, '.eps'), 
       height=20, width=6, units='cm')

data("paleo")
p <- plot_V(paleo, xlab="location", ylab = "genus", aspect.ratio=0.4)
dataset <- "paleo"
ggsave(p, filename = paste0(plots_file, '_', dataset, '.eps'), 
       height=5, width=10, units='cm')

data("catalanparliament")
p <- plot_V(catalanparliament, xlab="MP", ylab = "MP", aspect.ratio=1)
dataset <- "parliament"
ggsave(p, filename = paste0(plots_file, '_', dataset, '.eps'), 
       height=4.2, width=4.2, units='cm')

data("animals")
p <- plot_V(t(animals), xlab="animal", ylab = "attribute", aspect.ratio=1.15)
dataset <- "animals"
ggsave(p, filename = paste0(plots_file, '_', dataset, '.eps'), 
       height=3.93, width=3.5, units='cm')
