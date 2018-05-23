res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_CLs@DamageFoci/results/2018_05_23_12_26_proba-v_gNc.csv")

library(ggplot2)
library(RColorBrewer)

res$excludedVolumeCutOff <- as.factor(res$excludedVolumeCutOff)

{
  p <- ggplot(data = res, aes(x = excludedVolumeSpringConstant, y = repair_probability, color = excludedVolumeCutOff)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = excludedVolumeSpringConstant, 
                                                                                                                                                                             ymin = repair_probability - repair_CIhalflength,
                                                                                                                                                                             ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
  
  colourCount = 15
  getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  
  p + geom_vline(xintercept = 0.6) + scale_color_manual(values = getPalette(colourCount)) + coord_cartesian()
}



####
{
  w = c(rep(rep(1,5),3),rep(9,5),rep(rep(4,5),10))
  p <- ggplot(data = res, aes(x = genomicDistance, y = Nc, width = w)) + 
   geom_tile(aes(fill=repair_probability))
  
  p +  scale_fill_gradientn(colours= heat.colors(100))
}
####

library(plotly)

#######
pp <- ggplotly(p)
ggplotly(p) %>% add_surface(z = pp$z)
######

########
proba <- matrix(res$repair_probability, ncol = 14)
mfet <- matrix(res$meanFET, ncol = 14)
gs <- c(seq(2,4),seq(9,49,4))
nc <- seq(2,19,4)
py <- plot_ly(x = gs, y = nc, z = proba) %>% add_surface()
py
#########

########
proba.REMOVING <- matrix(res$repair_probability, ncol = 14)
proba.KEEPING <- matrix(res$repair_probability, ncol = 14)
py <- plot_ly(x = gs, y = nc, z = proba.REMOVING) %>% add_surface()

py <- plot_ly(x = gs, y = nc) %>% add_surface(z = proba.KEEPING, cmin = min(proba.KEEPING), cmax = max(proba.KEEPING), colorscale = list(c(0,1),c("rgb(255,112,184)","rgb(128,0,64)"))) %>% 
  add_surface(z = proba.REMOVING,cmin = min(proba.REMOVING), cmax = max(proba.REMOVING), colorscale = list(c(0,1),c("rgb(107,184,214)","rgb(0,90,124)")))
py
#########