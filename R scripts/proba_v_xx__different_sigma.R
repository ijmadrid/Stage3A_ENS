res <- read.csv("./Stage 3A ENS/DNA_Religation/projects/RCLPolymer_ExcludedVolumeExperiments/results/2018_05_31_17_39proba_VE_vs_Nc_NcinDF_adptEPSILONandSIGMA.csv")


library(ggplot2)

res$v = as.factor(round(res$excludedVolumeCutOff/res$encounterDistance))


res$v = as.factor(res$v)

### PROBA VS Nc
{
  p <- ggplot(data = res, aes(x = Nc, y = repair_probability, color = v)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = Nc, 
                                                                                                                                                   ymin = repair_probability - repair_CIhalflength,
                                                                                                                                                   ymax = repair_probability + repair_CIhalflength), alpha = 0.2) + theme_bw()
  p
}

### MFET VS Nc
{
  p <- ggplot(data = res, aes(x = Nc, y = meanFET, color = v)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = Nc, 
                                                                                                                                        ymin = meanFET - halfCI_FET,
                                                                                                                                        ymax = meanFET + halfCI_FET), alpha = 0.2) + theme_bw()
  p
}

### Repair MFET VS Nc
{
  p <- ggplot(data = res, aes(x = Nc, y = repairMFET, color = v)) + geom_line(size = 1.5) + theme_bw()
  p
}

### MSRG MFET VS Nc
{
  p <- ggplot(data = res, aes(x = Nc, y = Ensemble_MSRG, color = v)) + geom_line(size = 1.5) + geom_ribbon(data = res, aes(x = Nc, 
                                                                                                                                                   ymin = Ensemble_MSRG - Ensemble_MSRG_dx,
                                                                                                                                                   ymax = Ensemble_MSRG + Ensemble_MSRG_dx), alpha = 0.2) + theme_bw()
  p
}
