library(dplyr)
library(Seurat)

#####Load Data and DownSample
hum_met <- readRDS('Documents/davepOrtalData/AllenMTG/human_meta.RDS')
chi_met <- readRDS('Documents/davepOrtalData/AllenMTG/chimp_meta.RDS')
gor_met <- readRDS('Documents/davepOrtalData/AllenMTG/gorilla_meta.RDS')

n <- 50000
set.seed(42) #answer to all things

hum_mat <- readRDS('Documents/davepOrtalData/AllenMTG/human_mat.RDS')
hum_samps <- hum_mat %>% colnames %>% sample(size=n)
hum_mat <- hum_mat[,hum_samps]
hum_met <- hum_met %>% filter(sample_id %in% hum_samps)
hum_samps <- NULL

chi_mat <- readRDS('Documents/davepOrtalData/AllenMTG/chimp_mat.RDS')
chi_samps <- chi_mat %>% colnames %>% sample(size=n)
chi_mat <- chi_mat[,chi_samps]
chi_met <- chi_met %>% filter(sample_id %in% chi_samps)
chi_samps <- NULL

gor_mat <- readRDS('Documents/davepOrtalData/AllenMTG/gorilla_mat.RDS')
gor_samps <- gor_mat %>% colnames %>% sample(size=n)
gor_mat <- gor_mat[,gor_samps]
gor_met <- gor_met %>% filter(sample_id %in% gor_samps)
gor_samps <- NULL


#####Process gene_ids
orthos <- readRDS('Documents/davepOrtalData/AllenMTG/orthologous_genes.RDS') %>%
  mutate_all(as.character)

#Set to only ortho genes
hum_mat <- hum_mat[rownames(hum_mat) %in% orthos$human_symbol,]
chi_mat <- chi_mat[rownames(chi_mat) %in% orthos$chimp_symbol,]
gor_mat <- gor_mat[rownames(gor_mat) %in% orthos$gorilla_symbol,]

##Rename apes to human ids
#first enforce ortho order
hum_mat <- hum_mat[orthos$human_symbol,]
chi_mat <- chi_mat[orthos$chimp_symbol,]
gor_mat <- gor_mat[orthos$gorilla_symbol,]

#then rename
rownames(chi_mat) <- rownames(hum_mat)
rownames(gor_mat) <- rownames(hum_mat)

#Put the metadata cols in the same order
chi_met <- chi_met %>% select(one_of(hum_met %>% colnames))
gor_met <- gor_met %>% select(one_of(hum_met %>% colnames))


#####Join Data
comb_mat <- cbind(hum_mat,chi_mat) %>% cbind(gor_mat)
comb_met <- rbind(hum_met,chi_met) %>% rbind(gor_met) %>% as.data.frame
rownames(comb_met) <- comb_met$sample_id

#Make seurat
seu <- CreateSeuratObject(counts = comb_mat,
                          meta.data = comb_met,
                          assay = 'RNA',
                          project = "AllenMTG_v0.3")

saveRDS(seu,'Documents/davepOrtalData/AllenMTG/seu_v3_gor.RDS')

seu <- readRDS('Documents/davepOrtalData/AllenMTG/seu_v3_gor.RDS')

seu@meta.data$reg<- 'MTG'

seu@meta.data$ct<-''
seu@meta.data$ct[seu@meta.data$class == 'inh'] <- 'GABA'
seu@meta.data$ct[seu@meta.data$class == 'exc'] <- 'Glut'
seu@meta.data$ct[seu@meta.data$subclass == 'Astro'] <- 'Astr'
seu@meta.data$ct[seu@meta.data$subclass == 'Oligo'] <- 'Olig'
seu@meta.data$ct[seu@meta.data$subclass == 'OPC'] <- 'Olig'

#blank out for making bio_ct
seu@meta.data$st<-seu@meta.data$subclass
seu@meta.data$st[seu@meta.data$subclass == 'Astro'] <- ''
seu@meta.data$st[seu@meta.data$subclass == 'Oligo'] <- ''

#lets improve the inh subclass lables
seu@meta.data$st[seu@meta.data$subclass == 'Vip'] <- 'cge.Vip'
seu@meta.data$st[seu@meta.data$subclass == 'Sst Chodl'] <- 'mge.Sst-Chodl'
seu@meta.data$st[seu@meta.data$subclass == 'Sst'] <- 'mge.Sst'
seu@meta.data$st[seu@meta.data$subclass == 'Sncg'] <- 'cge.Sncg'
seu@meta.data$st[seu@meta.data$subclass == 'Pvalb'] <- 'mge.Pvalb'
seu@meta.data$st[seu@meta.data$subclass == 'Pax6'] <- 'cge.Pax6'
seu@meta.data$st[seu@meta.data$subclass == 'Lamp5_Lhx6'] <- 'llc.Lamp5-Lhx6'
seu@meta.data$st[seu@meta.data$subclass == 'Lamp5'] <- 'llc.Lamp5'
seu@meta.data$st[seu@meta.data$subclass == 'Chandelier'] <- 'llc.Chandelier'

seu@meta.data$subclass[seu@meta.data$subclass == 'Vip'] <- 'cge.Vip'
seu@meta.data$subclass[seu@meta.data$subclass == 'Sst Chodl'] <- 'mge.Sst-Chodl'
seu@meta.data$subclass[seu@meta.data$subclass == 'Sst'] <- 'mge.Sst'
seu@meta.data$subclass[seu@meta.data$subclass == 'Sncg'] <- 'cge.Sncg'
seu@meta.data$subclass[seu@meta.data$subclass == 'Pvalb'] <- 'mge.Pvalb'
seu@meta.data$subclass[seu@meta.data$subclass == 'Pax6'] <- 'cge.Pax6'
seu@meta.data$subclass[seu@meta.data$subclass == 'Lamp5_Lhx6'] <- 'llc.Lamp5-Lhx6'
seu@meta.data$subclass[seu@meta.data$subclass == 'Lamp5'] <- 'llc.Lamp5'
seu@meta.data$subclass[seu@meta.data$subclass == 'Chandelier'] <- 'llc.Chandelier'

#make bio_ct
seu <- seu %>% subset(ct != '')
seu@meta.data <- seu@meta.data %>% mutate(bio_ct =
                                            substring(species,1,3) %>% paste(ct,sep='') %>% paste(reg,sep='_') %>%
                                            paste(st,sep='.')) %>%
  mutate(bio_ct = bio_ct %>% stringr::str_replace_all(' ','-') %>% stringr::str_remove("\\.$")) %>%
  select(bio_ct,ct,subclass,reg,species,orig.ident)
seu@meta.data <- seu@meta.data %>%
  mutate(bio_ct = bio_ct %>% factor(levels = seu@meta.data$bio_ct %>% unique %>% sort))

saveRDS(seu,'Documents/davepOrtalData/AllenMTG/seu_comb_50k.rds')
