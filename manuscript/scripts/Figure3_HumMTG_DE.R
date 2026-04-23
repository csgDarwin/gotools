## de time
#setup
library(Seurat)
library(dplyr)

seu <- readRDS('Documents/davepOrtalData/AllenMTG/seu_comb_50k.rds')

seu <- seu %>% NormalizeData()

#lists of idents
idents.1 <- seu@meta.data %>% filter(species == 'human')
idents.1 <- idents.1$bio_ct %>% as.character %>% unique
idents.2 <- idents.1 %>% stringr::str_replace_all('hum','chi')
idents.3 <- idents.1 %>% stringr::str_replace_all('hum','gor')

#do testing
Idents(seu) <- 'bio_ct'

chi_de <- res <- purrr::map2(idents.1,idents.2,
                             .f = function(.x,.y){
                               res <- seu %>% FindMarkers(.x,.y) %>% tibble::rownames_to_column('gene_id')
                               res$ct <- .x
                               return(res)
                             })
gor_de <- res <- purrr::map2(idents.1,idents.3,
                             .f = function(.x,.y){
                               res <- seu %>% FindMarkers(.x,.y)  %>% tibble::rownames_to_column('gene_id')
                               res$ct <- .x
                               return(res)
                             })

#wrangle
chi_de2 <- chi_de %>% purrr::reduce(full_join)
gor_de2 <- gor_de %>% purrr::reduce(full_join)

chi_de3 <- chi_de2 %>% filter(p_val_adj < 0.05) %>%
  mutate(chi_dir = as.integer(avg_log2FC > 0)*2-1) %>% select(gene_id,chi_dir,ct)

gor_de3 <- gor_de2 %>% filter(p_val_adj <  0.05) %>%
  mutate(gor_dir = as.integer(avg_log2FC > 0)*2-1) %>% select(gene_id,gor_dir,ct)


de <- full_join(chi_de3,gor_de3)
de[is.na(de)] <- 0

de2 <- de %>% filter(chi_dir == gor_dir) %>% mutate(dir = chi_dir) %>% select(gene_id,dir,ct)

de_cnts <- de2 %>% group_by(ct) %>% summarise(cnt = n())

de_cnts %>% ggplot(aes(y=ct,x=cnt)) + geom_bar(stat='identity') + ylab('') + xlab('# Human DEGs')

write.csv(de2,'Documents/davepOrtalData/AllenMTG/human_de_50k.csv',row.names = FALSE)

###isssssha boy still got it. Lets cluster it and find "the wall"
#make matrix
de_mat <- de2 %>% tidyr::pivot_wider(names_from = ct,values_from = dir)

de_mat2 <- de_mat %>% select(!one_of('gene_id')) %>% as.matrix
de_mat2[is.na(de_mat2)] <- 0

colnames(de_mat2) <- colnames(de_mat)[2:length(colnames(de_mat))]
rownames(de_mat2) <- de_mat$gene_id

##hclust genes
hc <- de_mat2 %>%
  scale() %>% dist() %>% hclust()
gene_order <- hc$labels[hc$order]

#id clusters
gene_cut1000 <- cutree(hc, k = 1000) %>%
  as.data.frame %>% tibble::rownames_to_column('gene_id')
colnames(gene_cut1000) <- c('gene_id','cluster')

keeps <- gene_cut1000 %>% group_by(cluster) %>% summarize(cnt=n()) %>%
  filter(cnt >= 100)

gene_cut1000 <- gene_cut1000 %>% filter(cluster %in% keeps$cluster)

##hclust cts
hc<- de_mat2 %>% t() %>%
  scale() %>% dist() %>% hclust()
ct_order <- hc$labels[hc$order]


#now force the orders using the factor levels and plot
de3 <- de2 %>% mutate(gene_id = gene_id %>% factor(levels = gene_order)) %>%
  mutate(ct = ct %>% factor(levels=ct_order))

a <- de3 %>% filter(!is.na(ct)) %>% ggplot(aes(x=ct,y=gene_id,fill=dir)) +
  geom_tile() +
  scale_y_discrete('DEGs', drop=FALSE) +
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.position = 'none')

gene_cut1000 <- gene_cut1000 %>%
  mutate(gene_id = gene_id %>% factor(levels=gene_order))

de4 <- de2 %>% left_join(gene_cut1000)
write.csv(de4,'Documents/davepOrtalData/AllenMTG/human_de_50k_clust.csv',row.names = FALSE)



b <- gene_cut1000 %>% ggplot(aes(y=gene_id,fill=cluster %>% as.character,x='')) +
  geom_tile() + scale_y_discrete('DE Pattern Clusters', drop=FALSE) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = 'left')

egg::ggarrange(b,a,nrow = 1,widths = c(1,20))


ccres <- read.table('Documents/davepOrtalData/AllenMTG/ccre.9num1kb.lss.header.tsv',sep='\t',header=TRUE) %>%
  mutate(perc_fhss = Human / abs(start - end)) %>% mutate(gene_id = Gene_IDs)

gene_cut_ccres <- gene_cut %>% left_join(ccres,by= 'gene_id')



#write.csv(gene_cut_ccres,'Documents/davepOrtalData/AllenMTG/de_modules_ccre_merged.csv',)

#make some tables
clust <- 314

a <- gene_cut_ccres %>% filter(cluster == clust) %>%
  filter(!grepl('PLS',cCRE_class,fixed = TRUE)) %>% arrange(-perc_fhss) %>%
  select(gene_id,cCRE_ID,cCRE_class,perc_fhss,chr,start,end) %>% head(5)

b <- gene_cut_ccres %>% filter(cluster == clust) %>%
  filter(grepl('PLS',cCRE_class,fixed = TRUE)) %>% arrange(-perc_fhss) %>%
  select(gene_id,cCRE_ID,cCRE_class,perc_fhss,chr,start,end) %>% head(5)

rbind(a,b)


#a quick hyper geometric????

#lets use 50% fhss as a cutoff for enhancers
cut_off <- .95

enh_white_balls <- dim(ccres %>% filter(perc_fhss > cut_off) %>% filter(gene_id != '') %>% filter(cCRE_class == 'dELS'))[1]
enh_total_balls <- dim(ccres %>% filter(gene_id != '')%>% filter(cCRE_class == 'dELS'))[1]

mods <- gene_cut_ccres$cluster %>% unique
names(mods) <- paste('mod',mods,sep='_')

enh_enrich_res <- list()

enh_res <- mods %>% lapply(function(mod){
  tmp <- gene_cut_ccres %>% filter(cluster == mod)
  enh_pulled_white <- dim(tmp %>% filter(perc_fhss > cut_off) %>% filter(cCRE_class == 'dELS'))[1]
  enh_total_pulls <- dim(tmp %>% filter(cCRE_class == 'dELS'))[1]
  res <- phyper(enh_pulled_white,
                enh_white_balls,
                enh_total_balls-enh_white_balls,
                enh_total_pulls,lower.tail = FALSE)
  col <- data.frame(fe = (enh_pulled_white/enh_total_pulls)/(enh_white_balls/enh_total_balls),
                    p=res)
  return(col)
})

tmp <- enh_res %>% purrr::reduce(rbind)
tmp$cluster <- names(enh_res)

enh_res <- tmp %>% mutate(p_adj <- p %>% p.adjust(method = 'bonf'))

enh_res
