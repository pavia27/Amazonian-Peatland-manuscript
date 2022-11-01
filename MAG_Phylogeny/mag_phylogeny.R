rm(list = ls(all.names = TRUE))
setwd("/Users/michael_pavia/Dropbox (ASU)/03_Shared_Folders/MichaelPavia_coordinations_data/Data/00_Metagenome_Analysis/04_MAGs/01_gtdbtk/00_GToTree/")
library(treeio)
library(tidyverse)
library(ggtree)
library(egg)
library(Ckmeans.1d.dp)
library(scatterpie)
library(reshape2)
#bacteria
a<-read.newick("bacteria_GToTree_output.tre")
#tip and node labels
#ggtree(a,size=0.3,right=T)+geom_tiplab(size=1,align=F,linesize =0.15)+geom_text(aes(label=node),size=2)
b<-rootnode(a)
c<-ggtree(a,size=0.3,right=F)+geom_tiplab(size=0,align=F,linesize =0.15)
d <- c$data
d <- d[!d$isTip,]
d$label = as.numeric(d$label)
d[is.na(d)] <- 0
d$label <- round(d$label * 100)
e<- c + geom_point2(data=d,aes(fill=cut(label,c(0,70,100))),shape=21,size=1)+
  scale_fill_manual(values=c("black","white"),guide='none',breaks=c('(70,100]','(0,70]'))
e + geom_cladelab(node=559,label="Desulfobacterota_G",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=679,label="Desulfobacterota_B",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=784,label="Verrucomicrobiota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=401,label="Methylomirabilota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=691,label="Actinobacteriota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=742,label="Eremiobacterota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=409,label="Acidobacteriota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=579,label="Proteobacteria",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=482,label="Myxococcota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=747,label="Chloroflexota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=794,label="Patescibacteria",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=564,label="Myxococcota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=531,label="Desulfobacterota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=396,label="Dependentiae",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=385,label="Chlamydiota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=374,label="Eisenbacteria",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=357,label="Dormibacterota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=346,label="Cyanobacteria",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=772,label="Bacteroidota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=575,label="Myxococcota_A",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=496,label="Nitrospirota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=790,label="Spirochaetota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=781,label="Planctomycetota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=746,label="CSP1-3",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=765,label="Gemmatimonadota",align=T,offset.text=0, barsize=0.5, fontsize=1)

#archaea
f<-read.newick("archaea_GToTree.tre")
g<-root(f,outgroup ="SUPR01.1", edgelabel = TRUE)
#tip and node labels
#ggtree(g,size=0.3,right=T)+geom_tiplab(size=1,align=F,linesize =0.15)+geom_text(aes(label=node),size=2)
h <-ggtree(g,size=0.3,right=F)+geom_tiplab(size=0,align=F,linesize =0.15)
i <- h$data
i <- i[!i$isTip,]
i$label = as.numeric(i$label)
i[is.na(i)] <- 0
i$label <- round(i$label * 100)
j <- h + geom_point2(data=i,aes(fill=cut(label,c(0,70,100))),shape=21,size=1)+
  scale_fill_manual(values=c("black","white"),guide='none',breaks=c('(70,100]','(0,70]'))
j + geom_cladelab(node=148,label="Thermoproteota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=211,label="Thermoplasmatota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=188,label="Halobacteriota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=181,label="Asgardarchaeota",align=T,offset.text=0, barsize=0.5, fontsize=1)+
  geom_cladelab(node=184,label="Methanobacteriota",align=T,offset.text=0, barsize=0.5, fontsize=1)

#corresponding mag abundance
k<-read.csv("MAG_Metrics.csv") %>% 
  select("Bin.Id","Site","Phylum","Domain") %>%
  group_by(Site,Phylum) %>%
  summarize((count=n())) %>%
  `colnames<-`(c("Site","Phylum","Count")) %>% 
  mutate(Phylum=gsub("p__", "",Phylum)) %>%
  spread(Site, Count) %>% replace(is.na(.), 0) %>% 
  rowwise() %>% mutate(Total = sum(c_across(2:4)), holderx = 10) %>% 
  rowid_to_column(var="holdery") %>% as.data.frame()
m <- Ckmeans.1d.dp(k$Total,4)
k$cluster<-m[["cluster"]]
ggplot(k) + geom_scatterpie(aes(x=holderx, y=holdery*5, group=Phylum,r=cluster),color=NA, data=k,cols=c("Buena Vista","Quistococha","San Jorge" )) + 
  coord_equal() + geom_scatterpie_legend(k$cluster, x=-20, y=0) +
  scale_fill_manual(values=c("#252525","#969696","#f0f0f0"))+theme_article()+
  geom_text(aes(x=holderx, y=holdery*5,label = k$Phylum),nudge_x = -10,size = 2)

for (l in unique(k$cluster)) {
  m<-k[k$cluster==l,]
  print(paste0("cluster:",l))
  print(max(m$Total))
  print(min(m$Total))
}
