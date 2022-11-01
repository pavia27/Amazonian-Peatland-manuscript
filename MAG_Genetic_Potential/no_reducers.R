rm(list = ls(all.names = TRUE))
setwd("")
library(treeio)
library(tidyverse)
library(ggtree)
library(egg)
library(patchwork)
a<-read.newick("NO_reducer_GToTree.tre")
tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}
#tip and node labels
#ggtree(a,size=0.3,right=T)+geom_tiplab(size=1,align=F,linesize =0.15)+geom_text(aes(label=node),size=2)
#add supplementary information
b <- read.csv("heatmap_to_tree.csv")
c<-ggtree(a,size=0.3,right=TRUE)+geom_tiplab(size=0,align=F,linesize =0.15)
d<-c$data
d <- d[!d$isTip,]
d$label = as.numeric(d$label)
d[is.na(d)] <- 0
d$label <- round(d$label * 100)
e <- c + geom_point2(data=d,aes(fill=cut(label,c(0,70,100))),shape=21,size=1) +
  scale_fill_manual(values=c("black","white"),guide='none',breaks=c('(70,100]','(0,70]'))
f <- e %<+% b + geom_tippoint(aes(shape=Site,color=as.factor(copynumber)),size=3) + 
  scale_shape_manual(values=c(15,16,17)) + scale_color_grey()
g<-b[,c(1,9:25)] %>% reshape2::melt()
h <- ggplot(g, aes( tree_y(f, g), variable)) + geom_tile(aes(fill=value),color = "black") + 
  theme(legend.position="none") + scale_fill_gradient(low = "white", high = "black")+
  coord_flip() +theme_article() + scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0))
i<-ggtree(a,size=0.3,right=TRUE, branch.length='none')+geom_tiplab(size=0,align=F,linesize =0.15)
j <- i %<+% b + geom_tippoint(aes(shape=nortype),size=3) + 
  scale_shape_manual(values=c(10,1,19))
f + j + h + plot_annotation(tag_levels="A")

#phylogeny
f <- e %<+% b + geom_tippoint(aes(,color=class),size=3) +scale_color_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000',"red","blue"))
