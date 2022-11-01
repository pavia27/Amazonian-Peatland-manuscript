merge_gene_counts<-{
  rm(list = ls(all.names = TRUE))
  setwd("")
  library(tidyverse)
  b<-Sys.glob("*KO.csv")
  for (c in b) {
    d<-read.csv(c,quote = "",row.names = NULL,stringsAsFactors = FALSE)
    d<-d %>% group_by(KO) %>% summarize((count=n()))
    colnames(d)<- c("KO",substr(c,1,17))      a<-full_join(a,d, by = "KO")
    }
    a[is.na(a)] <- 0
    write.csv(a,"PMFB_Functional_Potential.csv")
  }
pco_stats<-{
  rm(list = ls(all.names = TRUE))
  setwd("")
  library(scatterplot3d)
  library(car)
  library(vegan)
  a<-read.csv("PMFB_Functional_Potential.csv")
  rownames(a)<-a[,2]
  a<-a[,-c(1:2)]
  sitegene <- as.factor(c("SJ","SJ","SJ","SJ","BV","BV","BV","BV","QUI","QUI","QUI","QUI","SJ","SJ","BV","BV","BV","BV","QUI","QUI","QUI","QUI","SJ","SJ"))
  pcoa.res <- capscale(t(a)~1, distance = "bray")
  summary(pcoa.res)
  b<-as.data.frame(scores(pcoa.res,choices = 1:3)$sites)
  b$site<-as.factor(c("SJO","SJO","SJO","SJO","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","SJO","SJO","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","SJO","SJO"))
  c = c(0, 1, 2) 
  c <- c[as.numeric(b$site)]
  scatter3d(b$MDS1,b$MDS2,b$MDS3, groups = b$site,grid = FALSE, surface = FALSE, ellipsoid = TRUE,radius=4)
  }
hclust<-{
  rm(list = ls(all.names = TRUE))
  setwd("")
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  a<-read.csv("PMFB_Functional_Potential.csv")
  rownames(a)<-a[,2]
  a<-a[,-c(1:2)]
  b<- as.data.frame(t(as.matrix(a)))
  hclustfunc <- function(x) hclust(x, method="ward.D")
  distfunc <- function(x) as.dist((1-cor(t(x)))/2)
  c <- distfunc(b)
  fit <- hclustfunc(c)
  dend <- as.dendrogram(fit)
  sitegene <- as.factor(c("SJ","SJ","SJ","SJ","BV","BV","BV","BV","QUI","QUI","QUI","QUI","SJ","SJ","BV","BV","BV","BV","QUI","QUI","QUI","QUI","SJ","SJ"))
  ggd1 <- as.ggdend(dend)
  ggplot(ggd1)
  ggplot(ggd1, horiz = TRUE, theme = NULL)+ scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
  dend1<-dend %>%set("leaves_pch",c(0,0,0,0,0,0,0,0,3,3,3,1,1,1,1,1,1,1,1,3,3,3,3,3))%>% set("leaves_cex",2)%>%plot()
  }
DESEQ<-{
    rm(list = ls(all.names = TRUE))
    library(DESeq2)
    library(tidyverse)
    library(genefilter)
    library(gplots)
    library(RColorBrewer)
    library(pheatmap)
    library(tidyverse)
    library(gridExtra)
    library(reshape2)
    library(vegan)
    setwd("")
    a<-read.csv("PMFB_Functional_Potential.csv")
    rownames(a)<-a[,2]
    a<-a[,-c(1:2)]
    col.order = data.frame(row.names = c("X0116_SJ02_MP02_10","X0116_SJ02_MP02_20","X0116_SJ02_MP15_10","X0116_SJ02_MP15_20","X0216_BV02_MP05_10","X0216_BV02_MP05_20","X0216_BV02_MP12_10","X0216_BV02_MP12_20","X0216_QUI2_MP05_10","X0216_QUI2_MP05_20","X0216_QUI2_MP10_10","X0216_QUI2_MP10_20","X0715_SJ02_MP02_20","X0715_SJ02_MP15_20","X0815_BVA2_MP04_20","X0815_BVA2_MP10_20","X1015_BV01_MP06_20","X1015_BV01_MP10_20","X1015_QUI2_SP05_10","X1015_QUI2_SP05_20","X1015_QUI2_SP10_10","X1015_QUI2_SP10_20","X1015_SJ02_MP02_20","X1015_SJ02_MP15_20"),treatment = c("SJO","SJO","SJO","SJO","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","SJO","SJO","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","SJO","SJO"))
    dds.ko = DESeqDataSetFromMatrix(countData = a, colData = col.order, design = ~treatment)
    dds<-DESeq(dds.ko)
    rld<-rlogTransformation(dds, blind=FALSE)
    #check data for heteroscedacisity
    dds_sizeF<-estimateSizeFactors(dds)
    cdsFull = estimateSizeFactors(dds)
    sizeFactors(dds_sizeF)
    normalized.func.data=counts(dds_sizeF, normalized = T)
    head(normalized.func.data)
    write.csv(as.data.frame(normalized.func.data),file="normlized-deseq-KO-data.csv")
    dds_disp = estimateDispersions(dds_sizeF)
    plotDispEsts(dds_disp)
    sampleDist<-dist(t(assay(rld)))
    dist.mat<-as.matrix(sampleDist)
    plot(log2(counts(dds_sizeF,normalize = T)[,1:8]+1), pch = 16, cex = 0.3 )
    dds.full = varianceStabilizingTransformation(cdsFull, blind = T)
    print(plotPCA(dds.full, intgroup= ("treatment")))
    res <- results(dds)
    head(res)
    res <- res[order(res$padj),]
    summary(res)
    plotMA(dds,main="DESeq2")
    write.csv(as.data.frame(res),file="KO_DESeq2-results.csv")
    hist(as.data.frame(res$pvalue))
    #To pull the LFC against the grand mean for generating the heatmaps
    TopVarGenes<-order(rowVars(assay(rld)), decreasing = TRUE)
    mat<- assay(rld)[TopVarGenes,]
    mat <- mat - rowMeans(mat)
    write.csv(as.data.frame(mat),file="KO_DESeq2-log2FC_against_grandMean.csv")
    #signficant genes were explored and chosen based on biogeochemcical significane 
    KO.heat<-read.csv("KO_DESeq2-log2FC_against_grandMean.csv", header = T)
    colnames(KO.heat)<-c("KO","SJO_0116_MP02_10cm","SJO_0116_MP02_20cm","SJO_0116_MP15_10cm","SJO_0116_MP15_20cm","BVA_0216_MP05_10cm","BVA_0216_MP05_20cm","BVA_0216_MP12_10cm","BVA_0216_MP12_20cm","QUI_0216_MP05_10cm","QUI_0216_MP05_20cm","QUI_0216_MP10_10cm","QUI_0216_MP10_20cm","SJO_0715_MP02_20cm","SJO_0715_MP15_20cm","BVA_0815_MP10_20cm","BVA_0815_MP04_20cm","BVA_1015_MP06_20cm","BVA_1015_MP10_20cm","QUI_1015_MP10_10cm","QUI_1015_MP10_20cm","QUI_1015_MP05_10cm","QUI_1015_MP05_20cm","SJO_1015_MP02_20cm","SJO_1015_MP15_20cm")
    #the "full" is what all genes are and "select" are only those in the heatmap
    #a<-read.csv("sig_genes_anno_full.csv")
    a<-read.csv("sig_genes_anno_select.csv")
    b<-semi_join(KO.heat,a, by = "KO") #keep only observations in df1 that match in df2.
    c<-full_join(a,b, by = "KO")
    d<-c[,-c(1,2,4:11)]
    d<-d%>% group_by(Complex) %>% summarise_all(list(mean))
    e<-unique(c[,c(3,5)])
    d<-as.data.frame(d)
    f<-full_join(d,e, by = "Complex")
    f<-as.data.frame(f)
    i<-as.data.frame(f)
    rownames(i) <- i[,1]
    i<-i[-1,-c(1,26)]
    col.order<- c("BVA_0216_MP05_10cm","BVA_0216_MP05_20cm","BVA_0216_MP12_10cm","BVA_0216_MP12_20cm","BVA_0815_MP10_20cm","BVA_0815_MP04_20cm","BVA_1015_MP06_20cm","BVA_1015_MP10_20cm","QUI_0216_MP05_10cm","QUI_0216_MP05_20cm","QUI_0216_MP10_10cm","QUI_0216_MP10_20cm","QUI_1015_MP10_10cm","QUI_1015_MP10_20cm","QUI_1015_MP05_10cm","QUI_1015_MP05_20cm","SJO_1015_MP02_20cm","SJO_1015_MP15_20cm","SJO_0116_MP02_10cm","SJO_0116_MP02_20cm","SJO_0116_MP15_10cm","SJO_0116_MP15_20cm","SJO_0715_MP02_20cm","SJO_0715_MP15_20cm")
    i<-i[,col.order]
    i$BVA<-rowMeans(i[,c(1:8)])
    i$SJO<-rowMeans(i[,c(9:16)])
    i$QUI<-rowMeans(i[,c(17:24)])
    i<-i[,-c(1:24)]
    i<-as.matrix(i)
    h<-data.frame(row.names = f$Complex, Metabolism = f$Label)
    #d<-read.csv("sig_genes_anno_select_ordered.csv")
    my_colour=list(Metabolism = c(Oligosaccharide_Degredation="#023858",Aromatic_Degredation="#0570b0",Saccharide_Degredation="#74a9cf",Pentose_phosphate_pathway="#d0d1e6",Methanogensis="#fff7fb",Nitrogen_Assimilation="#00441b",Denitrification="#5aae61",Sulfur_Assimilation="#fed976",Dissimilatory_Sulfur="#ffffcc",Phosphorous_Aquisition="#67000d",Iron_Acquisition="#cb181d",Inorganic_Ion_Transport="#fb6a4a",Biofilm_Formation="#d9d9d9",Low_pH_growth="#000000",AA_pH_Homeostasis="#252525",NH4_pH_Homeostasis="#525252",CRISPR_System="#969696",Toxin_Antitoxin="#bdbdbd",Cellular_Stress="#d9d9d9"))
    #full heatmap
    pheatmap(i,cluster_cols=F,cluster_rows=F,show_rownames=T,color=c("#b2182b","#d6604d","#f4a582","#f7f7f7","#92c5de","#4393c3","#2166ac"),clustering_method = "ward.D2",scale = "row",fontsize_col=6,fontsize_row=6,border_color=NA,annotation_legend=T,cellwidth=7,cellheight =7,annotation_colors = my_colour,annotation_row=h)
    items=unique(e$Label)
    plot_list=list()
    for (z in items){
      print(z)
      g<-f[(f$Label==z),]
      h<-unique(g[,c(1,26)])
      g<-g[,-c(26)]
      rownames(g) <- g[,1]
      g<-g[,-1]
      g<-as.matrix(g)
      col.order<- c("BVA_0216_MP05_10cm","BVA_0216_MP05_20cm","BVA_0216_MP12_10cm","BVA_0216_MP12_20cm","BVA_0815_MP10_20cm","BVA_0815_MP04_20cm","BVA_1015_MP06_20cm","BVA_1015_MP10_20cm","QUI_0216_MP05_10cm","QUI_0216_MP05_20cm","QUI_0216_MP10_10cm","QUI_0216_MP10_20cm","QUI_1015_MP10_10cm","QUI_1015_MP10_20cm","QUI_1015_MP05_10cm","QUI_1015_MP05_20cm","SJO_1015_MP02_20cm","SJO_1015_MP15_20cm","SJO_0116_MP02_10cm","SJO_0116_MP02_20cm","SJO_0116_MP15_10cm","SJO_0116_MP15_20cm","SJO_0715_MP02_20cm","SJO_0715_MP15_20cm")
      g<-g[,col.order]
      h<-data.frame(row.names = f$Complex, Metabolism = f$Label)
      g<-pheatmap(g,cluster_cols=F,cluster_rows=T,show_rownames=T,color=c("#b2182b","#d6604d","#f4a582","#f7f7f7","#92c5de","#4393c3","#2166ac"),clustering_method = "ward.D2",scale = "row",gaps_col=c(8,16),fontsize_col=6,fontsize_row=6,border_color=NA,annotation_legend=T,cellwidth=7,cellheight =7,annotation_colors = my_colour,annotation_row =h,treeheight_row = 0,main=z)
      plot_list[[z]] = g[[4]]
    }
    do.call(grid.arrange,plot_list[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)])
    #stats on gene subsets
    plot_list2=list()
    for (z in items){
      g<-f[(f$Label==z),]
      g<-g[,-c(2,26)]
      g<-melt(g)
      g<-dcast(g, variable ~ Complex)
      g$variable <- as.character(g$variable)
      g$variable[grepl("SJO", g$variable)] <- "SJO"
      g$variable[grepl("BVA", g$variable)] <- "BVA"
      g$variable[grepl("QUI", g$variable)] <- "QUI"
      h<-as.matrix(g[,-1])
      i<-vegdist(h, method='bray')
      j<-adonis2(i~variable, data=g, permutations = 999, method="bray")
      print(z)
      print(j)
    }
  }
