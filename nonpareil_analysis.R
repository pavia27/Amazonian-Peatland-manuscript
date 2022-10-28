rm(list = ls(all.names = TRUE))
library(Nonpareil)
library(gridExtra)
library(ggplot2)
setwd("")
all.filesQUI<-c("0216_QUI02_MP05_10_1.fastq_output.npo","1015_Q2_SP10_10_R1.fastq_output.npo","0216_QUI02_MP05_20_MG_1.fastq_output.npo","1015_Q2_SP10_20_R1.fastq_output.npo","0216_QUI02_MP10_10_MG_1.fastq_output.npo","1015_Q2_SP5_10_R1.fastq_output.npo","0216_QUI02_MP10_20_MG_1.fastq_output.npo","1015_Q2_SP5_20_R1.fastq_output.npo")
Nonpareil.curve.batch(all.filesQUI)
npsQUI<-Nonpareil.set(all.filesQUI)
print(npsQUI)
# Extract Nd diversity index
summary(npsQUI)[,"diversity"]
# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(npsQUI)[,"LRstar"]/1e9
# Show current coverage (as %)
summary(npsQUI)[,"C"]*100
# Predict coverage for a sequencing effort of 10Gbp
sapply(npsQUI$np.curves, predict, 10e9)
npsQUI.1<-Nonpareil.set(all.filesQUI, plot.opts = list(plot.observed =F, model.lwd = 1.2, legend.opts = FALSE,xlim=c(1e+6, 1e+15)))
Nonpareil.legend(npsQUI.1,x = 'bottomright',cex = 0.7)
###
all.filesBVA<-c('0815_BV2_10_20-1.fastq_output.npo','0815_BV2_4_20-1.fastq_output.npo','1015_BV01_MP06_20-1.fastq_output.npo','1015_BV01_MP10_20-1.fastq_output.npo','0216_BV02_MP05_10-1.fastq_output.npo','0216_BV02_MP05_20-1.fastq_output.npo','0216_BV02_MP12_10-1.fastq_output.npo','0216_BV02_MP12_20-1.fastq_output.npo')
Nonpareil.curve.batch(all.filesBVA)
npsBVA<-Nonpareil.set(all.filesBVA)
print(npsBVA)
# Extract Nd diversity index
summary(npsBVA)[,"diversity"]
# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(npsBVA)[,"LRstar"]/1e9
# Show current coverage (as %)
summary(npsBVA)[,"C"]*100
# Predict coverage for a sequencing effort of 10Gbp
sapply(npsBVA$np.curves, predict, 10e9)
npsBVA.1<-Nonpareil.set(all.filesBVA, plot.opts = list(plot.observed =F, model.lwd = 1.2, legend.opts = FALSE,xlim=c(1e+6, 1e+15)))
Nonpareil.legend(npsBVA.1,x = 'bottomright',cex = 0.7)
###
all.filesSJO<-c("0116_SJ02_MP02_10-1.fastq_output.npo","0116_SJ02_MP02_20-1.fastq_output.npo","0116_SJ02_MP15_10-1.fastq_output.npo","0116_SJ02_MP15_20-1.fastq_output.npo",'0715_SJ02_MP02_20-1.fastq_output.npo','0715_SJ02_MP15_20-1.fastq_output.npo','1015_SJ02_MP02_20-1.fastq_output.npo','1015_SJ02_MP15_20-1.fastq_output.npo')
Nonpareil.curve.batch(all.filesSJO)
npsSJO<-Nonpareil.set(all.filesSJO)
print(npsSJO)
# Extract Nd diversity index
summary(npsSJO)[,"diversity"]
# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(npsSJO)[,"LRstar"]/1e9
# Show current coverage (as %)
summary(npsSJO)[,"C"]*100
# Predict coverage for a sequencing effort of 10Gbp
sapply(npsSJO$np.curves, predict, 10e9)
npsSJO.1<-Nonpareil.set(all.filesSJO, mainplot.opts = list(plot.observed =F, model.lwd = 1.2, legend.opts = FALSE,xlim=c(1e+6, 1e+15)))
Nonpareil.legend(npsSJO.1,x = 'bottomright',cex = 0.7)

par(mfrow=c(3,1))
plot(npsBVA.1, main = "Buena Vista", legend.opts = T, xlim=c(1e+8, 1e+14))
Nonpareil.legend(npsBVA,x = 'bottomright',cex = 0.7)
plot(npsQUI.1, main = "Quistococha", legend.opts = FALSE, xlim=c(1e+8, 1e+14))
Nonpareil.legend(npsQUI,x = 'bottomright',cex = 0.7)
plot(npsSJO.1, main = "San Jorge", legend.opts = FALSE, xlim=c(1e+8, 1e+14))
Nonpareil.legend(npsSJO,x = 'bottomright',cex = 0.7)

plot_list<-list(npsBVA.1,npsQUI.1)
do.call(grid.arrange,plot_list)

a<-as.data.frame(summary(npsSJO)[,"diversity"])
b<-as.data.frame(summary(npsBVA)[,"diversity"])
c<-as.data.frame(summary(npsQUI)[,"diversity"])
colnames(a)<-"Diveristy"
colnames(b)<-"Diveristy"
colnames(c)<-"Diveristy"
d<-rbind(a,b)
e<-rbind(d,c)
e$site<-c("SJO","SJO","SJO","SJO","SJO","SJO","SJO","SJO","BVA","BVA","BVA","BVA","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","QUI","QUI","QUI","QUI")
ggplot(e, aes(x=site, y=Diveristy,fill=site))+geom_boxplot()+scale_fill_manual(values=c("#238443", "#fdae61", "#a50026"))+theme_classic()+theme(axis.text.x = element_text(angle=45, hjust = 1))
model<-aov(Diveristy ~ site, data = e)
summary(model)
TukeyHSD(model)
a<-as.data.frame(summary(npsSJO)[,"C"]*100)
b<-as.data.frame(summary(npsBVA)[,"C"]*100)
c<-as.data.frame(summary(npsQUI)[,"C"]*100)
colnames(a)<-"Coverage"
colnames(b)<-"Coverage"
colnames(c)<-"Coverage"
d<-rbind(a,b)
e<-rbind(d,c)
e$site<-c("SJO","SJO","SJO","SJO","SJO","SJO","SJO","SJO","BVA","BVA","BVA","BVA","BVA","BVA","BVA","BVA","QUI","QUI","QUI","QUI","QUI","QUI","QUI","QUI")
model<-aov(Coverage ~ site, data = e)
summary(model)
TukeyHSD(model)
ggplot(e, aes(x=site, y=Coverage,fill=site))+geom_boxplot()+scale_fill_manual(values=c("#238443", "#fdae61", "#a50026"))+theme_classic()+theme(axis.text.x = element_text(angle=45, hjust = 1))
}