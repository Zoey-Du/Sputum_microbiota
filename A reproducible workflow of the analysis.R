#A reproducible workflow of the analysis

#Title: Azithromycin Exposure Induces Transient Microbial Composition Shifts and Decreasing Ability of The Airway Microbiota Resilient from PM2.5 Stress in Healthy Adults: A Randomized, Double-blind, Placebo-controlled Trial

#Authors: Sisi Du

#Date:2020/8/3

#Baseline characteristics
library(plyr)
library(reshape2)
#Table 1. Demographics of healthy volunteers
baselinedo <- baseline[baseline$Date1=='D0',]
table(baselinedo$Group,baselinedo$Gender)
table(baselinedo$Group,baselinedo$Race)
table(baselinedo$Group,baselinedo$Residence)
table(baselinedo$Group,baselinedo$Occupation)
tapply(X = baselinedo$Age,INDEX = baselinedo$Group,FUN = mean )
tapply(X = baselinedo$Height,INDEX = baselinedo$Group,FUN = mean )
tapply(X = baselinedo$Weight,INDEX = baselinedo$Group,FUN = mean )
tapply(X = baselinedo$Age,INDEX = baselinedo$Group,FUN = sd )
tapply(X = baselinedo$Height,INDEX = baselinedo$Group,FUN = sd )
tapply(X = baselinedo$Weight,INDEX = baselinedo$Group,FUN = sd )

#Difference test
baselinedo1 <- baselinedo[baselinedo$Group=='one',]
baselinedo0 <- baselinedo[baselinedo$Group=='zero',]
#Age compared by two-sided Wilcoxon rank sum test
wilcox.test(baselinedo1$Age,baselinedo0$Age)
#Height compared by two-sided Student's t-test
t.test(baselinedo1$Height,baselinedo0$Height)
#Weight Compared by two-sided Student's t-test
t.test(baselinedo1$Weight,baselinedo0$Weight)
#Residence compared by two-sided Chi-squared test with continuity correction
chisq.test(table(baselinedo$Group,baselinedo$Residence),correct = T)
#Race compared by two-sided Chi-squared test with continuity correction.
chisq.test(table(baselinedo$Group,baselinedo$Race),correct = T)
#Occupation compared by two-sided Pearson's Chi-squared test.
chisq.test(table(baselinedo$Group,baselinedo$Occupation))

#Minimal influence of procedural contaminations on sputum samples sequencing
#Figure S1a PCoA plot based on unweighted UniFrac distance between sputum samples and negative control samples.
#PERMANOVA between sputum samples and negative control samples.
library(vegan)
library(ggplot2)
library(permute)
library(lattice)
path="D:/myfunction"
setwd(path)  
source("D:/myfunction/adonis5.R")
library(lmPerm)
library(ape)
source("D:/myfunction/vegdist2.R")
otu_tree <- read.tree("D:/sputum/otus.tree")
adonisdata <-read.csv("D:/sputum/adonisam.csv",header = T,row.names = 1)
adonisgroup =read.table("D:/sputum/adonisgroup.txt", header=T, row.names= 1,sep="\t")
adonis5(adonisdata~Sample,data=adonisgroup,permutations = 999, method="du",tree=otu_tree) -> au
au

#PCoA plot
unifracpcoasam <- vegdist2(adonisdata,method = 'du',tree =otu_tree )
unifracpcoasam1 <- as.matrix(unifracpcoasam)
idx =rownames(adonisgroup) %in% colnames(unifracpcoasam1)
sub_design =adonisgroup[idx,]
unifracpcoasam1 =unifracpcoasam1[rownames(sub_design), rownames(sub_design)]
pcoasam =cmdscale(unifracpcoasam1, k=2, eig=T)
pointsams = as.data.frame(pcoasam$points) 
eigsam = pcoasam$eig
levels(sub_design$Sample)=c("neg","sam")
pointsams = cbind(pointsams, sub_design$Sample)
colnames(pointsams) = c("PC1", "PC2","Sample") 
p = ggplot(pointsams, aes(x=PC1, y=PC2, color=Sample)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigsam[1] / sum(eigsam), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigsam[2] / sum(eigsam), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Samplecontrol') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S1b PCoA plot based on weighted UniFrac distance between sputum samples and negative control samples.
#PERMANOVA between sputum samples and negative control samples.
adonis5(adonisdata~Sample,data=adonisgroup,permutations = 999, method="dw",tree=otu_tree) -> aw
aw
#PCoA plot
wunifracpcoasam <- vegdist2(adonisdata,method = 'dw',tree =otu_tree )
wunifracpcoasam1 <- as.matrix(wunifracpcoasam)
idx =rownames(adonisgroup) %in% colnames(wunifracpcoasam1)
sub_design =adonisgroup[idx,]
wunifracpcoasam1 =wunifracpcoasam1[rownames(sub_design), rownames(sub_design)]
pcoasam1 =cmdscale(wunifracpcoasam1, k=2, eig=T)
pointsam1s = as.data.frame(pcoasam1$points) 
eigsam1 = pcoasam1$eig
levels(sub_design$Sample)=c("neg","sam")
pointsam1s = cbind(pointsam1s, sub_design$Sample)
colnames(pointsam1s) = c("PC1", "PC2","Sample") 
p = ggplot(pointsam1s, aes(x=PC1, y=PC2, color=Sample)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigsam1[1] / sum(eigsam1), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigsam1[2] / sum(eigsam1), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Samplecontrol') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1) 

#Figure S1c Top 15 ZOTUs are ranked in descending order of mean relative abundance of sputum samples.
otua <-read.table("D:/sputum/otua.txt", header=T,row.names= 1, sep="\t")
otua_relative <- otua/rowSums(otua)
write.table(otua_relative,file='D:/sputum/otua_relative.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)

tapply(X = sdsam$OTU_1,INDEX = sdsam$Group,FUN = sd)
tapply(X = sdsam$OTU_3,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_2,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_4,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_7,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_5,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_6,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_10,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_8,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_13,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_15,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_11,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_9,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_14,INDEX = sdsam$Group,FUN = sd )
tapply(X = sdsam$OTU_17,INDEX = sdsam$Group,FUN = sd )

View(samcon1)
samcon$Sam <- factor(samcon$Sam)
samcon$Sam <- relevel(samcon$Sam,ref = "Sputum samples")
ggplot(samcon,aes(OR,R,fill=Phylum,ymin=R,ymax=R+SD))+
  geom_bar(stat = 'identity',width = 0.45,position = 'identity')+
  geom_errorbar(width = 0.45,color='black')+
  scale_x_continuous(limits=c(0,8),breaks=seq(0.5,7.5,by=0.5),expand=c(0,0),labels = reorder(samcon1$Lablel,samcon1$OR))+
  facet_grid(Sam~.)+
  theme_classic()+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.20),expand=c(0,0))+
  theme(axis.text.x=element_text(size=6,angle=45,color="black",hjust = 1))+
  ylab("Relative abundance (%of all sequences)")+
  xlab('')+
  theme(axis.text.y=element_text(size=6), axis.title.y=element_text(size=6))+
  theme(strip.text = element_text(colour = 'white'), strip.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.text = element_text(size=6),legend.title = element_text(size = 6))+
  theme(panel.spacing = unit(2,'lines'))+
  scale_fill_manual(values = c('#FF3333','#330000','#660000','#990000'))

#Figure S1d  Top 15 ZOTUs are ranked in descending order of mean relative abundance of control samples.
tapply(X = sdneg$OTU_36,INDEX = sdneg$Group,FUN = sd)
tapply(X = sdneg$OTU_41,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_29,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_7,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_1,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_8,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_2,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_3,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_4,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_594,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_5,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_15,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_6,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_57,INDEX = sdneg$Group,FUN = sd )
tapply(X = sdneg$OTU_794,INDEX = sdneg$Group,FUN = sd )

ggplot(consam,aes(OR,R,fill=Phylum,ymin=R,ymax=R+SD))+
  geom_bar(stat = 'identity',width = 0.45,position = 'identity')+
  geom_errorbar(width = 0.45,color='black')+
  scale_x_continuous(limits=c(0,8),breaks=seq(0.5,7.5,by=0.5),expand=c(0,0),labels = reorder(consam1$Lablel,consam1$OR))+
  facet_grid(Sam~.)+
  theme_classic()+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.35),expand=c(0,0))+
  theme(axis.text.x=element_text(size=6,angle=45,color="black",hjust = 1))+
  ylab("Relative abundance (%of all sequences)")+
  xlab('')+
  theme(axis.text.y=element_text(size=6), axis.title.y=element_text(size=6))+
  theme(strip.text = element_text(colour = 'white'), strip.background = element_rect(fill = 'white', colour = 'white'))+
  theme(legend.text = element_text(size=6),legend.title = element_text(size = 6))+
  theme(panel.spacing = unit(2,'lines'))+
  scale_fill_manual(values = c('#FF9966','#FF3333','#330000','#660000','#990000'))

#Figure S1e-f Decontam package confirms contaminants ZOTUs in our data using both the frequency and prevalence methods with threshold 0.1 and 0.5 respectively.
install.packages("phyloseqGraphTest")
BiocManager::install("phyloseq",ask = F,update = F)
BiocManager::install("decontam",ask = F,update = F)
library(phyloseq)
library(decontam)
otudec<-read.csv("D:/sputum/otu_table.csv",header = T,row.names = 1)
meta<-read.csv("D:/sputum/sampledata.csv",header = T,row.names = 1)
OTUDEC <-otu_table(otudec,taxa_are_rows = F)
META <-sample_data(meta)
phy <-phyloseq(OTUDEC,META)
sample_data(phy)$is.neg <- sample_data(phy)$SampleControl== "neg"
contamdf.both <- isContaminant(phy, conc="Quan", neg="is.neg", method="both", threshold=c(0.1,0.5))
table(contamdf.both$contaminant)
which(contamdf.both$contaminant)
head(contamdf.both)
write.table(contamdf.both,file='D:/sputum/contamdf.both.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
#Figure S1e
plot_frequency(phy, taxa_names(phy)[c(762,1344,1515,1)], conc="Quan") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")+
  theme_classic()+
  theme(axis.text.y=element_text(size=15), axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15))+
  theme(strip.text = element_text(colour = 'black',size = 30))

#Figure S1f
phy.pa <- transform_sample_counts(phy, function(abund) 1*(abund>0))
phy.pa.neg <- prune_samples(sample_data(phy.pa)$SampleControl== "neg", phy.pa)
phy.pa.pos <- prune_samples(sample_data(phy.pa)$SampleControl== "sam", phy.pa)
df.pb <- data.frame(pa.pos=taxa_sums(phy.pa.pos), pa.neg=taxa_sums(phy.pa.neg),contaminant=contamdf.both$contaminant)
write.table(df.pb,file='D:/sputum/df.pb.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
df.pb =read.table("D:/sputum/df.pb.txt", header=T, row.names= 1,sep="\t")
ggplot(data=df.pb, aes(x=pa.neg, y=pa.pos)) + 
  geom_jitter(aes(size=contamdf.both$freq,color=contaminant),alpha=1/3)+
  geom_text(aes(y = pa.pos+15, label = Name))+
  scale_size_continuous(name='Relative abundance',breaks=c(0.000,0.025,0.050,0.075),labels=c('0%','2.5%','5%','7.5%'))+
  theme_classic()+
  xlab("Prevalence (Control samples)") + ylab("Prevalence (Sputum samples)")+
  theme(axis.text.y=element_text(size=15), axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15), axis.title.x=element_text(size=15))+
  scale_color_manual(values = c('#FF3333','#330000'))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))

#Stable bacterial DNA burden in healthy volunteers’ airways after azithromycin administration
#Figure 2a Sputum bacterial DNA quantification using droplet digital PCR of the 16S rRNA gene in azithromycin group.
tapply(X=DNAzero$DNA,INDEX = DNAzero$Date1,FUN = median)
tapply(X=DNAzero$DNA,INDEX = DNAzero$Date1,FUN = quantile)
View(DNAzerot)
ggplot(DNAzerot,aes(reorder(Date1,OR),Median,group = Group,ymin=Se1,ymax=Se2))+
  geom_pointrange(stat = 'identity')+
  theme_classic()+
  geom_line()+
  geom_point(DNAzerot,mapping=aes(fill=Date2),size=18,shape=21)+
  scale_fill_manual(values = c('#CC3300','grey'))+
  scale_y_continuous(limits=c(0,2e+08))+
  xlab("Timepoints") + ylab("Bacterial burden (16S rRNA gene copies/g))")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))

#Difference test
#Wilcoxon signed-rank test
##p.adjust(method = "bonferroni",n=4)

#Figure 2b Sputum bacterial DNA quantification using droplet digital PCR of the 16S rRNA gene in placebo group.
tapply(X=DNAone$DNA,INDEX = DNAone$Date1,FUN = median)
tapply(X=DNAone$DNA,INDEX = DNAone$Date1,FUN = quantile)
ggplot(DNAonet,aes(reorder(Date1,OR),Median,group = Group,ymin=Se1,ymax=Se2))+
  geom_pointrange(stat = 'identity')+
  theme_classic()+
  geom_line()+
  geom_point(DNAonet,mapping=aes(fill=Date2),size=18,shape=21)+
  scale_fill_manual(values = c('#FF9900','grey'))+
  scale_y_continuous(limits=c(0,2e+08))+
  xlab("Timepoints") + ylab("Bacterial burden (16S rRNA gene copies/g))")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))

#Difference test
#Wilcoxon signed-rank test
##p.adjust(method = "bonferroni",n=4)

#The shifts in sputum phylogenetic diversity after antibiotic exposure
#Alpha diversity
#Figure 3a-b Richness and Shannon index in azithromycin group
tapply(X=richnesszero$Richness,INDEX = richnesszero$Date1,FUN = max)
tapply(X=richnesszero$Richness,INDEX = richnesszero$Date1,FUN = min)
tapply(X=richnesszero$Richness,INDEX = richnesszero$Date1,FUN = mean)
tapply(X=richnesszero$Richness,INDEX = richnesszero$Date1,FUN = sd)

tapply(X=shannonzero$Shannon,INDEX = shannonzero$Date1,FUN = max)
tapply(X=shannonzero$Shannon,INDEX = shannonzero$Date1,FUN = min)
tapply(X=shannonzero$Shannon,INDEX = shannonzero$Date1,FUN = mean)
tapply(X=shannonzero$Shannon,INDEX = shannonzero$Date1,FUN = sd)

#Difference test
#zero
###Richness
#D0 D4
shapiro.test(pairzeroD0D4D14$D0D4richness)
t.test(pairzeroD0D4D14$D0Richness,pairzeroD0D4D14$D4Richness,paired=T)
#D0 D14
shapiro.test(pairzeroD0D4D14$D0D14richness)
t.test(pairzeroD0D4D14$D0Richness,pairzeroD0D4D14$D14Richness,paired=T)
#D4 D14
shapiro.test(pairzeroD0D4D14$D4D14richness)
t.test(pairzeroD0D4D14$D4Richness,pairzeroD0D4D14$D14Richness,paired=T)
#D0 D30
shapiro.test(pairzeroD0D4D14D30D60$D0D30richness)
t.test(pairzeroD0D4D14D30D60$D0Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D0 D60
shapiro.test(pairzeroD0D4D14D30D60$D0D60richness)
t.test(pairzeroD0D4D14D30D60$D0Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D4 D30
shapiro.test(pairzeroD0D4D14D30D60$D4D30richness)
t.test(pairzeroD0D4D14D30D60$D4Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D4 D60
shapiro.test(pairzeroD0D4D14D30D60$D4D60richness)
t.test(pairzeroD0D4D14D30D60$D4Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D14 D30
shapiro.test(pairzeroD0D4D14D30D60$D14D30richness)
t.test(pairzeroD0D4D14D30D60$D14Richness,pairzeroD0D4D14D30D60$D30Richness,paired = T)
#D14 D60
shapiro.test(pairzeroD0D4D14D30D60$D14D60richness)
t.test(pairzeroD0D4D14D30D60$D14Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)
#D30 D60
shapiro.test(pairzeroD0D4D14D30D60$D30D60richness)
t.test(pairzeroD0D4D14D30D60$D30Richness,pairzeroD0D4D14D30D60$D60Richness,paired = T)

p.adjust(c(9.799e-07,1.659e-06,0.6268,0.0008609,0.03872,0.001103,2.586e-05,4.366e-05,6.519e-06,0.009802),method = 'bonferroni',n=10)

###Shannon
#D0 D4
shapiro.test(pairzeroD0D4D14$D0D4shannon)
t.test(pairzeroD0D4D14$D0Shannon,pairzeroD0D4D14$D4Shannon,paired=T)
#D0 D14
shapiro.test(pairzeroD0D4D14$D0D14shannon)
t.test(pairzeroD0D4D14$D0Shannon,pairzeroD0D4D14$D14Shannon,paired=T)
#D4 D14
shapiro.test(pairzeroD0D4D14$D4D14shannon)
t.test(pairzeroD0D4D14$D4Shannon,pairzeroD0D4D14$D14Shannon,paired=T)
#D0 D30
shapiro.test(pairzeroD0D4D14D30D60$D0D30shannon)
t.test(pairzeroD0D4D14D30D60$D0Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D0 D60
shapiro.test(pairzeroD0D4D14D30D60$D0D60shannon)
t.test(pairzeroD0D4D14D30D60$D0Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D4 D30
shapiro.test(pairzeroD0D4D14D30D60$D4D30shannon)
t.test(pairzeroD0D4D14D30D60$D4Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D4 D60
shapiro.test(pairzeroD0D4D14D30D60$D4D60shannon)
t.test(pairzeroD0D4D14D30D60$D4Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D14 D30
shapiro.test(pairzeroD0D4D14D30D60$D14D30shannon)
t.test(pairzeroD0D4D14D30D60$D14Shannon,pairzeroD0D4D14D30D60$D30Shannon,paired = T)
#D14 D60
shapiro.test(pairzeroD0D4D14D30D60$D14D60shannon)
t.test(pairzeroD0D4D14D30D60$D14Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
#D30 D60
shapiro.test(pairzeroD0D4D14D30D60$D30D60shannon)
t.test(pairzeroD0D4D14D30D60$D30Shannon,pairzeroD0D4D14D30D60$D60Shannon,paired = T)
p.adjust(c(3.92e-07,7.106e-06,0.0001062,0.001765,0.4486,0.001766,2.782e-05,0.4506,0.002278,0.03951),method = 'bonferroni',n=10)

#Richness plot
ggplot()+
  geom_point(richnesszero,mapping=aes(reorder(Date1,OR),Richness,group=Name,colour=Date1),size=5) +
  labs(color="Timepoints")+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme_classic()+
  geom_line(richnesszero,mapping=aes(reorder(Date1,OR),Richness,group=Name),colour='grey',alpha=0.2,size=1)+
  geom_crossbar(richnesszerot,mapping=aes(reorder(Date1,OR),Richness,ymin=Richness-SD,ymax=Richness+SD,colour=Date1),width = 0.6,size=1, position = position_dodge(0.5))+
  geom_line(richnesszeropoint, mapping=aes(reorder(Date1,OR), Richness, group = Group,colour=Date1),size=1) +
  scale_y_continuous(limits=c(0,850))+
  xlab("Timepoints") + ylab("Richness")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))+
  geom_segment(aes(x=1, y=650, xend=1, yend=670),size=1)+
  geom_segment(aes(x=2, y=650, xend=2, yend=670),size=1)+
  geom_segment(aes(x=1, y=660, xend=1.5, yend=660))+
  geom_segment(aes(x=1.5, y=660, xend=2, yend=660))+
  annotate("text", x=1.5, y=690, label="q=9.799e-06", size=5)+
  geom_segment(aes(x=3, y=650, xend=3, yend=670),size=1)+
  geom_segment(aes(x=4, y=650, xend=4, yend=670),size=1)+
  geom_segment(aes(x=3, y=660, xend=3.5, yend=660))+
  geom_segment(aes(x=3.5, y=660, xend=4, yend=660))+
  annotate("text", x=3.5, y=690, label="q=4.366e-04", size=5)+
  geom_segment(aes(x=1, y=710, xend=1, yend=730),size=1)+
  geom_segment(aes(x=3, y=710, xend=3, yend=730),size=1)+
  geom_segment(aes(x=1, y=720, xend=2.5, yend=720))+
  geom_segment(aes(x=2.5, y=720, xend=3, yend=720))+
  annotate("text", x=2, y=750, label="q=1.659e-05", size=5)+
  geom_segment(aes(x=1, y=770, xend=1, yend=790),size=1)+
  geom_segment(aes(x=4, y=770, xend=4, yend=790),size=1)+
  geom_segment(aes(x=1, y=780, xend=2.5, yend=780))+
  geom_segment(aes(x=2.5, y=780, xend=4, yend=780))+
  annotate("text", x=2.5, y=810, label="q=8.609e-03", size=5)

#Shannon index plot
ggplot()+
  geom_point(shannonzero,mapping=aes(reorder(Date1,OR),Shannon,group=Name,colour=Date1),size=5) +
  labs(color="Timepoints")+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme_classic()+
  geom_line(shannonzero,mapping=aes(reorder(Date1,OR),Shannon,group=Name),colour='grey',alpha=0.2,size=1)+
  geom_crossbar(shannonzerot,mapping=aes(reorder(Date1,OR),Shannon,ymin=Shannon-SD,ymax=Shannon+SD,colour=Date1),width = 0.6, size=1,position = position_dodge(0.5))+
  geom_line(shannonzeropoint, mapping=aes(reorder(Date1,OR), Shannon, group = Group,colour=Date1),size=1)+ 
  xlab("Timepoints") + ylab("Shannon index")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))+
  geom_segment(aes(x=1, y=7, xend=1, yend=7.1),size=1)+
  geom_segment(aes(x=2, y=7, xend=2, yend=7.1),size=1)+
  geom_segment(aes(x=1, y=7.05, xend=1.5, yend=7.05))+
  geom_segment(aes(x=1.5, y=7.05, xend=2, yend=7.05))+
  annotate("text", x=1.5, y=7.15, label="q=3.92e-06", size=5)+
  geom_segment(aes(x=2, y=7, xend=2, yend=7.1),size=1)+
  geom_segment(aes(x=3, y=7, xend=3, yend=7.1),size=1)+
  geom_segment(aes(x=2, y=7.05, xend=2.5, yend=7.05))+
  geom_segment(aes(x=2.5, y=7.05, xend=3, yend=7.05))+
  annotate("text", x=2.5, y=7.15, label="q=1.062e-03", size=5)+
  geom_segment(aes(x=1, y=7.2, xend=1, yend=7.3),size=1)+
  geom_segment(aes(x=3, y=7.2, xend=3, yend=7.3),size=1)+
  geom_segment(aes(x=1, y=7.25, xend=2.5, yend=7.25))+
  geom_segment(aes(x=2.5, y=7.25, xend=3, yend=7.25))+
  annotate("text", x=2, y=7.35, label="q=7.106e-05", size=5)+
  geom_segment(aes(x=1, y=7.4, xend=1, yend=7.5),size=1)+
  geom_segment(aes(x=4, y=7.4, xend=4, yend=7.5),size=1)+
  geom_segment(aes(x=1, y=7.45, xend=2.5, yend=7.45))+
  geom_segment(aes(x=2.5, y=7.45, xend=4, yend=7.45))+
  annotate("text", x=2.5, y=7.55, label="q=1.765e-02", size=5)

#Figure S2a-b Richness and Shannon index in placebo group
tapply(X=richnessone$Richness,INDEX = richnessone$Date1,FUN = max)
tapply(X=richnessone$Richness,INDEX = richnessone$Date1,FUN = min)
tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = mean)
tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = sd)

tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = max)
tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = min)
tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = mean)
tapply(X=shannonone$Shannon,INDEX = shannonone$Date1,FUN = sd)

#Difference test
#one
###RICHNESS
#D0 D4
shapiro.test(paironeD0D4D14$D0D4richness)
t.test(paironeD0D4D14$D0Richness,paironeD0D4D14$D4Richness,paired=T)
#D0 D14
shapiro.test(paironeD0D4D14$D0D14richness)
t.test(paironeD0D4D14$D0Richness,paironeD0D4D14$D14Richness,paired=T)
#D4 D14
shapiro.test(paironeD0D4D14$D4D14richness)
t.test(paironeD0D4D14$D4Richness,paironeD0D4D14$D14Richness,paired=T)
#D0 D30
shapiro.test(paironeD0D4D14D30$D0D30richness)
t.test(paironeD0D4D14D30$D0Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D4 D30
shapiro.test(paironeD0D4D14D30$D4D30richness)
t.test(paironeD0D4D14D30$D4Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D14 D30
shapiro.test(paironeD0D4D14D30$D14D30richness)
t.test(paironeD0D4D14D30$D14Richness,paironeD0D4D14D30$D30Richness,paired = T)
#D0 D60
shapiro.test(paironeD0D4D14D30D60$D0D60richness)
t.test(paironeD0D4D14D30D60$D0Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D4 D60
shapiro.test(paironeD0D4D14D30D60$D4D60richness)
t.test(paironeD0D4D14D30D60$D4Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D14 D60
shapiro.test(paironeD0D4D14D30D60$D14D60richness)
t.test(paironeD0D4D14D30D60$D14Richness,paironeD0D4D14D30D60$D60Richness,paired = T)
#D30 D60
shapiro.test(paironeD0D4D14D30D60$D30D60richness)
t.test(paironeD0D4D14D30D60$D30Richness,paironeD0D4D14D30D60$D60Richness,paired = T)

p.adjust(c(0.0594,0.08711,0.8442,0.09092,0.3446,0.4468,0.7597,0.1593,0.1278,0.04169),method = 'bonferroni',n=10)

###shannon
#D0 D4
shapiro.test(paironeD0D4D14$D0D4shannon)
t.test(paironeD0D4D14$D0Shannon,paironeD0D4D14$D4Shannon,paired=T)
#D0 D14
shapiro.test(paironeD0D4D14$D0D14shannon)
t.test(paironeD0D4D14$D0Shannon,paironeD0D4D14$D14Shannon,paired=T)
#D4 D14
shapiro.test(paironeD0D4D14$D4D14shannon)
t.test(paironeD0D4D14$D4Shannon,paironeD0D4D14$D14Shannon,paired=T)
#D0 D30
shapiro.test(paironeD0D4D14D30$D0D30shannon)
t.test(paironeD0D4D14D30$D0Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D4 D30
shapiro.test(paironeD0D4D14D30$D4D30shannon)
t.test(paironeD0D4D14D30$D4Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D14 D30
shapiro.test(paironeD0D4D14D30$D14D30shannon)
t.test(paironeD0D4D14D30$D14Shannon,paironeD0D4D14D30$D30Shannon,paired = T)
#D0 D60
shapiro.test(paironeD0D4D14D30D60$D0D60shannon)
t.test(paironeD0D4D14D30D60$D0Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D4 D60
shapiro.test(paironeD0D4D14D30D60$D4D60shannon)
t.test(paironeD0D4D14D30D60$D4Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D14 D60
shapiro.test(paironeD0D4D14D30D60$D14D60shannon)
t.test(paironeD0D4D14D30D60$D14Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)
#D30 D60
shapiro.test(paironeD0D4D14D30D60$D30D60shannon)
t.test(paironeD0D4D14D30D60$D30Shannon,paironeD0D4D14D30D60$D60Shannon,paired = T)

p.adjust(c(0.4779,0.2411,0.4133,0.5112,0.7783,0.7604,0.09514,0.06938,0.05103,0.09506),method = 'bonferroni',n=10)

#Richness plot
ggplot()+
  geom_point(richnessone,mapping=aes(reorder(Date1,OR),Richness,group=Name,colour=Date1),size=5) +
  labs(color="Timepoints")+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme_classic()+
  geom_line(richnessone,mapping=aes(reorder(Date1,OR),Richness,group=Name),colour='grey',alpha=0.2,size=1)+
  geom_crossbar(richnessonet,mapping=aes(reorder(Date1,OR),Richness,ymin=Richness-SD,ymax=Richness+SD,colour=Date1),width = 0.6, size=1,position = position_dodge(0.5))+
  geom_line(richnessonepoint, mapping=aes(reorder(Date1,OR), Richness, group = Group,colour=Date1),size=1) +
  scale_y_continuous(limits=c(0,850))+
  xlab("Timepoints") + ylab("Richness")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))

#Shannon index plot
ggplot()+
  geom_point(shannonone,mapping=aes(reorder(Date1,OR),Shannon,group=Name,colour=Date1),size=5) +
  labs(color="Timepoints")+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme_classic()+
  geom_line(shannonone,mapping=aes(reorder(Date1,OR),Shannon,group=Name),colour='grey',alpha=0.2,size=1)+
  geom_crossbar(shannononet,mapping=aes(reorder(Date1,OR),Shannon,ymin=Shannon-SD,ymax=Shannon+SD,colour=Date1),width = 0.6,size=1, position = position_dodge(0.5))+
  geom_line(shannononepoint, mapping=aes(reorder(Date1,OR), Shannon, group = Group,colour=Date1),size=1)+
  xlab("Timepoints") + ylab("Shannon index")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15))

#Beta diversity
pairwise.adonis5 <-function(x,factors, sim.method,tree=NULL, p.adjust.m)
  
{
  library(vegan)
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    ad = adonis5(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                   
                   factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method,tree=tree);
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  
}
#Figure 3c PCoA plot based on unweighted UniFrac distance in azithromycin group
otu_tree <- read.tree("D:/sputum/otus.tree")
adoniszero <- read.csv("D:/sputum/adoniszero.csv",header = T,row.names = 1)
samplezero =read.table("D:/sputum/samplezero.txt", header=T, row.names= 1,sep="\t")
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="du",tree=otu_tree) -> adt
adt
a <- pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="du",tree=otu_tree, p.adjust.m= "bonferroni")

unifracpcoa <- vegdist2(adoniszero,method = 'du',tree =otu_tree )
unifracpcoa1 <- as.matrix(unifracpcoa)
idx =rownames(samplezero) %in% colnames(unifracpcoa1)
sub_design =samplezero[idx,]
unifracpcoa1 =unifracpcoa1[rownames(sub_design), rownames(sub_design)]
pcoa =cmdscale(unifracpcoa1, k=2, eig=T)
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
levels(sub_design$Date1)=c("D0","D14","D30","D4","D60")
points = cbind(points, sub_design$Date1)
colnames(points) = c("PC1", "PC2","Date1") 
p = ggplot(points, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S2c PCoA plot based on unweighted UniFrac distance in placebo group
adonisone <- read.csv("D:/sputum/adonisone.csv",header = T,row.names = 1)
sampleone =read.table("D:/sputum/sampleone.txt", header=T, row.names= 1,sep="\t")
adonis5(adonisone~Date1,data=sampleone,permutations = 999, method="du",tree=otu_tree) -> adl
adl
pairwise.adonis5(adonisone,sampleone$Date1, sim.method ="du",tree=otu_tree, p.adjust.m= "bonferroni")

unifracpcoat <- vegdist2(adonisone,method = 'du',tree =otu_tree )
unifracpcoa1t <- as.matrix(unifracpcoat)
idx =rownames(sampleone) %in% colnames(unifracpcoa1t)
sub_designt =sampleone[idx,]
unifracpcoa1t =unifracpcoa1t[rownames(sub_designt), rownames(sub_designt)]
pcoat =cmdscale(unifracpcoa1t, k=2, eig=T)
pointst = as.data.frame(pcoat$points) 
eigt = pcoat$eig
levels(sub_designt$Date1)=c("D0","D14","D30","D4","D60")
pointst = cbind(pointst, sub_designt$Date1)
colnames(pointst) = c("PC1", "PC2","Date1") 
p = ggplot(pointst, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigt[1] / sum(eigt), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigt[2] / sum(eigt), digits=4), "%)", sep=""),
       title="UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure 3d Unweighted unifrac distance dissimilarity between the two groups
beta1 <- beta[beta$Date1=='D0_D4',]
beta2 <- beta[beta$Date1=='D0_D14',]
beta3 <- beta[beta$Date1=='D0_D30',]
beta4 <- beta[beta$Date1=='D0_D60',]

tapply(X = beta1$unifrac,INDEX = beta1$Group,FUN = max )
tapply(X = beta1$unifrac,INDEX = beta1$Group,FUN = min )
tapply(X = beta1$unifrac,INDEX = beta1$Group,FUN = mean )
tapply(X = beta1$unifrac,INDEX = beta1$Group,FUN = sd )

tapply(X = beta2$unifrac,INDEX = beta2$Group,FUN = max )
tapply(X = beta2$unifrac,INDEX = beta2$Group,FUN = min )
tapply(X = beta2$unifrac,INDEX = beta2$Group,FUN = mean )
tapply(X = beta2$unifrac,INDEX = beta2$Group,FUN = sd )

tapply(X = beta3$unifrac,INDEX = beta3$Group,FUN = max )
tapply(X = beta3$unifrac,INDEX = beta3$Group,FUN = min )
tapply(X = beta3$unifrac,INDEX = beta3$Group,FUN = mean )
tapply(X = beta3$unifrac,INDEX = beta3$Group,FUN = sd )

tapply(X = beta4$unifrac,INDEX = beta4$Group,FUN = max )
tapply(X = beta4$unifrac,INDEX = beta4$Group,FUN = min )
tapply(X = beta4$unifrac,INDEX = beta4$Group,FUN = mean )
tapply(X = beta4$unifrac,INDEX = beta4$Group,FUN = sd )  

#Difference test
#t-test:t.test(paired=F)
#p.adjust(method = "bonferroni",n=4)
shapiro.test(beta1$unifrac)
t.test(beta1$unifrac~beta1$Group,paired=F)
p.adjust(c(8.28e-09,1.868e-06,0.008359,0.07482),method = "bonferroni",n=4)

t.test(beta4$weight_unifrac~beta4$Group,paired=F)
p.adjust(c(0.002116,0.3121,0.5818,0.7259),method = "bonferroni",n=4)


unifrac$Date1 <- factor(unifrac$Date1)
unifracpoint$Date1 <- factor(unifracpoint$Date1)
unifract$Date1 <- factor(unifract$Date1)
unifrac$Date1 <- relevel(unifrac$Date1,ref = 'D0_D4 3.312e-08')
unifracpoint$Date1 <- relevel(unifracpoint$Date1,ref = 'D0_D4 3.312e-08')
unifract$Date1 <- relevel(unifract$Date1,ref = 'D0_D4 3.312e-08')
levels(unifrac$Date1) <- c("D0_D4\nq=3.312e-08", "D0_D14\nq=3.7.472e-06", "D0_D30\nq=3.3436e-02","D0_D60")
levels(unifract$Date1) <- c("D0_D4\nq=3.312e-08", "D0_D14\nq=3.7.472e-06", "D0_D30\nq=3.3436e-02","D0_D60")
levels(unifracpoint$Date1) <- c("D0_D4\nq=3.312e-08", "D0_D14\nq=3.7.472e-06", "D0_D30\nq=3.3436e-02","D0_D60")

ggplot()+
  geom_point(unifrac,mapping=aes(Group,unifrac,colour=Group),size=5) +
  geom_crossbar(unifract,mapping=aes(Group,unifrac,ymin=unifrac-SD,ymax=unifrac+SD,colour=Group),width = 0.6, size=1,position = position_dodge(0.5))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  geom_line(unifracpoint, mapping=aes(Group,unifrac, group = Position,colour=Group),size=1) +
  facet_grid(.~Date1 ) +
  theme(strip.text = element_text(colour = 'black',size = 15), strip.background = element_rect(fill = 'white',size = rel(2)))+
  xlab("") + ylab("UniFrac distance dissimilarity")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_blank())+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

#Figure 3e PCoA plot based on weighted UniFrac distance in azithromycin group
adonis5(adoniszero~Date1,data=samplezero,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
pairwise.adonis5(adoniszero,samplezero$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")

wunifracpcoa <- vegdist2(adoniszero,method = 'dw',tree =otu_tree )
wunifracpcoa1 <- as.matrix(wunifracpcoa)
idx =rownames(samplezero) %in% colnames(wunifracpcoa1)
sub_design =samplezero[idx,]
wunifracpcoa1 =wunifracpcoa1[rownames(sub_design), rownames(sub_design)]
pcoa =cmdscale(wunifracpcoa1, k=2, eig=T)
points = as.data.frame(pcoa$points) 
eig = pcoa$eig
levels(sub_design$Date1)=c("D0","D14","D30","D4","D60")
points = cbind(points, sub_design$Date1)
colnames(points) = c("PC1", "PC2","Date1") 
p = ggplot(points, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)

#Figure S2e PCoA plot based on weighted UniFrac distance in placebo group
adonis5(adonisone~Date1,data=sampleone,permutations = 999, method="dw",tree=otu_tree) -> adt
adt
pairwise.adonis5(adonisone,sampleone$Date1, sim.method ="dw",tree=otu_tree, p.adjust.m= "bonferroni")
wunifracpcoat <- vegdist2(adonisone,method = 'dw',tree =otu_tree )
wunifracpcoa1t <- as.matrix(wunifracpcoat)
idx =rownames(sampleone) %in% colnames(wunifracpcoa1t)
sub_designt =sampleone[idx,]
wunifracpcoa1t =wunifracpcoa1t[rownames(sub_designt), rownames(sub_designt)]
pcoat =cmdscale(wunifracpcoa1t, k=2, eig=T)
pointst = as.data.frame(pcoat$points) 
eigt = pcoat$eig
levels(sub_designt$Date1)=c("D0","D14","D30","D4","D60")
pointst = cbind(pointst, sub_designt$Date1)
colnames(pointst) = c("PC1", "PC2","Date1") 
p = ggplot(pointst, aes(x=PC1, y=PC2, color=Date1)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigt[1] / sum(eigt), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigt[2] / sum(eigt), digits=4), "%)", sep=""),
       title="Weight UniFrac distance PCoA",color='Timepoints') + theme_classic()+
  scale_color_manual(values = c('#003399','#FF3333','#FF9966','#CC6600','#330000'))+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

p + stat_ellipse(level=0.68,size=1)


#Figure 3f Weighted unifrac distance dissimilarity between the two groups
tapply(X = beta1$weight_unifrac,INDEX = beta1$Group,FUN = max )
tapply(X = beta1$weight_unifrac,INDEX = beta1$Group,FUN = min )
tapply(X = beta1$weight_unifrac,INDEX = beta1$Group,FUN = mean )
tapply(X = beta1$weight_unifrac,INDEX = beta1$Group,FUN = sd )

tapply(X = beta2$weight_unifrac,INDEX = beta2$Group,FUN = max )
tapply(X = beta2$weight_unifrac,INDEX = beta2$Group,FUN = min )
tapply(X = beta2$weight_unifrac,INDEX = beta2$Group,FUN = mean )
tapply(X = beta2$weight_unifrac,INDEX = beta2$Group,FUN = sd )

tapply(X = beta3$weight_unifrac,INDEX = beta3$Group,FUN = max )
tapply(X = beta3$weight_unifrac,INDEX = beta3$Group,FUN = min )
tapply(X = beta3$weight_unifrac,INDEX = beta3$Group,FUN = mean )
tapply(X = beta3$weight_unifrac,INDEX = beta3$Group,FUN = sd )

tapply(X = beta4$weight_unifrac,INDEX = beta4$Group,FUN = max )
tapply(X = beta4$weight_unifrac,INDEX = beta4$Group,FUN = min )
tapply(X = beta4$weight_unifrac,INDEX = beta4$Group,FUN = mean )
tapply(X = beta4$weight_unifrac,INDEX = beta4$Group,FUN = sd )

#Difference test
#t-test:t.test(paired=F)
#p.adjust(method = "bonferroni",n=4)

wunifrac$Date1 <- factor(wunifrac$Date1)
wunifracpoint$Date1 <- factor(wunifracpoint$Date1)
wunifract$Date1 <- factor(wunifract$Date1)
wunifrac$Date1 <- relevel(wunifrac$Date1,ref = 'D0_D4')
wunifracpoint$Date1 <- relevel(wunifracpoint$Date1,ref = 'D0_D4')
wunifract$Date1 <- relevel(wunifract$Date1,ref = 'D0_D4')
levels(wunifrac$Date1) <- c("D0_D4\nq=0.0085", "D0_D14", "D0_D30","D0_D60")
levels(wunifract$Date1) <- c("D0_D4\nq=0.0085", "D0_D14", "D0_D30","D0_D60")
levels(wunifracpoint$Date1) <- c("D0_D4\nq=0.0085", "D0_D14", "D0_D30","D0_D60")


ggplot()+
  geom_point(wunifrac,mapping=aes(Group,weight_unifrac,colour=Group),size=5) +
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  geom_crossbar(wunifract,mapping=aes(Group,weight_unifrac,ymin=weight_unifrac-SD,ymax=weight_unifrac+SD,colour=Group),size=1,width = 0.6, position = position_dodge(0.5))+
  geom_line(wunifracpoint, mapping=aes(Group,weight_unifrac, group = Position,colour=Group),size=1) +
  facet_grid(.~Date1) +
  theme(strip.text = element_text(colour = 'black',size = 15), strip.background = element_rect(fill = 'white',size = rel(2)))+
  xlab("") + ylab("Weight UniFrac distance dissimilarity")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_blank())+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

#Microbial taxonomic variation during 60 days’ follow-up 
otur <-read.table("D:/sputum/otur.txt", header=T,row.names= 1, sep="\t")
otur_relative <- otur/rowSums(otur)
write.table(otur_relative,file='D:/sputum/otur_relative.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)
View(otur)

#Difference test in detection rate
#McNemar-Bowker test 
###ZERO
#D0 VS D4
P <- c(rep(0,ncol(D01)))
for(i in 5:ncol(D01))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D01[,i]==1&D4[,i]==1))
matrix1[1,2]=length(which(D01[,i]==1&D4[,i]==0))
matrix1[2,1]=length(which(D01[,i]==0&D4[,i]==1))
matrix1[2,2]=length(which(D01[,i]==0&D4[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P1.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D14
P <- c(rep(0,ncol(D01)))
for(i in 5:ncol(D01))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D01[,i]==1&D14[,i]==1))
matrix1[1,2]=length(which(D01[,i]==1&D14[,i]==0))
matrix1[2,1]=length(which(D01[,i]==0&D14[,i]==1))
matrix1[2,2]=length(which(D01[,i]==0&D14[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P2.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)


#D0 VS D30
P <- c(rep(0,ncol(D03)))
for(i in 5:ncol(D03))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D03[,i]==1&D30[,i]==1))
matrix1[1,2]=length(which(D03[,i]==1&D30[,i]==0))
matrix1[2,1]=length(which(D03[,i]==0&D30[,i]==1))
matrix1[2,2]=length(which(D03[,i]==0&D30[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P3.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D60
P <- c(rep(0,ncol(D03)))
for(i in 5:ncol(D03))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D03[,i]==1&D60[,i]==1))
matrix1[1,2]=length(which(D03[,i]==1&D60[,i]==0))
matrix1[2,1]=length(which(D03[,i]==0&D60[,i]==1))
matrix1[2,2]=length(which(D03[,i]==0&D60[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P4.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

###ONE
#D0 VS D4
P <- c(rep(0,ncol(D010)))
for(i in 5:ncol(D010))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D010[,i]==1&D41[,i]==1))
matrix1[1,2]=length(which(D010[,i]==1&D41[,i]==0))
matrix1[2,1]=length(which(D010[,i]==0&D41[,i]==1))
matrix1[2,2]=length(which(D010[,i]==0&D41[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P11.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D14
P <- c(rep(0,ncol(D010)))
for(i in 5:ncol(D010))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D010[,i]==1&D141[,i]==1))
matrix1[1,2]=length(which(D010[,i]==1&D141[,i]==0))
matrix1[2,1]=length(which(D010[,i]==0&D141[,i]==1))
matrix1[2,2]=length(which(D010[,i]==0&D141[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P21.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)


#D0 VS D30
P <- c(rep(0,ncol(D0d30)))
for(i in 5:ncol(D0d30))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D0d30[,i]==1&D301[,i]==1))
matrix1[1,2]=length(which(D0d30[,i]==1&D301[,i]==0))
matrix1[2,1]=length(which(D0d30[,i]==0&D301[,i]==1))
matrix1[2,2]=length(which(D0d30[,i]==0&D301[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P31.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)

#D0 VS D60
P <- c(rep(0,ncol(D0d60)))
for(i in 5:ncol(D0d60))
{ matrix1 <- matrix(data=0,nrow=2,ncol=2,byrow=FALSE,dimnames=NULL)
matrix1[1,1]=length(which(D0d60[,i]==1&D601[,i]==1))
matrix1[1,2]=length(which(D0d60[,i]==1&D601[,i]==0))
matrix1[2,1]=length(which(D0d60[,i]==0&D601[,i]==1))
matrix1[2,2]=length(which(D0d60[,i]==0&D601[,i]==0))
y= mcnemar.test(matrix1)
P[i] <- y$p.value}
write.table(P,file='D:/sputum/P41.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)


#Difference test in relative abundance
#linear mixed model 
library(lmerTest)
glmgenzero$Date1 <- factor(glmgenzero$Date1)
glmgenzero$Date1 <- relevel(glmgenzero$Date1,ref = 'D30')
fit1=lmer(Phocaeicola~Date1+(1|Name),data = glmgenzero)
summary(fit1)
glmgenzero$Date1 <- relevel(glmgenzero$Date1,ref = 'D30')
#Figure 4 Heat map
library(ComplexHeatmap)
library(GetoptLong)
library(dendextend)
library(circlize)
library(grid)
eg = read.table("D:/sputum/pvalue.txt",header=T,sep='\t',row.names=1)
SSP= read.table("D:/sputum/SSP.txt",header=T,sep='\t',row.names=1)
det = read.table("D:/sputum/det.txt",header=T,sep='\t',row.names=1)
det <- as.matrix(det)
detz = apply(det, 1, scale)
rownames(detz) = colnames(det)
write.table(detz,file='D:/sputum/detz.txt',quote=FALSE,sep='\t',row.names=F,col.names = T)
detz = read.table("D:/sputum/detz.txt",header=T,sep='\t',row.names=1)
detzt <- t(detz)
H1=Heatmap(detzt,name="Detection rate\n     Z-score", col = colorRamp2(c(-1,1), c("white","#0000CC")),show_column_names=T, row_names_side='left',show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(250,'mm'),heatmap_legend_param=list(legend_direction='horizontal'),split=SSP$Family)
rel = read.table("D:/sputum/rel.txt",header=T,sep='\t',row.names=1)
rel <- as.matrix(rel)
relz = apply(rel, 1, scale)
rownames(relz) = colnames(rel)
relzt <- t(relz)
pvalue = row_anno_points(eg$D0_D14)
ha_mix_right = HeatmapAnnotation(p = pvalue,which = 'row',width  = unit( 2, "cm"))
H2=Heatmap(relzt,name="Relative abundance\n     Z-score", col = colorRamp2(c(-1,1), c("white","#FF0000")),right_annotation = ha_mix_right,show_column_names=T, show_row_names = F,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(250,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

H3=add_heatmap(H1,H2,direction =  "horizontal")

###Figure S4
SSP= read.table("D:/sputum/SSPone.txt",header=T,sep='\t',row.names=1)
det = read.table("D:/sputum/detone.txt",header=T,sep='\t',row.names=1)
det <- as.matrix(det)
H1=Heatmap(det,name="Detection rate", col = colorRamp2(c(0.1,1), c("white","#0000CC")),show_column_names=T, row_names_side='left',show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(300,'mm'),heatmap_legend_param=list(legend_direction='horizontal'),split=SSP$Family)
rel = read.table("D:/sputum/relone.txt",header=T,sep='\t',row.names=1)
rel <- as.matrix(rel)
H2=Heatmap(rel,name="Relative abundance", col = colorRamp2(c(0.0001,0.01), c("white","#FF0000")),show_column_names=T, show_row_names = F,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(40, "mm"),height = unit(300,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

H3=add_heatmap(H1,H2,direction =  "horizontal")

#Figure S3a Sankey plots to describe the number of ZOTUs whose detection rate shift or return to the baseline level across the five timepoints. 
library(riverplot)
library(RColorBrewer) 
edges = data.frame(ID =c('F','G','H','I','J','K','L','M','N','O','P','Q','A1','A2','A3','A4'),N1 = c('A','A','Cno','A','Dno','A','A','B','B','C','C','Cno','Cno','D','D','Dno'),  
                   N2 = c('Eno','E','D','D','E','C','B','Cno','C','Dno','D','Eno','E','Eno','E','Eno'),  
                   Value = c(17,25,11,32,3,38,182,16,166,5,199,1,4,8,234,2),  
                   stringsAsFactors = F)  
nodes = data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = FALSE)  
nodes$x = c(1,5,7,3,5,7,9,9)  
nodes$y = c(3,-2,-2,9,7,5,0,4)  
rownames(nodes) = nodes$ID  
head(nodes) 
palette = paste0(brewer.pal(9,"Blues"), "70") 
palette1 = paste0(brewer.pal(9,"BuPu"), "70")
j <- makeRiver(nodes, edges,node_labels= c( A= "D0", Cno='No return',Dno='No return',B= "D4", C= "D14",D= "D30",Eno='No return',E= "D60" ),
               node_styles= list( A= list( col= '#0066CC'),Cno=list(col='#990066'),Dno=list(col='#990066'),B=list(col= "#0066CC"),C=list(col= "#0066CC"),D=list(col= "#0066CC"),Eno=list(col='#990066'),E=list(col= "#0066CC")), 
               edge_styles = list(F=list(col=palette1[6]),G=list(col=palette[5]),H=list(col=palette[5]),I=list(col=palette[5]),J=list(col=palette[5]),
                                  K=list(col=palette[5]), L=list(col=palette[5]),M=list(col=palette1[6]),N=list(col=palette[5]),O=list(col=palette1[6]),
                                  P=list(col=palette[5]),Q=list(col=palette1[6]),A1=list(col=palette[5]),A2=list(col=palette1[6]),A3=list(col=palette[5]),A4=list(col=palette1[6])))

ds <- default.style()
ds[["edgecol"]] <- "col"
ds[['col']]<-"white"

plot(j,yscale=0.02,srt='0',default_style=ds)


#Figure S3b Sankey plots to describe the number of ZOTUs whose relative abundance shift or return to the baseline level across the five timepoints. 

edges1 = data.frame(ID =c('F','G','H','I','J','K','L','M','N','O','P','Q','A1','A2','A3','A4'),N1 = c('A','A','Cno','A','Dno','A','A','B','B','C','C','Cno','Cno','D','D','Dno'),  
                    N2 = c('Eno','E','D','D','E','C','B','Cno','C','Dno','D','Eno','E','Eno','E','Eno'),  
                    Value = c(12,11,16,12,19,32,227,31,196,23,205,6,9,14,219,4),  
                    stringsAsFactors = F)  
nodes1 = data.frame(ID = unique(c(edges1$N1, edges1$N2)), stringsAsFactors = FALSE) 
nodes1$x = c(1,5,7,3,5,7,9,9)  
nodes1$y = c(3,-2,-2,9,8,9,0,4)  
rownames(nodes1) = nodes1$ID  
head(nodes1)  
palette2 = paste0(brewer.pal(9,"Reds"), "70") #####
palette1 = paste0(brewer.pal(9,"BuPu"), "70")####
j1 <- makeRiver(nodes1, edges1,node_labels= c( A= "D0", Cno='No return',Dno='No return',B= "D4", C= "D14",D= "D30",Eno='No return',E= "D60" ),
                node_styles= list( A= list( col= '#FF3333'),Cno=list(col='#990066'),Dno=list(col='#990066'),B=list(col= "#FF3333"),C=list(col= "#FF3333"),D=list(col= "#FF3333"),Eno=list(col='#990066'),E=list(col= "#FF3333")), 
                edge_styles = list(F=list(col=palette1[6]),G=list(col=palette2[5]),H=list(col=palette2[5]),I=list(col=palette2[5]),J=list(col=palette2[5]),
                                   K=list(col=palette2[5]), L=list(col=palette2[5]),M=list(col=palette1[6]),N=list(col=palette2[5]),O=list(col=palette1[6]),
                                   P=list(col=palette2[5]),Q=list(col=palette1[6]),A1=list(col=palette2[5]),A2=list(col=palette1[6]),A3=list(col=palette2[5]),A4=list(col=palette1[6])))

ds <- default.style()
ds[["edgecol"]] <- "col"
ds[['col']]<-'white'
plot(j1,yscale=0.02,srt='0',default_style=ds)

#Figure S3c For the taxa whose relative abundance have been of no significant difference from the D0 level at D14, novel shifts in their relative abundance at timepoint D30 are shown at ZOTU levels.
OTUen1$log2R <- log2(OTUen1$R)
OTUen2$log2R <- log2(OTUen2$R)
OTUen1$P <- factor(OTUen1$P)
levels(OTUen1$P) <- c("D0 VS D30", "D14 VS D30")
OTUen2$P <- factor(OTUen2$P)
levels(OTUen2$P) <- c("D0 VS D30", "D14 VS D30")
ggplot() +
  geom_point(OTUen1,mapping=aes(log2R,reorder(Name,log2R),color=P),size=10)+
  scale_color_manual(values = c('#FF0000','#0000CC'))+
  geom_point(OTUen2,mapping=aes(log2R,reorder(Name,log2R),colour=P),size=5)+
  scale_color_manual(values = c('#FF0000','#0000CC'))+
  scale_x_continuous(limits=c(-3,4),breaks = seq(-3,4,by=1))+
  guides(color=guide_legend(title=NULL))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log2(D30/D14)") + ylab("ZOTU")+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size = 15),legend.key.size = unit(1,'cm'))

#Figure S3d For the taxa whose relative abundance have been of no significant difference from the D0 level at D14, novel shifts in their relative abundance at timepoint D30 are shown at family and genus levels.
gen1$log2R <- log2(gen1$R)
gen2$log2R <- log2(gen2$R)
gen1$P <- factor(gen1$P)
levels(gen1$P) <- c("D0 VS D30", "D14 VS D30")
gen2$P <- factor(gen2$P)
levels(gen2$P) <- c("D0 VS D30", "D14 VS D30")
ggplot() +
  geom_point(gen1,mapping=aes(log2R,reorder(Genus,log2R),color=P),size=10)+
  scale_color_manual(values = c('#FF0000','#0000CC'))+
  geom_point(gen2,mapping=aes(log2R,reorder(Genus,log2R),colour=P),size=5)+
  scale_color_manual(values = c('#FF0000','#0000CC'))+
  scale_x_continuous(limits=c(-3,3),breaks = seq(-3,3,by=1))+
  guides(color=guide_legend(title=NULL))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log2(D30/D14)") + ylab("")+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x=element_text(size=15))+
  theme(legend.text = element_text(size = 15),legend.key.size = unit(1,'cm'))

#Attenuated resilience in response to environmental factors disturbance
#Figure S5 The daily mean air pollution data and temperature and humidity during the sputum samples collecting.
library(gtable)
library(grid) 
ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  # Get the location of the plot panel in g1.
  # These are used later when transformed elements of g2 are put back into g1
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from g2
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  # Put the transformed label on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from g2 (axis line, tick marks, and tick mark labels)
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of g1
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}
factordata$MonthDate <- as.Date(factordata$MonthDate)
factordata1$start <- as.Date(factordata1$start)
factordata1$end <- as.Date(factordata1$end)
factordata1$Date1 <- factor(factordata1$Date1)
levels(factordata1$Date1)=c("Day2","Day1","Day")
g1 <- ggplot() + 
  geom_rect(mapping = aes(xmin=start,xmax=end,fill=Date1),ymin=-Inf,ymax=Inf,alpha=0.3,data = factordata1)+
  guides(fill=F)+
  scale_fill_manual(values = c('white','#FFCCCC','#FFFF99'))+
  geom_line(factordata, mapping=aes(MonthDate, AQI,colour = 'AQI'),size=2)+ 
  geom_line(factordata, mapping=aes(MonthDate, PM10,colour = 'PM10'),size=2)+ 
  geom_line(factordata, mapping=aes(MonthDate, PM2.5,colour = 'PM2.5'),size=2)+
  labs(y = 'PM2.5/PM10/AQI(μg/m3)',x='') + 
  scale_x_date(date_labels="%m/%d/%Y",date_breaks  ="7 day")+
  theme_classic()+
  theme(axis.text.x=element_text(size=15,angle=30,color="black",hjust = 1))+
  theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 15))+
  scale_color_manual(values =c('#FF0000','#990066','#000066') )+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'),legend.position = 'left')%+replace%
  theme(panel.background = element_rect(fill = "transparent",colour = NA))
g2 <- ggplot() + 
  geom_line(factordata, mapping=aes(MonthDate, Temperature,colour='Temperature'),size=2)+
  geom_line(factordata, mapping=aes(MonthDate, Humidity1,colour='Humidity'),size=2) +
  labs(y = 'Temperature/Humidity(℃/%)',x='')+
  scale_color_manual(values = c('#FF9900','#336600'))+
  scale_x_date(date_labels="%m/%d/%Y",date_breaks  ="7 day")+
  scale_y_continuous(limits=c(-10,15),breaks=seq(-10,15,by=5),labels = c('-10','-5','0','5/50','10/100','15/150'))+
  theme_classic()+
  theme(axis.text.y=element_text(size=15),axis.title.y = element_text(size = 15))+
  theme(axis.text.x=element_text(size=15,angle=30,color="black",hjust = 1),legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm')) 
ggplot2.two_y_axis(g2, g1)

#Difference test
tapply(X = factordata$PM2.5,INDEX = factordata$Date1,FUN = median )
tapply(X = factordata$PM2.5,INDEX = factordata$Date1,FUN = quantile )
tapply(X = factordata$PM10,INDEX = factordata$Date1,FUN = median )
tapply(X = factordata$PM10,INDEX = factordata$Date1,FUN = quantile )
tapply(X = factordata$AQI,INDEX = factordata$Date1,FUN = median )
tapply(X = factordata$AQI,INDEX = factordata$Date1,FUN = quantile )
tapply(X = factordata$Humidity,INDEX = factordata$Date1,FUN = median )
tapply(X = factordata$Humidity,INDEX = factordata$Date1,FUN = quantile )
tapply(X = factordata$Temperature,INDEX = factordata$Date1,FUN = median )
tapply(X = factordata$Temperature,INDEX = factordata$Date1,FUN = quantile )
kruskal.test(PM2.5~Date1, data=factordata) 
kruskal.test(PM10~Date1, data=factordata) 
kruskal.test(AQI~Date1, data=factordata)
kruskal.test(Humidity~Date1, data=factordata)
kruskal.test(Temperature~Date1, data=factordata)
posthoc.kruskal.nemenyi.test(PM2.5~Date1, data=factordatar,dist="Tukey",method="bonferroni")
p.adjust(c(2.0e-08,4.0e-14,0.0514,0.1330,0.0018,0.0096,0.0081,6.2e-11,1.5e-10,0.9994),method = "bonferroni",n=10)
posthoc.kruskal.nemenyi.test(PM10~Date1, data=factordatar,dist="Tukey")
p.adjust(c(0.0016,3.4e-14,0.9984,0.4869,1.1e-06,0.0052,0.3051,4.3e-14,5.4e-11,0.6701),method = "bonferroni",n=10)
tapply(X = factordatar$AQI,INDEX = factordatar$Date1,FUN = IQR )
posthoc.kruskal.nemenyi.test(AQI~Date1, data=factordatar,dist="Tukey")
p.adjust(c(8.3e-06,3.9e-14,0.157,0.378,1.6e-06,0.062,0.035,9.2e-14,1.3e-13,0.997),method = "bonferroni",n=10)

#Figure 5 Distance-based redundancy analysis (db-RDA) based on weighted UniFrac distance. 
#Figure 5a The relation of the environmental factors, antibiotic event and the sputum microbiota variation at D4. 
library('ggrepel')
library('dplyr')
library(devtools)
library(usethis)
library(vegan)
library(permute)
library(lattice)
library(ade4)
path="D:/myfunction"
setwd(path)  
source("D:/myfunction/vegdist2.R")
source("D:/myfunction/adonis5.R")
library(ape)
otu_tree <- read.tree("D:/sputum/otus.tree")

otu1 <- read.csv("D:/sputum/RDAD0D4.csv", row.names=1)
RDAfactor1 <- read.csv("D:/sputum/factorD0D4.csv", row.names=1)
weiu1 <- vegdist2(otu1, method = 'dw',tree = otu_tree)
pcoa1 <- cmdscale(weiu1, k = nrow(otu1) - 1, eig = TRUE, add = TRUE)
pcoa_site1 <- pcoa1$point
rda_db1 <- rda(pcoa_site1~Humidity+PM10+Temperature+Date1,RDAfactor1, scale = FALSE)
rda_db_test1 <- anova(rda_db1, permutations = 999)
rda_db_test_axis1 <- anova(rda_db1, by = 'term', permutations = 999)
rda_db_test_axis1$`Pr(>F)` <- p.adjust(rda_db_test_axis1$`Pr(>F)`, method = 'bonferroni')
RsquareAdj(rda(pcoa_site1~ Humidity, data=RDAfactor1))$adj.r.squared
RsquareAdj(rda(pcoa_site1~ Humidity+Date1, data=RDAfactor1))$adj.r.squared
rda_db_forward_vp <- varpart(pcoa_site1, RDAfactor1['Date1'], RDAfactor1['Humidity'])
rda_db_forward_vp
anova(rda(pcoa_site1, RDAfactor1['Date1']), permutations = 999)
anova(rda(pcoa_site1, RDAfactor1['Humidity']), permutations = 999)

total1=summary(rda_db1)
sites1 = total1$sites[,1:2] %>% data.frame() %>% merge(RDAfactor1[,1:2],by = 'row.names') 
biplot1 =total1$biplot[,1:2] %>% data.frame()
biplot1$env = c("Humidity","PM10","Temperature","Antibiotic event") 
sites1$Date1 = factor(sites1$Date1, levels = c('D0','D4')) 
levels(sites1$Date1) = c('Antibiotic-before','Antibiotic-after') 
ggplot(sites1, aes(x = RDA1, y =RDA2, color = Date1)) +
  theme_classic() + 
  scale_color_manual(values = c('#000099','#FF6633'))+ 
  labs( x= "RDA1 (10.25%)", y = "RDA2 (2.85%)", color = '') + 
  guides(color = guide_legend(override.aes = list(size=5)))+ 
  geom_hline(yintercept=0, linetype=2,color='grey') + 
  geom_vline(xintercept=0, linetype=2,color='grey') + 
  geom_point(size = 5) +
  stat_ellipse(show.legend = F,size=1) + 
  geom_segment(data = biplot1,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(1/2, 'picas')), lwd = 1,
               colour = "#993300") +  
  
  geom_text_repel(data = biplot1, 
                  aes(x=RDA1,y=RDA2,label=env),
                  size= 5, fontface='bold',color='black')+ 
  
  theme(legend.position = 'right',
        legend.background = element_blank(),
        legend.text = element_text(face = 'bold',color='black',size=15),
        axis.title = element_text(face = 'bold',color='black',size=15),
        axis.text = element_text(face = 'bold',color='black',size=15),
        panel.grid = element_blank())

#Figure 5b The relation of the environmental factors, antibiotic event and the sputum microbiota variation at D30. 
otu <- read.csv("D:/sputum/RDAD0D30.csv", row.names=1)
RDAfactor <- read.csv("D:/sputum/factorD0D30.csv", row.names=1)
weiu <- vegdist2(otu, method = 'dw',tree = otu_tree)
pcoa <- cmdscale(weiu, k = nrow(otu) - 1, eig = TRUE, add = TRUE)
pcoa_site <- pcoa$point
rda_dbb <- rda(pcoa_site~Humidity+PM10+PM2.5+Date1, RDAfactor, scale = FALSE)
rda_db_testb <- anova(rda_dbb, permutations = 999)
rda_db_test_axisb <- anova(rda_dbb, by = 'term', permutations = 999)
rda_db_test_axisb$`Pr(>F)` <- p.adjust(rda_db_test_axisb$`Pr(>F)`, method = 'bonferroni')
RsquareAdj(rda(pcoa_site~ PM2.5 , data=RDAfactor))$adj.r.squared
RsquareAdj(rda(pcoa_site~ PM2.5+Date1 , data=RDAfactor))$adj.r.squared
total=summary(rda_dbb)
sites = total$sites[,1:2] %>% data.frame() %>% merge(RDAfactor[,1:2],by = 'row.names') 
biplot =total$biplot[,1:2] %>% data.frame()
biplot$env = c("Humidity","PM10","PM2.5","Antibiotic event")
sites$Date1 = factor(sites$Date1, levels = c('D0','D30')) 
levels(sites$Date1) = c('Antibiotic-before','Antibiotic-after') 

ggplot(sites, aes(x = RDA1, y =RDA2, color = Date1)) +
  theme_classic() + 
  scale_color_manual(values = c('#000099','#FF6633'))+
  labs( x= "RDA1 (9.03%)", y = "RDA2 (5.05%)", color = '') + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_hline(yintercept=0, linetype=2,color='grey') + 
  geom_vline(xintercept=0, linetype=2,color='grey') + 
  geom_point(size = 5) +
  stat_ellipse(show.legend = F,size=1) + 
  geom_segment(data = biplot,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(1/2, 'picas')), lwd = 1,
               colour = "#993300") +  
  
  geom_text_repel(data = biplot, 
                  aes(x=RDA1,y=RDA2,label=env),
                  size= 5, color='black')+ 
  
  theme(legend.position = 'right',
        legend.background = element_blank(),
        legend.text = element_text(color='black',size=15),
        axis.title = element_text(color='black',size=15),
        axis.text = element_text(color='black',size=15),
        panel.grid = element_blank())

###Figure S6 The relation of the environmental factors and the sputum microbiota variation in placebo group.
#D0 VS D4
otu1 <- read.csv("D:/sputum/RDAD0D4one.csv", row.names=1)
RDAfactor1 <- read.csv("D:/sputum/factorD0D4one.csv", row.names=1)
weiu1 <- vegdist2(otu1, method = 'dw',tree = otu_tree)
pcoa1 <- cmdscale(weiu1, k = nrow(otu1) - 1, eig = TRUE, add = TRUE)
pcoa_site1 <- pcoa1$point
rda_db1 <- rda(pcoa_site1~Humidity+PM10+Temperature,RDAfactor1, scale = FALSE)
vif.cca(rda_db1)
rda_db_test1 <- anova(rda_db1, permutations = 999)
rda_db_test_axis1 <- anova(rda_db1, by = 'term', permutations = 999)
rda_db_test_axis1$`Pr(>F)` <- p.adjust(rda_db_test_axis1$`Pr(>F)`, method = 'bonferroni')
RsquareAdj(rda(pcoa_site1~ Humidity, data=RDAfactor1))$adj.r.squared
RsquareAdj(rda(pcoa_site1~ Humidity+PM10, data=RDAfactor1))$adj.r.squared
RsquareAdj(rda(pcoa_site1~ Humidity+PM10+Temperature, data=RDAfactor1))$adj.r.squared


total1=summary(rda_db1)
sites1 = total1$sites[,1:2] %>% data.frame() %>% merge(RDAfactor1[,1:2],by = 'row.names') 
biplot1 =total1$biplot[,1:2] %>% data.frame()
biplot1$env = c("Humidity","PM10","Temperature") 
sites1$Date1 = factor(sites1$Date1, levels = c('D0','D4')) 
levels(sites1$Date1) = c('Antibiotic-before','Antibiotic-after') 
ggplot(sites1, aes(x = RDA1, y =RDA2, color = Date1)) +
  theme_classic() + 
  scale_color_manual(values = c('#000099','#FF6633'))+ 
  labs( x= "RDA1 (2.55%)", y = "RDA2 (1.97%)", color = '') + 
  guides(color = guide_legend(override.aes = list(size=5)))+ 
  geom_hline(yintercept=0, linetype=2,color='grey') + 
  geom_vline(xintercept=0, linetype=2,color='grey') + 
  geom_point(size = 5) +
  stat_ellipse(show.legend = F,size=1) + 
  geom_segment(data = biplot1,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(1/2, 'picas')), lwd = 1,
               colour = "#993300") +  
  
  geom_text_repel(data = biplot1, 
                  aes(x=RDA1,y=RDA2,label=env),
                  size= 5, fontface='bold',color='black')+ 
  
  theme(legend.position = 'right',
        legend.background = element_blank(),
        legend.text = element_text(face = 'bold',color='black',size=15),
        axis.title = element_text(face = 'bold',color='black',size=15),
        axis.text = element_text(face = 'bold',color='black',size=15),
        panel.grid = element_blank())

#D0 VS D30
otu <- read.csv("D:/sputum/RDAD0D30one.csv", row.names=1)
RDAfactor <- read.csv("D:/sputum/factorD0D30one.csv", row.names=1)
weiu <- vegdist2(otu, method = 'dw',tree = otu_tree)
pcoa <- cmdscale(weiu, k = nrow(otu) - 1, eig = TRUE, add = TRUE)
pcoa_site <- pcoa$point
rda_dbb <- rda(pcoa_site~Humidity+PM10+PM2.5, RDAfactor, scale = FALSE)
vif.cca(rda_dbb)
rda_db_testb <- anova(rda_dbb, permutations = 999)
rda_db_test_axisb <- anova(rda_dbb, by = 'term', permutations = 999)
rda_db_test_axisb$`Pr(>F)` <- p.adjust(rda_db_test_axisb$`Pr(>F)`, method = 'bonferroni')
RsquareAdj(rda(pcoa_site~ Humidity , data=RDAfactor))$adj.r.squared
RsquareAdj(rda(pcoa_site~ Humidity+PM10 , data=RDAfactor))$adj.r.squared
RsquareAdj(rda(pcoa_site~ Humidity+PM10+PM2.5 , data=RDAfactor))$adj.r.squared

total=summary(rda_dbb)
sites = total$sites[,1:2] %>% data.frame() %>% merge(RDAfactor[,1:2],by = 'row.names') 
biplot =total$biplot[,1:2] %>% data.frame()
biplot$env = c("Humidity","PM10","PM2.5")
sites$Date1 = factor(sites$Date1, levels = c('D0','D30')) 
levels(sites$Date1) = c('Antibiotic-before','Antibiotic-after') 

ggplot(sites, aes(x = RDA1, y =RDA2, color = Date1)) +
  theme_classic() + 
  scale_color_manual(values = c('#000099','#FF6633'))+
  labs( x= "RDA1 (3.48%)", y = "RDA2 (1.76%)", color = '') + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_hline(yintercept=0, linetype=2,color='grey') + 
  geom_vline(xintercept=0, linetype=2,color='grey') + 
  geom_point(size = 5) +
  stat_ellipse(show.legend = F,size=1) + 
  geom_segment(data = biplot,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(1/2, 'picas')), lwd = 1,
               colour = "#993300") +  
  
  geom_text_repel(data = biplot, 
                  aes(x=RDA1,y=RDA2,label=env),
                  size= 5, color='black')+ 
  
  theme(legend.position = 'right',
        legend.background = element_blank(),
        legend.text = element_text(color='black',size=15),
        axis.title = element_text(color='black',size=15),
        axis.text = element_text(color='black',size=15),
        panel.grid = element_blank())


#The effects of azithromycin on airway microbial interactions within the microbial community network
#Figure 6a The network analysis in azithromycin group. 
library(igraph)
library(reshape2)
library(bindrcpp)
library(VennDiagram)
library(grid)
library(futile.logger)
library(ggthemes)
library(shiny)
library(agricolae)
library(vegan)
library(LearnBayes)
library(dplyr)
library(psych)
library(sqldf)
library(gsubfn)
library(proto)
library(RSQLite)
library(digest)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(parallel)
library(Biobase)
library(IRanges)
library(S4Vectors)
library(impute)
library(GO.db)
library(preprocessCore)
library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)
library(multtest)
source("D:/myfunction/matrix2igraph.R")
source("D:/myfunction/net_pro.R")
source("D:/myfunction/node_pro.R")
#D0 network
otu_sample_file <- "D:/sputum/netpre0.txt"
otu_tax_file<-"D:/sputum/netpretax0.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre0 <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre0<-t(otunetpre0)
otunetpretax0 <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax0<-as.data.frame(otunetpretax0[colnames(otunetpre0),])
otu_abundancepre0 <- colSums(otunetpre0)
otu_pronetpre0 <- cbind(otunetpretax0,otu_abundancepre0)
igraphpre0<-matrix2igraph(otunetpre0,r.threshold,p.threshold)
igraphpre0.weight <- E(igraphpre0)$weight
E(igraphpre0)$weight <- NA
sum(igraphpre0.weight>0)# number of postive correlation
sum(igraphpre0.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results0")
#The topological properties of microbial networks
netpro_result<-net_pro(igraphpre0)
write.csv(netpro_result,"D:/sputum/network_results0/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraphpre0)
write.csv(nodepro_result,"D:/sputum/network_results0/igraphpre.node.pro.csv")

#D4 network
otu_sample_file <- "D:/sputum/net0D4.txt"
otu_tax_file<-"D:/sputum/net0D4tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph0D4<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph0D4.weight <- E(igraph0D4)$weight
E(igraph0D4)$weight <- NA
sum(igraph0D4.weight>0)# number of postive correlation
sum(igraph0D4.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results0D4")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph0D4)
write.csv(netpro_result,"D:/sputum/network_results0D4/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph0D4)
write.csv(nodepro_result,"D:/sputum/network_results0D4/igraphpre.node.pro.csv")

#D14 network
otu_sample_file <- "D:/sputum/net0D14.txt"
otu_tax_file<-"D:/sputum/net0D14tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph0D14<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph0D14.weight <- E(igraph0D14)$weight
E(igraph0D14)$weight <- NA
sum(igraph0D14.weight>0)# number of postive correlation
sum(igraph0D14.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results0D14")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph0D14)
write.csv(netpro_result,"D:/sputum/network_results0D14/igraphpre.network.pro.csv")

#The topological properties of microbial notes
nodepro_result<-node_pro(igraph0D14)
write.csv(nodepro_result,"D:/sputum/network_results0D14/igraphpre.node.pro.csv")

#D30 network
otu_sample_file <- "D:/sputum/net0D30.txt"
otu_tax_file<-"D:/sputum/net0D30tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph0D30<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph0D30.weight <- E(igraph0D30)$weight
E(igraph0D30)$weight <- NA
sum(igraph0D30.weight>0)# number of postive correlation
sum(igraph0D30.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results0D30")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph0D30)
write.csv(netpro_result,"D:/sputum/network_results0D30/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph0D30)
write.csv(nodepro_result,"D:/sputum/network_results0D30/igraphpre.node.pro.csv")

#D60 network
otu_sample_file <- "D:/sputum/net0D60.txt"
otu_tax_file<-"D:/sputum/net0D60tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph0D60<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph0D60.weight <- E(igraph0D60)$weight
E(igraph0D60)$weight <- NA
sum(igraph0D60.weight>0)# number of postive correlation
sum(igraph0D60.weight<0)# number of negative correlation


dir.create("D:/sputum/network_results0D60")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph0D60)
write.csv(netpro_result,"D:/sputum/network_results0D60/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph0D60)
write.csv(nodepro_result,"D:/sputum/network_results0D60/igraphpre.node.pro.csv")

#Figure S8a The network analysis in placebo group.
#D0 network
otu_sample_file <- "D:/sputum/netpre1.txt"
otu_tax_file<-"D:/sputum/netpretax1.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre1 <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre1<-t(otunetpre1)
otunetpretax1 <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax1<-as.data.frame(otunetpretax1[colnames(otunetpre1),])
otu_abundancepre1 <- colSums(otunetpre1)
otu_pronetpre1 <- cbind(otunetpretax1,otu_abundancepre1)
igraphpre1<-matrix2igraph(otunetpre1,r.threshold,p.threshold)

igraphpre1.weight <- E(igraphpre1)$weight
E(igraphpre1)$weight <- NA
sum(igraphpre1.weight>0)# number of postive correlation
sum(igraphpre1.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results1")
#The topological properties of microbial networks
netpro_result<-net_pro(igraphpre1)
write.csv(netpro_result,"D:/sputum/network_results1/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraphpre1)
write.csv(nodepro_result,"D:/sputum/network_results1/igraphpre.node.pro.csv")

#D4 network
otu_sample_file <- "D:/sputum/net1D4.txt"
otu_tax_file<-"D:/sputum/net1D4tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph1D4<-matrix2igraph(otunetpre,r.threshold,p.threshold)

igraph1D4.weight <- E(igraph1D4)$weight

E(igraph1D4)$weight <- NA

sum(igraph1D4.weight>0)# number of postive correlation
sum(igraph1D4.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results1D4")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph1D4)
write.csv(netpro_result,"D:/sputum/network_results1D4/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph1D4)
write.csv(nodepro_result,"D:/sputum/network_results1D4/igraphpre.node.pro.csv")


#D14 network
otu_sample_file <- "D:/sputum/net1D14.txt"
otu_tax_file<-"D:/sputum/net1D14tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph1D14<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph1D14.weight <- E(igraph1D14)$weight

E(igraph1D14)$weight <- NA
sum(igraph1D14.weight>0)# number of postive correlation
sum(igraph1D14.weight<0)# number of negative correlation


dir.create("D:/sputum/network_results1D14")

#The topological properties of microbial networks
netpro_result<-net_pro(igraph1D14)
write.csv(netpro_result,"D:/sputum/network_results1D14/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph1D14)
write.csv(nodepro_result,"D:/sputum/network_results1D14/igraphpre.node.pro.csv")

#D30
otu_sample_file <- "D:/sputum/net1D30.txt"
otu_tax_file<-"D:/sputum/net1D30tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)

otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph1D30<-matrix2igraph(otunetpre,r.threshold,p.threshold)

igraph1D30.weight <- E(igraph1D30)$weight
E(igraph1D30)$weight <- NA
sum(igraph1D30.weight>0)# number of postive correlation
sum(igraph1D30.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results1D30")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph1D30)
write.csv(netpro_result,"D:/sputum/network_results1D30/igraphpre.network.pro.csv")

#The topological properties of microbial notes
nodepro_result<-node_pro(igraph1D30)
write.csv(nodepro_result,"D:/sputum/network_results1D30/igraphpre.node.pro.csv")

#D60
otu_sample_file <- "D:/sputum/net1D60.txt"
otu_tax_file<-"D:/sputum/net1D60tax.txt"
r.threshold=0.6
p.threshold=0.05
size=3
gcol=2
glab=3
otunetpre <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
otunetpre<-t(otunetpre)
otunetpretax <- read.table(otu_tax_file, header=T, row.names= 1, sep="\t", comment.char="")
otunetpretax<-as.data.frame(otunetpretax[colnames(otunetpre),])
otu_abundancepre <- colSums(otunetpre)
otu_pronetpre <- cbind(otunetpretax,otu_abundancepre)
igraph1D60<-matrix2igraph(otunetpre,r.threshold,p.threshold)
igraph1D60.weight <- E(igraph1D60)$weight

E(igraph1D60)$weight <- NA
sum(igraph1D60.weight>0)# number of postive correlation
sum(igraph1D60.weight<0)# number of negative correlation

dir.create("D:/sputum/network_results1D60")
#The topological properties of microbial networks
netpro_result<-net_pro(igraph1D60)
write.csv(netpro_result,"D:/sputum/network_results1D60/igraphpre.network.pro.csv")
#The topological properties of microbial notes
nodepro_result<-node_pro(igraph1D60)
write.csv(nodepro_result,"D:/sputum/network_results1D60/igraphpre.node.pro.csv")

#Figure 6b The number of shared edges between network D0 and timepoints D4, D14, D30, D60 in azithromycin

D40 <- length(E(intersection(igraphpre0,igraph0D4)))
D140 <- length(E(intersection(igraphpre0,igraph0D14)))
D300 <- length(E(intersection(igraphpre0,igraph0D30)))
D600 <- length(E(intersection(igraphpre0,igraph0D60)))

edgenumpre0 <- length(E(igraphpre0))
edgenum0D4<-length(E(igraph0D4))
edgenum0D14<-length(E(igraph0D14))
edgenum0D30<-length(E(igraph0D30))
edgenum0D60<-length(E(igraph0D60))

###D0 D4
grid.newpage()
draw.pairwise.venn(area1=edgenumpre0,area2=edgenum0D4,cross.area=D40
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

###D0 D14
grid.newpage()
draw.pairwise.venn(area1=edgenumpre0,area2=edgenum0D14,cross.area=D140
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

####D0 D30
grid.newpage()
draw.pairwise.venn(area1=edgenumpre0,area2=edgenum0D30,cross.area=D300
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

####D0 D60
grid.newpage()
draw.pairwise.venn(area1=edgenumpre0,area2=edgenum0D60,cross.area=D600
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=360, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

#Figure S8b The number of shared edges between network D0 and timepoints D4, D14, D30, D60 in placebo

D4pre<-length(E(intersection(igraphpre1,igraph1D4)))
D14pre<-length(E(intersection(igraphpre1,igraph1D14)))
D30pre<-length(E(intersection(igraphpre1,igraph1D30)))
D60pre<-length(E(intersection(igraphpre1,igraph1D60)))

edgenum1D4<-length(E(igraph1D4))
edgenum1D14<-length(E(igraph1D14))
edgenum1D30<-length(E(igraph1D30))
edgenum1D60<-length(E(igraph1D60))
edgenumpre1<-length(E(igraphpre1))
###D0 D4
grid.newpage()
draw.pairwise.venn(area1=edgenumpre1,area2=edgenum1D4,cross.area=D4pre
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=180, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

###D0 D14
grid.newpage()
draw.pairwise.venn(area1=edgenumpre1,area2=edgenum1D14,cross.area=D14pre
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=180, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

####D0 D30
grid.newpage()
draw.pairwise.venn(area1=edgenumpre1,area2=edgenum1D30,cross.area=D30pre
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=180, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

####D0 D60
grid.newpage()
draw.pairwise.venn(area1=edgenumpre1,area2=edgenum1D60,cross.area=D60pre
                   ,lwd=1,lty=1
                   ,col=c('#FF3399','#009999'),fill=c('#FF3399','#33FF33')
                   ,cat.col=c('#FF3399','#009999')
                   ,rotation.degree=180, scaled = FALSE,cex = 2,cat.cex =1, cat.pos = c(400, 200))

#Figure 6c The closeness centralization of shared nodes between network D0 and timepoints D4, D14, D30, D60 in azithromycin
library(ComplexHeatmap)
library(circlize)
zeroD4 = read.table("D:/sputum/zeroD4.txt",header=T,sep='\t',row.names=1)
zeroD4 <- as.matrix(zeroD4)
Heatmap(zeroD4,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

zeroD14 = read.table("D:/sputum/zeroD14.txt",header=T,sep='\t',row.names=1)
zeroD14 <- as.matrix(zeroD14)
Heatmap(zeroD14,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

zeroD30 = read.table("D:/sputum/zeroD30.txt",header=T,sep='\t',row.names=1)
zeroD30 <- as.matrix(zeroD30)
Heatmap(zeroD30,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

zeroD60 = read.table("D:/sputum/zeroD60.txt",header=T,sep='\t',row.names=1)
zeroD60 <- as.matrix(zeroD60)
Heatmap(zeroD60,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

#Figure S8c The closeness centralization of shared nodes between network D0 and timepoints D4, D14, D30, D60 in placebo
oneD4 = read.table("D:/sputum/oneD4.txt",header=T,sep='\t',row.names=1)
oneD4 <- as.matrix(oneD4)
Heatmap(oneD4,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

oneD14 = read.table("D:/sputum/oneD14.txt",header=T,sep='\t',row.names=1)
oneD14 <- as.matrix(oneD14)
Heatmap(oneD14,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

oneD30 = read.table("D:/sputum/oneD30.txt",header=T,sep='\t',row.names=1)
oneD30 <- as.matrix(oneD30)
Heatmap(oneD30,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

oneD60 = read.table("D:/sputum/oneD60.txt",header=T,sep='\t',row.names=1)
oneD60 <- as.matrix(oneD60)
Heatmap(oneD60,name='Rank of the closeness',col = colorRamp2(c(0,0.05), c("#FFFF99",'#3366FF')), show_column_names = F,show_row_names = T,cluster_rows = FALSE, cluster_columns = FALSE,width =  unit(150, "mm"),height = unit(40,'mm'),heatmap_legend_param=list(legend_direction='horizontal'))

#Figure S7 The topological properties of microbial networks in azithromycin and placebo groups.
library(ggplot2)
ggplot()+
  geom_point(tra1,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra1,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Vertices (n)")

ggplot()+
  geom_point(tra2,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra2,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Edges (n)")

ggplot()+
  geom_point(tra3,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra3,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Positive edges (n)")

ggplot()+
  geom_point(tra4,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra4,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Negative edges (n)")

ggplot()+
  geom_point(tra5,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra5,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Average degree")

ggplot()+
  geom_point(tra6,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra6,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Connectance")

ggplot()+
  geom_point(tra7,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra7,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Average clustering coefficient")

ggplot()+
  geom_point(tra8,mapping = aes(reorder(Date1,OR),Parameter,colour=Group),size=5)+
  geom_line(tra8,mapping=aes(reorder(Date1,OR),Parameter,group=Group,colour=Group))+
  scale_color_manual(values = c('#CC3300','#FF9900'))+
  theme_classic()+
  xlab("Timepoints") + ylab("Average centralization closeness")

