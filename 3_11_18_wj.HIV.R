# raw data analysis was done on 1_23_2017 R code, cd4_igg_d0 and microbiome data paper (final).r is the final data analysis pulled from "1_23_2017 R code" for the first submission.
# This code is based on cd4_igg_d0 and microbiome data paper (final).r and the comments from the reviewers.
#this code was used for the paper submitted to Frontiers and the revisions from the reviewr
setwd("~/R_directory/Wei Jiang/")
library(vegan)
library(ecodist)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(scatterplot3d)
library(lattice)
library(reshape2)
library(RSvgDevice)
library(reshape)
library(indicspecies)
library(devEMF)
library(tidyr)
library(car)

#citation
citation()
citation(package="vegan")

#use the method to remove the absolute value from the original absolute value.

#remove the previous list
rm(list=ls())


#### first, lets read in data clinical and bacteria microbiome
otu.raw <- read.table(file="./2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/bacteria/120215JW515F-pr.fasta.otus.fa.OTU.txt", sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)
otu.raw <- otu.raw[,(7:42)]
colSums (otu.raw, na.rm = FALSE, dims = 1) 
otu.raw <- as.data.frame(t(otu.raw)) 
#remove the absolute values in negative control
neg <- otu.raw[(35:36),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw-negmean
otu <- otu2[-(35:36),]
otu[otu < 0] <- 0 

# calculate simpson diversity
otu$simpson <- diversity(otu,"simp")
otu$shannon <- diversity(otu, index = "shannon")
otu$Sample_ID <- rownames(otu)
a.div <- otu[,c("Sample_ID","simpson","shannon")]

#test two clinical files.
a <- read.csv(file="sample information-microbiome analysis.csv")
b <- read.csv(file="updated table I.csv") 
ab <- inner_join (a,b) %>% filter(Sample_ID %in% a.div$Sample_ID)
ab$X16sdna
ab$X16S.DNA..copies.ul.
ab$X16S.Watercont
ab$LPS.plasma
ab$LPS.Plasma.pg.ml.
ab$cd4.igg.d0
ab$CD4.autoIgG..ng.ml.


##rarefy to 100
otu.per <- prop.table(as.matrix(otu[,(1:2413)]), margin = 1) %>%data.frame()
otu.per <- otu.per*100
rowSums(otu.per)

## taxonomy based on OTUs
taxa <- read.csv(file="taxa.csv")

# removed count
negmean <- t(negmean)
negmean <- data.frame(negmean)
negmean$OTU <- row.names(negmean)
removed <- inner_join(negmean, taxa) %>% filter(negmean>0)
a <- data.frame(table(removed$kingdom))
# phyla level of water contamination.
a <- table(removed$phyla) %>% data.frame() %>% filter(Freq>0) %>% mutate (percent= Freq/439)
a <- table(removed$class) %>% data.frame()
a <- table(removed$order) %>% data.frame()%>%dplyr::rename(negative_control=Freq)


# subjects count
subjectmean <- colMeans(otu.raw[(1:34),])
subjectmean <-  data.frame(subjectmean)
subjectmean$OTU <- row.names(subjectmean)
subject <- inner_join(subjectmean,taxa)%>% filter(subjectmean>0)

b <- data.frame(table(subject$kingdom))
# phyla level of water contamination.
b <- table(subject$phyla) %>% data.frame() %>% filter(Freq>0) %>% mutate (percent= Freq/2408)
b <- table(subject$class) %>% data.frame() %>% mutate (percent= Freq/2408)
b <- table(subject$order) %>% data.frame()%>%dplyr::rename(subject=Freq)
c <- full_join(a,b)%>%dplyr::rename(order=Var1)
write.csv(c, file="order level comparision.csv")

a <- table(removed$family) %>% data.frame()%>%dplyr::rename(negative_control=Freq)
b <- table(subject$family) %>% data.frame()%>%dplyr::rename(subject=Freq)
c <- full_join(a,b)%>%dplyr::rename(family=Var1)
write.csv(c, file="family level comparision.csv")

a <- table(removed$genus) %>% data.frame()%>%dplyr::rename(negative_control=Freq)
b <- table(subject$genus) %>% data.frame()%>%dplyr::rename(subject=Freq)
c <- full_join(a,b)%>%dplyr::rename(genus=Var1)
c <- inner_join(c, taxa)
write.csv(c, file="genus level comparision.csv")

a <- table(removed$class) %>% data.frame()%>%dplyr::rename(negative_control=Freq)
b <- table(subject$class) %>% data.frame()%>%dplyr::rename(subject=Freq)
c <- full_join(a,b)%>%dplyr::rename(class=Var1)
write.csv(c, file="class level comparision.csv")

# average otu per subject
c <- rowSums(otu.raw[1:34,] != 0)%>% data.frame()
mean(c$.)
sd(c$.)

# average otu per negative control
c <- rowSums(otu.raw[35:36,] != 0)%>% data.frame()
mean(c$.)
sd(c$.)

rm(a, ab, b, c, neg, negmean,  subjectmean, otu2)

# L7 after remove OTU
otu2 <- otu[,(1:2413)]
otu2 <- t(otu2)
otu2 <- data.frame(otu2)
otu2$OTU <- rownames(otu2)

L7 <- inner_join(taxa[,c("OTU","species")],otu2)
colSums(L7[,(3:36)])
L7 <- L7[,(2:36)] %>% group_by(species) %>% summarise_each(funs(sum))
colSums(L7[,(2:35)])
L7.per <- cbind(L7[,1], prop.table(as.matrix(L7[,-1]), margin = 2))%>% data.frame()
L7.per[,(2:35)] <- L7.per[,(2:35)]*100
colSums(L7.per[,2:35])
rownames(L7.per) <- L7.per$species
L7.per <- L7.per[,-1]
L7.per <- t(L7.per) %>%data.frame()
rownames(L7) <- L7$species
L7.t <- L7[,2:35]
rownames(L7.t) <- row.names(L7)
L7 <- t(L7.t) %>% data.frame()

# L6 after remove OTU
L6 <- inner_join(taxa[,c("OTU","genus")],otu2)
colSums(L6[,(3:36)])
L6 <- L6[,(2:36)] %>% group_by(genus) %>% summarise_each(funs(sum))
colSums(L6[,(2:35)])
L6.per <- cbind(L6[,1], prop.table(as.matrix(L6[,-1]), margin = 2))%>% data.frame()
L6.per[,(2:35)] <- L6.per[,(2:35)]*100
colSums(L6.per[,2:35])
rownames(L6.per) <- L6.per$genus
L6.per <- L6.per[,-1]
L6.per <- t(L6.per) %>%data.frame()
rownames(L6) <- L6$genus
L6.t <- L6[,2:35]
rownames(L6.t) <- row.names(L6)
L6 <- t(L6.t)
L6 <- data.frame(L6)

# L5 after remove OTU
L5 <- inner_join(taxa[,c("OTU","family")],otu2)
colSums(L5[,(3:36)])
L5 <- L5[,(2:36)] %>% group_by(family) %>% summarise_each(funs(sum))
colSums(L5[,(2:35)])
L5.per <- cbind(L5[,1], prop.table(as.matrix(L5[,-1]), margin = 2))%>% data.frame()
L5.per[,(2:35)] <- L5.per[,(2:35)]*100
colSums(L5.per[,2:35])
rownames(L5.per) <- L5.per$family
L5.per <- L5.per[,-1]
L5.per <- t(L5.per) %>%data.frame()
rownames(L5) <- L5$family
L5.t <- L5[,2:35]
rownames(L5.t) <- row.names(L5)
L5 <- t(L5.t) %>% data.frame()

# L4 after remove OTU
L4 <- inner_join(taxa[,c("OTU","order")],otu2)
colSums(L4[,(3:36)])
L4 <- L4[,(2:36)] %>% group_by(order) %>% summarise_each(funs(sum))
colSums(L4[,(2:35)])
L4.per <- cbind(L4[,1], prop.table(as.matrix(L4[,-1]), margin = 2))%>% data.frame()
L4.per[,(2:35)] <- L4.per[,(2:35)]*100
colSums(L4.per[,2:35])
rownames(L4.per) <- L4.per$order
L4.per <- L4.per[,-1]
L4.per <- t(L4.per) %>%data.frame()
rownames(L4) <- L4$order
L4.t <- L4[,2:35]
rownames(L4.t) <- row.names(L4)
L4 <- t(L4.t) %>% data.frame()

# L3 after remove OTU
L3 <- inner_join(taxa[,c("OTU","class")],otu2)
colSums(L3[,(3:36)])
L3 <- L3[,(2:36)] %>% group_by(class) %>% summarise_each(funs(sum))
colSums(L3[,(2:35)])
L3.per <- cbind(L3[,1], prop.table(as.matrix(L3[,-1]), margin = 2))%>% data.frame()
L3.per[,(2:35)] <- L3.per[,(2:35)]*100
colSums(L3.per[,2:35])
rownames(L3.per) <- L3.per$class
L3.per <- L3.per[,-1]
L3.per <- t(L3.per) %>% data.frame()
rownames(L3) <- L3$class
L3.t <- L3[,-1]
rownames(L3.t) <- rownames(L3)
L3 <- t(L3.t)
L3 <- data.frame(L3)

# L2 after remove OTU
L2 <- inner_join(taxa[,c("OTU","phyla")],otu2)
colSums(L2[,(3:36)])
L2 <- L2[,(2:36)] %>% group_by(phyla) %>% summarise_each(funs(sum))
colSums(L2[,(2:35)])
L2.per <- cbind(L2[,1], prop.table(as.matrix(L2[,-1]), margin = 2))%>% data.frame()
L2.per[,(2:35)] <- L2.per[,(2:35)]*100
colSums(L2.per[,2:35])
rownames(L2.per) <- L2.per$phyla
L2.per <- L2.per[,-1]
L2.per <- t(L2.per) %>% data.frame()
colnames(L2.per) <- L2$phyla
rownames(L2) <- L2$phyla
L2.t <- L2[,-1]
rownames(L2.t) <- rownames(L2)
L2 <- t(L2.t)
L2 <- data.frame(L2)

rm(L3.t, L2.t, L4.t, L5.t, L6.t, L7.t)

clinical <- read.csv(file="updated table I.csv", sep=",", header=T)
setdiff(otu$Sample_ID, clinical$Sample_ID)
setdiff(clinical$Sample_ID, otu$Sample_ID)

a.div <- inner_join(a.div, clinical)

#demo
demo <-read.csv(file="sample information-microbiome analysis.csv") %>% filter(Sample_ID %in% a.div$Sample_ID) %>% select(Sample_ID,group, age,gender,Race)
summary(demo[demo$group=="Patient",])
sd(demo[demo$group=="Patient",]$age)
summary(demo[demo$group=="Healthy control",])
sd(demo[demo$group=="Healthy control",]$age)
a.div <- inner_join(a.div, demo)

### check the total reads count
total <- read.table(file="./2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/bacteria/120215JW515F-pr.fasta.otus.fa.OTU.txt", sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)

a <- colSums(total[, (7:42)])%>%data.frame()
colSums(a)
# subject
mean(a[(1:34),1])
sd(a[(1:34),1])
# water control
mean(a[(35:36),1])
sd(a[(35:36),1])

pt.a.div <- filter(a.div, group=="Patient")

pt.a.div$cd4_igg.cut <- cut(pt.a.div$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("low", "high"))




#figure 2

#L4
### pt vs. healthy
L4.per$Sample_ID <- rownames(L4.per)
L4 <- inner_join(L4.per, a.div[,c("Sample_ID", "group")])
L4$group <- revalue(L4$group, c("Patient"="HIV Patient"))
colMax <- function(data) sapply(data, max, na.rm=TRUE)
test <- cbind(sum=colSums(L4[,1:58]),mean=colMeans(L4[,1:58]),max=colMax(L4[,1:58]))%>%data.frame()

test <- data.frame(sapply(test, function(x) as.numeric(as.character(x))))
test$order <- colnames(L4[,1:58])
# classify any bacteria with a mean < 3% to "other bacteria"
colSums(test[,1:3])
test$order.other <-test$order
test$order.other[test$mean<=3] <- "bacteria_other"
test2 <- test[test$mean<=3,]
L4$bacteria_other <- rowSums(L4[,(colnames(L4)%in% test2$order)])
str(test2$order)
L4 <-L4[,-which(names(L4)%in%test2$order)]
str(L4)

order.level<- melt(L4, id.vars=c("group","Sample_ID"))


order.level$variable <- revalue(order.level$variable, c("o__actinomycetales"= "Actinomycetales","o__bacillales"= "Bacillales", "o__bacteroidales"="Bacteroidales","o__burkholderiales"="Burkholderiales","o__caulobacterales"="Caulobacterales","o__desulfuromonadales"="Desulfuromonadales","o__enterobacteriales"= "Enterobacteriales","o__flavobacteriales"="Flavobacteriales","o__hydrogenophilales"="Hydrogenophilales","o__lactobacillales"="Lactobacillales","o__myxococcales"="Myxococcales","o__oscillatoriales"="Oscillatoriales", "o__pseudomonadales"="Pseudomonadales","o__rhizobiales"="Rhizobiales","o__rhodobacterales"="Rhodobacterales", "o__sphingomonadales"="Sphingomonadales", "o__thermales"="Thermales","o__xanthomonadales"="Xanthomonadales", "bacteria_other"="Bacteria Other"))


order.color<-c("Actinomycetales" = "#FFBBFF", 
               "Bacillales" = "#8A2BE2", 
               "Bacteroidales"="blue",
               "Burkholderiales"="black",
               "Caulobacterales" = "#FFFF00", 
               "Desulfuromonadales" = "#458B00",
               "Enterobacteriales" = "#FF69B4", 
               "Flavobacteriales" = "#FF4040", 
               "Hydrogenophilales" = "#D2691E", 
               "Lactobacillales" = "#8A8A8A", 
               "Myxococcales" = "#53868B", 
               "Oscillatoriales" = "#000080", 
               "Pseudomonadales" = "#32CD32", 
               "Rhizobiales" = "#C0FF3E", 
               "Sphingomonadales" ="#0000FF" , 
               "Thermales" = "#B8860B", 
               "Bacteria Other" = "#FFE4C4",
               "Xanthomonadales" = "#A52A2A", 
               "Rhodobacterales" = "#00FFFF", 
               "gammaproteobacteria" = "#FFB90F", 
               "flavobacteriia" = "#C0FF3E")

order.level.2 <- aggregate(value~variable+group, data=order.level, FUN=mean) 

figure3E <- ggplot(order.level.2, aes(x = factor(group), y = value)) +scale_fill_manual(values=order.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3E


### only focus on pt group
pt.L4 <- filter(L4, group=="HIV Patient")
pt.L4 <- inner_join(pt.a.div[, c("Sample_ID","CD4.autoIgG..ng.ml.")],pt.L4)
pt.L4$cd4_igg.cut <- cut(pt.L4$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
order.level<- melt(pt.L4, id.vars=c("group","Sample_ID", "CD4.autoIgG..ng.ml.","cd4_igg.cut" ))
order.level$variable <- revalue(order.level$variable, c("o__actinomycetales"= "Actinomycetales","o__bacillales"= "Bacillales", "o__bacteroidales"="Bacteroidales","o__burkholderiales"="Burkholderiales","o__caulobacterales"="Caulobacterales","o__desulfuromonadales"="Desulfuromonadales","o__enterobacteriales"= "Enterobacteriales","o__flavobacteriales"="Flavobacteriales","o__hydrogenophilales"="Hydrogenophilales","o__lactobacillales"="Lactobacillales","o__myxococcales"="Myxococcales","o__oscillatoriales"="Oscillatoriales", "o__pseudomonadales"="Pseudomonadales","o__rhizobiales"="Rhizobiales","o__rhodobacterales"="Rhodobacterales", "o__sphingomonadales"="Sphingomonadales", "o__thermales"="Thermales","o__xanthomonadales"="Xanthomonadales", "bacteria_other"="Bacteria Other"))

str(order.level)
order.level.2 <- aggregate(value~variable+cd4_igg.cut, data=order.level, FUN=mean) 


figure3F <- ggplot(order.level.2, aes(x = factor(cd4_igg.cut), y = value)) +scale_fill_manual(values=order.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3F

library(gridExtra)
grid.arrange(figure3A,figure3B,figure3C,figure3D, ncol=2)

## composition of class level
L3.per$Sample_ID <- row.names(L3.per)
pt.L3 <- inner_join(L3.per, pt.L4[,c("Sample_ID","cd4_igg.cut")])
L3.low <- pt.L3 %>% filter(cd4_igg.cut==c("Anti-CD4 IgG baseline <= 50"))
L3.low <- colSums(L3.low[,(1:29)])%>% data.frame()
colSums(L3.low)
L3.low$frequency <- L3.low$./1300
L3.low$Sample_ID <- row.names(L3.low)
L3.low[c("c__gammaproteobacteria"),2]+L3.low[c("c__betaproteobacteria"),2]+L3.low[c("c__bacilli"),2]+L3.low[c("c__alphaproteobacteria"),2]

L3.high <- pt.L3 %>% filter(cd4_igg.cut==c("Anti-CD4 IgG baseline > 50"))
L3.high <- colSums(L3.high[,(1:29)])%>% data.frame()
colSums(L3.high)
L3.high$frequency <- L3.high$./900
L3.high$Sample_ID <- row.names(L3.high)
L3.high[c("c__gammaproteobacteria"),2]+L3.high[c("c__betaproteobacteria"),2]+L3.high[c("c__bacilli"),2]+L3.high[c("c__actinobacteria"),2]

### 1.1 data analysis on a-diversity and bacteria count.
# t test for simpson between groups
leveneTest(pt.a.div$simpson, pt.a.div$cd4_igg.cut)
t.test(pt.a.div$simpson~pt.a.div$cd4_igg.cut)
a.div$cd4_igg.cut <- cut(a.div$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
t.test(a.div$simpson~a.div$cd4_igg.cut)
a.div$groups <- a.div$cd4_igg.cut %>% as.character()

a.div$groups[a.div$group=="Healthy control"] <- "Healthy control"

a <- a.div%>%filter(!groups=="Anti-CD4 IgG baseline <= 50")
t.test(a$simpson~a$groups)

a <- a.div%>%filter(!groups=="Anti-CD4 IgG baseline > 50")
t.test(a$simpson~a$groups)
write.csv(a.div, file="a.div.test.csv")


# 1.1.3 Simpson diversity figure

a.div$cd4_igg.cut <- cut(a.div$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
a.div$cd4_igg.cut <- as.character(a.div$cd4_igg.cut)
a.div$cd4_igg.cut[a.div$group=="Healthy control"] <- "Healthy Control"
ggplot(a.div, aes(cd4_igg.cut, simpson))+geom_boxplot(aes(fill=cd4_igg.cut))+geom_point()+ labs(fill = "Key",x="", y="Simpson Diversity Index")+theme_bw()

### 1.2 cd4_igg and beta diversity
otu.per$Sample_ID <- row.names(otu.per)
pt.otu.all<- inner_join(otu.per, pt.a.div)
rownames(pt.otu.all) <- pt.otu.all$Sample_ID
## 1.2.3 permanova bray-curtis 
#patient only
adonis(pt.otu.all[,(1:2413)] ~pt.otu.all$CD4.autoIgG..ng.ml.+pt.otu.all$LPS.Plasma.pg.ml.+pt.otu.all$ART.yr.+pt.otu.all$CD4.ul+pt.otu.all$age+pt.otu.all$X16S.Watercont +pt.otu.all$gender, na.rm=TRUE, permutations=999, method="bray")#sig
adonis(pt.otu.all[,(1:2413)] ~ pt.otu.all$cd4_igg.cut+pt.otu.all$ART.yr.+pt.otu.all$CD4.ul, na.rm=TRUE, permutations=999, method="bray")#not sig

# patient and healthy
otu.clinical <- inner_join(otu, a.div)
adonis(otu.clinical[,(1:2413)] ~otu.clinical$CD4.autoIgG..ng.ml.+otu.clinical$LPS.Plasma.pg.ml.+otu.clinical$ART.yr.+otu.clinical$CD4.ul+otu.clinical$age+otu.clinical$X16S.Watercont+otu.clinical$gender, na.rm=TRUE, permutations=999, method="bray")#not sig
adonis(otu.clinical[,(1:2413)] ~otu.clinical$cd4_igg.cut+otu.clinical$LPS.Plasma.pg.ml.+otu.clinical$ART.yr.+otu.clinical$CD4.ul, na.rm=TRUE, permutations=999, method="bray")

# groups and negative control
otu.raw$Sample_ID <- rownames(otu.raw)
otu.raw.clinical <- full_join(otu.raw, a.div[,c("Sample_ID","group")])
otu.raw.clinical$group <- as.character(otu.raw.clinical$group)
otu.raw.clinical$group[is.na(otu.raw.clinical$group)] <- "negcontrol"
otu.raw.clinical$group
adonis(otu.raw.clinical[,(1:2413)] ~otu.raw.clinical$group, na.rm=TRUE, permutations=999, method="bray")
otu.raw.clinical$group2 <- otu.raw.clinical$group
otu.raw.clinical$group2[!otu.raw.clinical$group=="negcontrol"] <- "subjects"
adonis(otu.raw.clinical[,(1:2413)] ~otu.raw.clinical$group2, na.rm=TRUE, permutations=999, method="bray")


## 1.2.4 permanova jaccard
#patient only
adonis(pt.otu.all[,(1:2413)] ~pt.otu.all$CD4.autoIgG..ng.ml.+pt.otu.all$LPS.Plasma.pg.ml.+pt.otu.all$ART.yr.+pt.otu.all$CD4.ul,  na.rm=TRUE, permutations=999, method="jaccard")#sig

# patient and healthy
adonis(otu.clinical[,(1:2413)] ~otu.clinical$CD4.autoIgG..ng.ml.+otu.clinical$LPS.Plasma.pg.ml.+otu.clinical$ART.yr.+otu.clinical$CD4.ul, na.rm=TRUE, permutations=999, method="jaccard")#not sig

## 1.2.6 nms with three groups, envfit
otu.clinical <- inner_join(otu, a.div[, c("Sample_ID","group","CD4.autoIgG..ng.ml.","age","gender", "ART.yr.", "CD4.ul","LPS.Plasma.pg.ml.", "X16S.Watercont")])

rownames(otu.clinical) <- otu.clinical$Sample_ID
otu.clinical$cd4_igg.cut <-cut(otu.clinical$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
otu.clinical$cd4_igg.cut <- as.character(otu.clinical$cd4_igg.cut)
otu.clinical$cd4_igg.cut[otu.clinical$group=="Healthy control"] <- "Healthy Control"

bc.nms <- metaMDS(otu.clinical[,1:1974], k=2, trymin=50, trymax=500, wascores=F)
scrs <- as.data.frame(scores(bc.nms, display ='sites'))
scrs <- cbind(scrs, `anti-CD4 IgG` = otu.clinical$CD4.autoIgG..ng.ml., `Simpson Diversity`=otu.clinical$simpson, Gender = otu.clinical$gender, `Duration of ART`=otu.clinical$ART.yr., Age = otu.clinical$age, "Plasma LPS"=otu.clinical$LPS.Plasma.pg.ml., "16s rDNA count" = otu.clinical$X16S.Watercont,CD4=otu.clinical$CD4.ul,group=otu.clinical$cd4_igg.cut)

scrs$Gender <- revalue(scrs$Gender, c("M"="1", "F"="0"))
scrs$Gender <- as.numeric(scrs$Gender)

#spec <- envfit(bc.nms,otu.clinical[,1:2490] , perm=999)
spec2<- envfit(bc.nms, scrs[,(4:10)], perm=999)
scores(spec2, "vectors")
scrs$groupcolor <- revalue(scrs$group,c("Anti-CD4 IgG baseline <= 50"="red", "Anti-CD4 IgG baseline > 50"="purple", "Healthy Control"="green"))
plot(scrs[,1:2], col=scrs$groupcolor, pch=16)
ordiellipse(scrs[,1:2], scrs$group,draw = "line", display="sites", kind = "sd", col =c("red","purple","green"))

plot(spec2,cex=0.7, p.max = 0.05)
legend("topleft", legend=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50","Healthy Control"),col=c("red","purple", "green"), lty=1, cex=0.8)

## 1.2.7 nms with two groups, envfit
otu.clinical <- inner_join(otu, pt.a.div[, c("Sample_ID","CD4.autoIgG..ng.ml.","age","gender","LPS.Plasma.pg.ml.","X16S.Watercont", "ART.yr.", "CD4.ul")])
rownames(otu.clinical) <- otu.clinical$Sample_ID
otu.clinical$cd4_igg.cut <-cut(otu.clinical$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))


bc.nms <- metaMDS(otu.clinical[,1:1974], k=2, trymin=50, trymax=500, wascores=F)
scrs <- as.data.frame(scores(bc.nms, display ='sites'))
scrs <- cbind(scrs, `Anti-CD4 IgG` = otu.clinical$CD4.autoIgG..ng.ml., `Simpson Diversity`=otu.clinical$simpson, "Gender" = otu.clinical$gender, Age = otu.clinical$age, "Plasma LPS"=otu.clinical$LPS.Plasma.pg.ml., "16s rDNA count" = otu.clinical$X16S.Watercont,`Duration of ART`=otu.clinical$ART.yr., `CD4 count`=otu.clinical$CD4.ul,group=otu.clinical$cd4_igg.cut )

scrs$`Gender` <- revalue(scrs$`Gender`, c("M"="1", "F"="0"))
scrs$`Gender`<- as.numeric(scrs$`Gender`)

#spec <- envfit(bc.nms,otu.clinical[,1:2490] , perm=999)
spec2<- envfit(bc.nms, scrs[,(4:10)], perm=999)
scores(spec2, "vectors")
scrs$groupcolor <- revalue(scrs$group,c("Anti-CD4 IgG baseline <= 50"="red", "Anti-CD4 IgG baseline > 50"="green"))
plot(scrs[,1:2], col=as.character(scrs$groupcolor), pch=16)
ordiellipse(scrs[,1:2], scrs$group,draw = "polygon", display="sites", kind="sd", col =c("red","green"))
#ordiellipse(scrs[,1:2], scrs$group,draw = "line", display="sites", kind = "sd", col =c("red","green"))
plot(spec2,cex=0.9)
#legend("topleft", legend=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"),col=c("red","green"), lty=1, cex=0.8)

# plot the phyla level
L2.per$Sample_ID <- rownames(L2.per)
L2.pt <- L2.per[rownames(L2.per)%in%rownames(otu.clinical),]
L2.pt <- rename(L2.pt, c("p__acidobacteria"="Acidobacteria","p__actinobacteria"="Actinobacteria","p__armatimonadetes"="Armatimonadetes" , "p__bacteroidetes"="Bacteroidetes", "p__candidatus saccharibacteria"="Candidatus Saccharibacteria", "p__chloroflexi"="Chloroflexi","p__cyanobacteria"="Cyanobacteria", "p__deinococcus_thermus"="Deinococcust Thermus" ,"p__firmicutes"="Firmicutes" , "p__gemmatimonadetes"="Gemmatimonadetes" ,"p__planctomycetes"="Planctomycetes","p__proteobacteria"="Proteobacteria","p__tenericutes"="Tenericutes", "p__verrucomicrobia"="Verrucomicrobia" ))
spec<- envfit(bc.nms, L2.pt[,1:14], perm=999)
head(scores(spec, "vectors"))
scrs$groupcolor <- revalue(scrs$group,c("Anti-CD4 IgG baseline <= 50"="red", "Anti-CD4 IgG baseline > 50"="green"))
plot(scrs[,1:2], col=as.character(scrs$groupcolor), pch=16)
#ordiellipse(scrs[,1:2], scrs$group,draw = "polygon", display="sites", kind = "se", col =c("red","green"))
ordiellipse(scrs[,1:2], scrs$group,draw = "polygon", display="sites", kind = "sd", col =c("red","green"))
plot(spec, cex=0.8, p.max=0.05)
#legend("topleft", legend=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"),col=c("red","green"), lty=1, cex=0.8)

### 1.3 indicator species
## 1.3.1 cd4_igg_d0 at genus level
# cut cd4_igg.d0 by 50 
L6$Sample_ID <- rownames(L6)
pt.L6 <- inner_join(L6, pt.a.div[,c("Sample_ID", "CD4.autoIgG..ng.ml.")])
pt.L6$cd4_igg.cut <- cut(pt.L6$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
indic <- multipatt(pt.L6[,(1:223)], pt.L6$cd4_igg.cut, control = how(nperm=999))
write.csv(file="indicator.species.csv",indic$sign%>%
            tibble::rownames_to_column(var = "genus")%>%
            mutate(p.fdr = round(p.adjust(p.value, "fdr"),3))%>%
            right_join(taxa, by = "genus")%>%
            filter(p.value < 0.05) %>%
            arrange(index))

summary(indic, alpha=1)
signassoc(pt.L6[,(1:223)], cluster=pt.L6$cd4_igg.cut,  alternative = "two.sided",control = how(nperm=999))


pt.otu.all <- inner_join(otu, pt.a.div[,c("Sample_ID", "CD4.autoIgG..ng.ml.")])
pt.otu.all$cd4_igg.cut <- cut(pt.otu.all$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
indic <- multipatt(pt.otu.all[,(1:2413)], pt.otu.all$cd4_igg.cut, control = how(nperm=999))
write.csv(file="indicator.species.csv",indic$sign%>%
            tibble::rownames_to_column(var = "OTU")%>%
            mutate(p.fdr = round(p.adjust(p.value, "fdr"),3))%>%
            right_join(taxa, by = "OTU")%>%
            filter(p.value < 0.05) %>%
            arrange(index))


write.csv(pt.L6, file="test.csv")

indicL6 <-  multipatt(pt.L6[,(1:223)], pt.L6$cd4_igg.cut, control = how(nperm=999))
summary(indicL6)


#L2
### pt vs. healthy
L2.per$Sample_ID <- rownames(L2.per)
L2 <- inner_join(L2.per, a.div[,c("Sample_ID", "group")])
L2$group <- revalue(L2$group, c("Patient"="HIV Patient"))
colMax <- function(data) sapply(data, max, na.rm=TRUE)
test <- cbind(sum=colSums(L2[,1:14]),mean=colMeans(L2[,1:14]),max=colMax(L2[,1:14]))%>%data.frame()

test <- data.frame(sapply(test, function(x) as.numeric(as.character(x))))
test$phyla <- colnames(L2[,1:14])
# classify any bacteria with a mean < 0.1 to "other bacteria"

test$phyla.other[test$mean<0.1] <- "bacteria_other"
test$phyla.other[test$mean>=0.1] <-test$phyla
test2 <- test[test$mean<=0.1,]
L2$bacteria_other <- rowSums(L2[,(colnames(L2)%in% test2$phyla)])
str(test2$phyla)
L2 <-L2[,-which(names(L2)%in%test2$phyla)]
str(L2)

phyla.level<- melt(L2, id.vars=c("group","Sample_ID"))


phyla.level$variable <- revalue(phyla.level$variable, c("p__acidobacteria"= "Acidobacteria","p__actinobacteria"= "Actinobacteria", "p__armatimonadetes"="Armatimonadetes","p__bacteroidetes"="Bacteroidetes","p__chloroflexi"="Chloroflexi","p__cyanobacteria"="Cyanobacteria","p__deinococcus_thermus"= "Deinococcus-Thermus","p__firmicutes"="Firmicutes","p__planctomycetes"="Planctomycetes","p__proteobacteria"="Proteobacteria","p__tenericutes"="Tenericutes", "bacteria_other"="Bacteria Other"))


phyla.color<-c("Acidobacteria" = "#FFBBFF", 
               "Actinobacteria" = "#8A2BE2", 
               "Armatimonadetes"="blue",
               "Bacteroidetes"="black",
               "Chloroflexi" = "#FFFF00", 
               "Cyanobacteria" = "#458B00",
               "Deinococcus-Thermus" = "#FF69B4", 
               "Firmicutes" = "#FF4040", 
               "Planctomycetes" = "#D2691E", 
               "Proteobacteria" = "#8A8A8A", 
               "Tenericutes" = "#53868B", 
               "B" = "#000080", 
               "Gammaproteobacteria" = "#32CD32", 
               "Mollicutes" = "#C0FF3E", 
               "Phycisphaerae" ="#0000FF" , 
               "Planctomycetia" = "#B8860B", 
               "Bacteria Other" = "#FFE4C4",
               "cyanobacteri" = "#A52A2A", 
               "clostridia" = "#00FFFF", 
               "gammaproteobacteria" = "#FFB90F", 
               "flavobacteriia" = "#C0FF3E")

phyla.level.2 <- aggregate(value~variable+group, data=phyla.level, FUN=mean) 

figure3A <- ggplot(phyla.level.2, aes(x = factor(group), y = value)) +scale_fill_manual(values=phyla.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3A


### only focus on pt group
pt.L2 <- filter(L2, group=="HIV Patient")
pt.L2 <- inner_join(pt.a.div[, c("Sample_ID","CD4.autoIgG..ng.ml.")],pt.L2)
pt.L2$cd4_igg.cut <- cut(pt.L2$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50"))
phyla.level<- melt(pt.L2, id.vars=c("group","Sample_ID", "CD4.autoIgG..ng.ml.","cd4_igg.cut" ))
phyla.level$variable <- revalue(phyla.level$variable, c("p__acidobacteria"= "Acidobacteria","p__actinobacteria"= "Actinobacteria", "p__armatimonadetes"="Armatimonadetes","p__bacteroidetes"="Bacteroidetes","p__chloroflexi"="Chloroflexi","p__cyanobacteria"="Cyanobacteria","p__deinococcus_thermus"= "Deinococcus-Thermus","p__firmicutes"="Firmicutes","p__planctomycetes"="Planctomycetes","p__proteobacteria"="Proteobacteria","p__tenericutes"="Tenericutes", "bacteria_other"="Bacteria Other"))

str(phyla.level)
phyla.level.2 <- aggregate(value~variable+cd4_igg.cut, data=phyla.level, FUN=mean) 


figure3B <- ggplot(phyla.level.2, aes(x = factor(cd4_igg.cut), y = value)) +scale_fill_manual(values=phyla.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3B

### heatmap
library(heatmap3)






L4.per$Sample_ID <- rownames(L4.per)
L4 <- inner_join(L4.per, a.div[,c("Sample_ID", "group")])
L4$group <- revalue(L4$group, c("Patient"="HIV Patient"))
colMax <- function(data) sapply(data, max, na.rm=TRUE)
test <- cbind(sum=colSums(L4[,1:58]),mean=colMeans(L4[,1:58]),max=colMax(L4[,1:58]))%>%data.frame()



#L3
### pt vs. healthy
L3.per$Sample_ID <- rownames(L3.per)
L3 <- inner_join(L3.per, a.div[,c("Sample_ID", "group")])
L3$group <- revalue(L3$group, c("Patient"="HIV+ subject (N=22)", "Healthy control"="Healthy control (N=12)"))
colMax <- function(data) sapply(data, max, na.rm=TRUE)
test <- cbind(colSums(L3[,1:29]),colMeans(L3[,1:29]),colMax(L3[,1:29]))%>% data.frame()

# classify any bacteria with a mean < 2 to "other bacteria"
test <- cbind(sum=colSums(L3[,1:29]),mean=colMeans(L3[,1:29]),max=colMax(L3[,1:29]))%>%data.frame()
test <- data.frame(sapply(test, function(x) as.numeric(as.character(x))))
test$class <- colnames(L3[,1:29])
# classify any bacteria 
colSums(test[,1:3])
test$class.other <-test$class
test$class.other[test$mean<=2] <- "bacteria_other"
test2 <- test[test$mean<=2,]
L3$bacteria_other <- rowSums(L3[,(colnames(L3)%in% test2$class)])
str(test2$class)
L3 <-L3[,-which(names(L3)%in%test2$class)]
str(L3)







# classify any bacteria with a mean < 0.1 to "other bacteria"

#L3$Bacteria_other <-  L3$c__caldilineae+ L3$c__candidatus.saccharibacteria+ L3$c__chloroflexi+ L3$c__chloroflexia+ L3$c__dehalococcoidia+ L3$c__erysipelotrichia +L3$c__gemmatimonadetes +L3$c__negativicutes+ L3$c__sphingobacteriia +L3$c__thermoleophilia +L3$c__thermomicrobia+ L3$c__verrucomicrobiae

#L3$c__acidobacteria <- L3$c__acidobacteria+L3$c__acidobacteriia

#L3 <- L3[,c ("c__acidobacteria","c__actinobacteria", "c__alphaproteobacteria","c__bacilli","c__bacteroidia","c__betaproteobacteria","c__chthonomonadetes","c__clostridia","c__cyanobacteria","c__deinococci","c__deltaproteobacteria","c__flavobacteriia","c__gammaproteobacteria", "c__mollicutes", "c__phycisphaerae", "c__planctomycetia", "Bacteria_other","Sample_ID","group")]





class.color<-c("Acidobacteria" = "blue", 
               "Actinobacteria" = "darkmagenta", 
               "Alphaproteobacteria"="darkorange",
               "Bacilli"="gold",
               "Bacteroidia" = "pink", 
               "Betaproteobacteria" = "#458B00",
               "Chthonomonadetes" = "#FF69B4", 
               "Clostridia" = "#FF4040", 
               "Cyanobacteria" = "darkolivegreen1", 
               "Deinococci" = "#8A8A8A", 
               "Deltaproteobacteria" = "#53868B", 
               "Flavobacteriia" = "#000080", 
               "Gammaproteobacteria" = "dodgerblue", 
               "Mollicutes" = "#B8860B", 
               "Phycisphaerae" ="#0000FF" , 
               "Planctomycetia" = "#C0FF3E", 
               "Bacteria Other" = "rosybrown4",
               "cyanobacteri" = "#A52A2A", 
               "clostridia" = "#00FFFF", 
               "flavobacteriia" = "#C0FF3E")
class.level<- melt(L3, id.vars=c("group","Sample_ID"))


class.level$variable <- revalue(class.level$variable, c("c__acidobacteria"= "Acidobacteria","c__actinobacteria"= "Actinobacteria", "c__alphaproteobacteria"="Alphaproteobacteria","c__bacilli"="Bacilli","c__bacteroidia"="Bacteroidia","c__betaproteobacteria"="Betaproteobacteria","c__chthonomonadetes"= "Chthonomonadetes","c__clostridia"="Clostridia","c__cyanobacteria"="Cyanobacteria","c__deinococci"="Deinococci","c__deltaproteobacteria"="Deltaproteobacteria","c__flavobacteriia"="Flavobacteriia","c__gammaproteobacteria"="Gammaproteobacteria", "c__mollicutes"="Mollicutes", "c__phycisphaerae"="Phycisphaerae", "c__planctomycetia"="Planctomycetia", "bacteria_other"="Bacteria Other"))


class.level.2 <- aggregate(value~variable+group, data=class.level, FUN=mean) 

figure3C <- ggplot(class.level.2, aes(x = factor(group), y = value)) +scale_fill_manual(values=class.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3C


### only focus on pt group
pt.L3 <- filter(L3, group=="HIV+ subject (N=22)")
pt.L3 <- inner_join(pt.a.div[, c("Sample_ID","CD4.autoIgG..ng.ml.")],pt.L3)
pt.L3$cd4_igg.cut <- cut(pt.L3$CD4.autoIgG..ng.ml., breaks=c(-Inf, 50, Inf), labels=c("HIV+/aCD4<50 (N=13)", "HIV+/aCD4>50 (N=9)"))
class.level<- melt(pt.L3, id.vars=c("group","Sample_ID", "CD4.autoIgG..ng.ml.","cd4_igg.cut" ))
class.level$variable <- revalue(class.level$variable, c("c__acidobacteria"= "Acidobacteria","c__actinobacteria"= "Actinobacteria", "c__alphaproteobacteria"="Alphaproteobacteria","c__bacilli"="Bacilli","c__bacteroidia"="Bacteroidia","c__betaproteobacteria"="Betaproteobacteria","c__chthonomonadetes"= "Chthonomonadetes","c__clostridia"="Clostridia","c__cyanobacteria"="Cyanobacteria","c__deinococci"="Deinococci","c__deltaproteobacteria"="Deltaproteobacteria","c__flavobacteriia"="Flavobacteriia","c__gammaproteobacteria"="Gammaproteobacteria", "c__mollicutes"="Mollicutes", "c__phycisphaerae"="Phycisphaerae", "c__planctomycetia"="Planctomycetia", "bacteria_other"="Bacteria Other"))
str(class.level)
class.level.2 <- aggregate(value~variable+cd4_igg.cut, data=class.level, FUN=mean) 


figure3D <- ggplot(class.level.2, aes(x = factor(cd4_igg.cut), y = value)) +scale_fill_manual(values=class.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figure3D



## heatmap
library(gplots)
library(RColorBrewer)

# version 1: calculate and remove the sum otu < 5, then group the otu based on genus
genus <- inner_join(taxa[,c("OTU","genus")], otu2)
genus.per <- cbind(genus[,1:2], prop.table(as.matrix(genus[,-(1:2)]), margin = 2))%>% data.frame()
genus.per[,(3:36)] <- genus.per[,(3:36)]*100
colSums(genus.per[,(3:36)])
genus.per$max <- apply(genus.per[,(3:36)], 1, max)

genus.per.2 <- genus.per[genus.per$max>=5,]
str(genus.per.2)
genus.heatmap <- genus.per.2[,(2:36)] %>% group_by(genus)%>% summarise_each(funs(sum))

row.names(genus.heatmap) <- genus.heatmap$genus
genus.heatmap <- t(genus.heatmap) %>% data.frame() 
genus.heatmap <- genus.heatmap[-1,]
names(genus.heatmap) <- gsub("g__", "", names(genus.heatmap), fixed = TRUE) #remove "g__"
names(genus.heatmap) <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",names(genus.heatmap), perl = TRUE) #repalce the first letter to upper case
str(genus.heatmap)
scaleyellowred <- colorRampPalette(c("white", "blue"), space="rgb")(100)

genus.heatmap2 <- data.frame(sapply(genus.heatmap, function(x) as.numeric(as.character(x))))
row.names(genus.heatmap2) <- rownames(genus.heatmap)
str(genus.heatmap2)

heatmap(as.matrix(genus.heatmap2), Rowv=NA, Colv=NA, col=scaleyellowred, margins = c(10,2))

# calculate the bray-curtis dissimilarity matrix on the full dataset
data.dist <- vegdist(genus.heatmap2, method="bray")
row.clus <- hclust(data.dist,"aver")
heatmap(as.matrix(genus.heatmap2), Rowv = as.dendrogram(row.clus), Colv = NA, col=scaleyellowred, margins= c(10, 2))

data.dist.g <- vegdist(t(genus.heatmap2), method="bray")
col.clus <- hclust(data.dist.g, "aver")
heatmap(as.matrix(genus.heatmap2), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col=scaleyellowred, margins= c(10, 2))

genus.heatmap2$Sample_ID <- row.names(genus.heatmap2)
genus.heatmap <- inner_join(genus.heatmap2, a.div[,c("Sample_ID", "cd4_igg.cut")])
genus.heatmap$cd4.color <- revalue(genus.heatmap$cd4_igg.cut, c("Anti-CD4 IgG baseline <= 50"= "red", "Anti-CD4 IgG baseline > 50"="green","Healthy Control"="blue"))
row.names(genus.heatmap)<-genus.heatmap$Sample_ID 
cd.color <- genus.heatmap[,62]
genus.heatmap <- genus.heatmap[,(1:59)]


namesotugenus <- colnames(genus.heatmap)%>% data.frame()
namesotugenus <- rename(namesotugenus, c("."="genus"))
taxa2 <- taxa[,c("genus","phyla")]
taxa2 <- taxa2[!duplicated(taxa2),]
taxa2$genus <- gsub("g__", "", taxa2$genus, fixed = TRUE) #remove "g__"
taxa2$genus <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",taxa2$genus, perl = TRUE) #repalce the first letter to upper case
taxa2$phyla <- gsub("p__", "", taxa2$phyla, fixed = TRUE) #remove "g__"
taxa2$phyla <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",taxa2$phyla, perl = TRUE) #repalce the first letter to upper case
namesotugenus <- inner_join(namesotugenus,taxa2)
table(namesotugenus$phyla)
namesotugenus$phyla.color <- revalue(namesotugenus$phyla, c("Actinobacteria"= "blue", "Bacteroidetes"="red","Cyanobacteria"="green", "Deinococcus_thermus"="yellow", "Firmicutes"="purple", "Proteobacteria"="orange"))


heatmap.2(as.matrix(genus.heatmap),Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = cd.color,ColSideColors=namesotugenus$phyla.color,trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", margins=c(8,5))

showLegend(legend = c("Anti-CD4 IgG baseline <= 50", "Anti-CD4 IgG baseline > 50","Healthy Control"),fill=c( "red","green","blue"), lty=0 )

showLegend(legend = c("Actinobacteria", "Bacteroidetes","Cyanobacteria", "Deinococcus_thermus", "Firmicutes", "Proteobacteria"),fill=c( "blue","red","green", "yellow", "purple", "orange"), lty=0 )

test <- inner_join(genus.heatmap2, a.div[, c("Sample_ID", "cd4_igg.cut")])
write.csv(test, file="test.csv")

#version 2, load all otu with genus name and remove any with max value < 5%
genus <- inner_join(taxa[,c("OTU","genus")],otu2)
genus.per <- cbind(genus[,1:2], prop.table(as.matrix(genus[,-(1:2)]), margin = 2))%>% data.frame()
genus.per[,(3:36)] <- genus.per[,(3:36)]*100
colSums(genus.per[,(3:36)])
genus.heatmap <- genus.per
row.names(genus.heatmap) <- genus.heatmap$OTU
genus.heatmap <- t(genus.heatmap) %>% data.frame() 
names(genus.heatmap) <- lapply(genus.heatmap[2, ], as.character)

names(genus.heatmap) <- gsub("g__", "", names(genus.heatmap), fixed = TRUE) #remove "g__"
names(genus.heatmap) <- gsub("^(\\w)(\\w+)", "\\U\\1\\L\\2",names(genus.heatmap), perl = TRUE) #repalce the first letter to upper case
namesotugenus <- genus.heatmap[(1:2),]
genus.heatmap <- genus.heatmap[-(1:2),]

scaleyellowred <- colorRampPalette(c("white", "blue"), space="rgb")(100)

genus.heatmap2 <- data.frame(sapply(genus.heatmap, function(x) as.numeric(as.character(x))))
row.names(genus.heatmap2) <- rownames(genus.heatmap)
str(genus.heatmap2)
heatmap(as.matrix(genus.heatmap2), Rowv=NA, Colv=NA, col=scaleyellowred)
maxab <- apply(genus.heatmap2, 2, max)
head(maxab)
rowSums(genus.heatmap2)
n1 <- names(which(maxab<5)) #anything less than 5 percent.
genus.heatmap <- genus.heatmap2[, -which(names(genus.heatmap2) %in% n1)]
heatmap(as.matrix(genus.heatmap), Rowv=NA, Colv=NA, col=scaleyellowred, margins = c(10,2))

# calculate the bray-curtis dissimilarity matrix on the full dataset
data.dist <- vegdist(genus.heatmap, method="bray")
row.clus <- hclust(data.dist,"aver")
heatmap(as.matrix(genus.heatmap), Rowv = as.dendrogram(row.clus), Colv = NA, col=scaleyellowred, margins= c(10, 3))

data.dist.g <- vegdist(t(genus.heatmap), method="bray")
col.clus <- hclust(data.dist.g, "aver")
heatmap(as.matrix(genus.heatmap), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col=scaleyellowred, margins= c(10, 3))

genus.heatmap$Sample_ID <- row.names(genus.heatmap)
genus.heatmap <- inner_join(genus.heatmap, a.div[,c("Sample_ID", "cd4_igg.cut")])
genus.heatmap$cd4.color <- revalue(genus.heatmap$cd4_igg.cut, c("Anti-CD4 IgG baseline <= 50"= "red", "Anti-CD4 IgG baseline > 50"="green","Healthy Control"="blue"))
row.names(genus.heatmap)<-genus.heatmap$Sample_ID 
cd.color <- genus.heatmap[,89]
genus.heatmap <- genus.heatmap[,(1:86)]
# needs to figure out
#namesotugenus2 <- namesotugenus[,colnames(namesotugenus) %in% colnames(genus.heatmap)] 
namesotugenus2 <- t(namesotugenus)%>%data.frame()
namesotugenus2$names2 <- colnames(namesotugenus)
namesotugenus2 <- inner_join(namesotugenus2[,(2:3)],taxa[,c("genus","phyla")])

namesotugenus2 <- namesotugenus2[namesotugenus2$names2 %in% colnames(genus.heatmap),]


heatmap.2(as.matrix(genus.heatmap),Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = cd.color, trace = "none", density.info = "none", xlab = "genera", ylab = "Samples")



#supplemental table1
sup1 <- read.table(file="./2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/bacteria/120215JW515F-pr.fasta.otus.fa.OTU.percentages.txt",sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)
sup1 <- sup1[,(2:37)]
colSums (sup1, na.rm = FALSE, dims = 1) 
sup1 <- as.data.frame(t(sup1)) 
sup1$Sample_ID <- rownames(sup1)
reads <- rowSums(otu.raw)%>% data.frame() 
reads$reads <- reads$. 
reads$Sample_ID <- rownames(reads) %>% as.character()
reads <- full_join(reads, a.div[, c("Sample_ID" ,"groups")])
sup1 <- inner_join(reads,sup1)
write.csv(sup1, file="suplemental 1.csv")


#comparison of OTUs and absolute abundance in negative and subjects
#abundance
L3 <- read.table(file="./2015_12_15_120215JW515F_Analysis_Pipeline/analysisfiles/bacteria/120215JW515F-pr.fasta.otus.fa.classs.txt",sep="\t", header=T,row.names = 1, stringsAsFactors = F, as.is=T)
L3 <- L3[,(1:36)]
L3 <- as.data.frame(t(L3)) 
L3$Sample_ID <- row.names(L3)
L3 <- full_join(L3,a.div[,c("Sample_ID", "group")] )
str(L3)
L3$group <- as.character(L3$group)
L3$group[is.na(L3$group)] <- "negative"
L3 <- L3 %>% group_by(group)%>% summarise_each(funs(mean))
d <- t(L3[,2:30])%>% data.frame()
colnames(d) <- L3$group
d$class <- rownames(d)
d <- d %>% dplyr::rename(`Healthy control abudance`=`Healthy control`, `negative abundance`=negative, `Patient abundance`=Patient)

# patient count
otu.raw2 <- otu.raw
otu.raw2$Sample_ID <- rownames(otu.raw2)
otu.raw2 <- full_join(otu.raw2,a.div[,c("Sample_ID", "group")] )
patient.raw <- otu.raw2 %>% filter(group=="Patient")
patientmean <- colMeans(patient.raw[,(1:2413)])
patientmean <- as.data.frame(patientmean)
patientmean$OTU <- row.names(patientmean)
patient <- inner_join(patientmean,taxa)%>% filter(patientmean>0)

# healthy count
healthy.raw <- otu.raw2 %>% filter(group=="Healthy control")
healthymean <- colMeans(healthy.raw[,(1:2413)])
healthymean <- as.data.frame(healthymean)
healthymean$OTU <- row.names(healthymean)
healthy <- inner_join(healthymean,taxa)%>% filter(healthymean>0)

#total 
totalmean <- colMeans(otu.raw2[,(1:2413)])
totalmean <- as.data.frame(totalmean)
totalmean$OTU <- row.names(totalmean)
total <- inner_join(totalmean,taxa)%>% filter(totalmean>0)



a <- table(removed$class) %>% data.frame()%>%dplyr::rename(`negative control number of OTU`=Freq)

b <- table(patient$class) %>% data.frame()%>%dplyr::rename(`HIV patient number of OTU` =Freq)
c <- table(healthy$class) %>% data.frame()%>%dplyr::rename(`Healthy control number of OTU` =Freq)
g <- table(total$class) %>% data.frame()%>%dplyr::rename(`total number of OTU` =Freq)

e <- full_join(a,b)
e <- full_join(e,c)
e <- full_join(e,g)%>%dplyr::rename(class=Var1)
e$class <- revalue(e$class,c("c__acidobacteria"="acidobacteria","c__acidobacteriia"= "acidobacteriia"  ,  "c__actinobacteria"=  "actinobacteria" , "c__alphaproteobacteria"= "alphaproteobacteria", "c__bacilli"=  "bacilli", "c__bacteroidia"= "bacteroidia", "c__betaproteobacteria"="betaproteobacteria" , "c__caldilineae"="caldilineae", "c__candidatus saccharibacteria"= "candidatus saccharibacteria",   "c__chloroflexi"="chloroflexi" ,"c__chloroflexia"="chloroflexia", "c__chthonomonadetes"="chthonomonadetes", "c__clostridia"="clostridia" , "c__cyanobacteria"="cyanobacteria" , "c__dehalococcoidia"="dehalococcoidia" , "c__deinococci"="deinococci" , "c__deltaproteobacteria"="deltaproteobacteria", "c__erysipelotrichia"="erysipelotrichia",  "c__flavobacteriia"="flavobacteriia"  , "c__gammaproteobacteria"="gammaproteobacteria","c__gemmatimonadetes"="gemmatimonadetes",  "c__mollicutes"="mollicutes" , "c__negativicutes"= "negativicutes" ,  "c__phycisphaerae"="phycisphaerae" , "c__planctomycetia"="planctomycetia",  "c__sphingobacteriia"="sphingobacteriia" ,  "c__thermoleophilia"="thermoleophilia",  "c__thermomicrobia"="thermomicrobia" ,  "c__verrucomicrobiae"= "verrucomicrobiae"))
e <- e %>%mutate(negpercent=(`negative control number of OTU`/`total number of OTU`),hivpercent=(`HIV patient number of OTU`/`total number of OTU`), healthypercent= (`Healthy control number of OTU`/`total number of OTU`))

f <- full_join(d,e)
write.csv(f, file="class level comparision.csv")

# diversity of the otu
h <- removed %>% filter(OTU %in% healthy$OTU) 
413/439*100 # percent of OTU in negative control that were found in healthy subject.

h <- removed %>% filter(OTU %in% patient$OTU) 
431/439*100 # percent of OTU in negative control that were found in HIV patient.

h <- healthy %>% filter(OTU %in% removed$OTU) 
413/1662*100 # percent of OTU in healthy subject that were found in negative control.

h <- healthy %>% filter(OTU %in% patient$OTU) 
1476/1662*100 # percent of OTU in healthy subject that were found in negative control.

h <- patient %>% filter(OTU %in% removed$OTU) 
431/2222*100 # percent of OTU in HIV patient that were found in negative control.

h <- patient %>% filter(OTU %in% healthy$OTU) 
1476/2222*100 # percent of OTU in HIV patient that were found in healthy control.

