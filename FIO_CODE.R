# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Wed Apr 12 14:40:13 EDT 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE43861", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15207", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "0011"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="none" ,lfc = 0.584,p.value = 0.05 , number=11000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","GI","Gene.Symbol","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE43861", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15207", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "0011"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("FIO","Treatment")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE43861", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")


library(EMA)
library(cluster)
library(devtools)
library(ggfortify)
library(ggplot2)
library(som)
library(gplots)
#setwd("C:\\Users\\ikhatri\\Desktop\\Panc-New\\DE\\24279\\24279-new")
data=read.table("C:\\Users\\alokum\\Desktop\\FIO\\ex_gene_ver2.txt", header=TRUE, sep="\t")
mv.genes<-genes.selection(data[,-c(1)], thres.num=200)
selected=data[mv.genes,]
sel=normalize(selected)
heatmap.2(as.matrix(sel),ColSideColors=c(rep("green",2),rep("blue",2)),col=redgreen(75),scale="none",key=TRUE,symkey=FALSE,density.info="none",trace="none",cexRow=0.5,dendrogram="row",Colv = F)


####plot selected genes
var=c("HIST1H3G"	,"TNFAIP8L1"	,"OTOP2"	,"C17orf78"	,"CTAGE15P"	,"F8A1 "	,"LOC285501"	,"SAMD7"	,"ARRDC5"	,"CGB"	,"C1orf173"	,"FAM86C"	,"HIST1H2BI"	,"HIST1H4E"	,"HIST1H2AJ"	,"HIST1H2BF"	,"C2CD4B"	,"HIST1H2BE"	,"HIST1H3A"	,"RAB3D"	,"PDZD9"	,"CACNG8"	,"MORC2"	,"HIST1H2BN"	,"GPR32"	,"PCDHGB5"	,"ZNF600"	,"KRTAP6-2"	,"KRTAP20-1"	,"HIST1H2BH"	,"PCDHGA1"	,"GLIPR1L2"	,"SNAPC5"	,"OR8U9"	,"SSX7"	,"IFNA10"	,"PCDHA13"	,"PCDHGB3"	,"GSC2"	,"OR9Q2"	,"KRTAP19-7"	,"KRTAP20-2"	,"C3orf15"	,"HIST2H2AB"	,"PCDHA4"	,"PRB3"	,"PCDHGB1"	,"ZNF675"	,"TUBB4Q "	,"BCAR1"	,"ITPKC"	,"SETD1A"	,"COMTD1"	,"SPIRE1"	,"PROS1"	,"ATP13A3"	,"ANAPC4"	,"ALKBH2"	,"VSTM2L"	,"RILP"	,"HMMR"	,"ASTN1"	,"C1orf96"	,"DAPK1"	,"SH3D19"	,"C1orf116"	,"DENND4B"	,"NRIP1"	,"ZC3H7A"	,"ACAD10"	,"GDA"	,"WDR19"	,"DENND3"	,"CXCL13"	,"VEGFC"	,"FKBPL"	,"ASXL1"	,"CDC42EP5"	,"SELL"	,"MRPS16"	,"SF1"	,"NOL6"	,"DDX28"	,"DRG1"	,"C17orf71"	,"EEF1B2"	,"MAD2L1"	,"SLC27A3"	,"NDUFA7"	,"AP3S2"	,"NDUFAF4"	,"CLPB"	,"KIAA1737"	,"TOMM40L"	,"RNF6"	,"DVL2"	,"LASS2"	,"MGMT"	,"GMPR"	,"AURKB")
#var=c("HIST1H3G,HIST1H3G,HIST1H3G,TNFAIP8L1,OTOP2,C17orf78,CTAGE15P /// CTAGE6P,F8A1 /// F8A2 /// F8A3,LOC285501,SAMD7,ARRDC5,CGB /// CGB1 /// CGB2 /// CGB5 /// CGB7 /// CGB8,C1orf173,FAM86C,FAM86C,HIST1H2BI,HIST1H4E,HIST1H2AJ,HIST1H2BF,C2CD4B,HIST1H2BE,HIST1H3A,RAB3D,PDZD9,CACNG8,MORC2,HIST1H2BN,GPR32,PCDHGB5,ZNF600,KRTAP6-2,KRTAP20-1,HIST1H2BH,PCDHGA1,GLIPR1L2,SNAPC5,OR8U9,SSX7 /// SSX9,CGB /// CGB1 /// CGB2 /// CGB5 /// CGB7 /// CGB8,IFNA10,PCDHA13,PCDHGB3,GSC2,OR9Q2,KRTAP19-7,KRTAP20-2,SSX7,C3orf15,HIST2H2AB,HIST2H2AB,NRIP1,NRIP1,ZC3H7A,ACAD10,GDA,WDR19,DENND3,DENND3,CXCL13,CXCL13,VEGFC,FKBPL,ASXL1,CDC42EP5,SELL,SELL,MRPS16,MRPS16,MRPS16,SF1,SF1,SF1,NOL6,NOL6,NOL6,DDX28,DDX28,DDX28,DRG1,C17orf71,EEF1B2,EEF1B2,MAD2L1,MAD2L1,SLC27A3,NDUFA7,AP3S2,AP3S2,NDUFAF4,NDUFAF4,NDUFAF4,CLPB,CLPB,KIAA1737,TOMM40L,RNF6,DVL2,LASS2,MGMT,GMPR")
#var=c("AFFX-r2-Bs-dap-M_at","AFFX-DapX-M_at","AFFX-r2-Bs-phe-M_at","AFFX-r2-Bs-lys-M_at","11716036_x_at","11716034_a_at","11750411_x_at","AFFX-LysX-M_at","11716035_at","11744413_x_at","AFFX-PheX-M_at","11755205_a_at","AFFX-DapX-3_at","11754060_a_at","11717580_a_at","11743721_at","11744904_a_at","11754142_x_at","11757758_x_at","AFFX-ThrX-3_at","AFFX-r2-Bs-thr-3_s_at","11720342_s_at","11738958_at","11754246_x_at","11718341_s_at","11715446_a_at","11725438_a_at","11757719_x_at","11721483_at","11758299_a_at","AFFX-r2-Bs-lys-5_at","11750884_a_at","11746188_s_at","AFFX-LysX-5_at","11753560_s_at","11715461_at","11729120_x_at","11754724_a_at","AFFX-r2-Bs-dap-3_at","11734745_a_at","11755206_x_at","AFFX-r2-Bs-lys-3_at","11729119_a_at","11756380_s_at","11758300_x_at","11720757_a_at","11744077_at","11741768_a_at","11717836_a_at","11715824_at","11731050_at","11723369_at","11754314_a_at","11757077_s_at","11740046_a_at","11734993_at","11752838_s_at","11722137_a_at","11759762_a_at","11722228_at","11764273_at","11762509_x_at","11760453_x_at","11716128_at","AFFX-M27830_3_at","11746708_at","11728110_at","11717717_at","11757549_s_at","11717572_a_at","11730719_x_at","11755261_a_at","11720277_a_at","11733459_at","11726255_x_at","11753789_a_at","11757339_x_at","11758234_s_at","11726772_a_at","11758606_s_at","11734026_at","11760302_a_at","11727783_s_at","11723785_s_at","11762525_s_at","11746086_x_at","11758633_s_at","11733706_s_at","11758668_s_at","11757297_s_at","11746622_a_at","11746085_at","11739011_s_at","11717889_s_at","11746592_s_at","11716426_x_at","11759203_a_at","11715293_s_at","11716424_a_at","11752762_s_at")
esetSel=data[var,]
sel=normalize(esetSel)
heatmap.2(as.matrix(sel),ColSideColors=c(rep("green",2),rep("blue",2)),col=redgreen(75),scale="none",key=TRUE,symkey=FALSE,density.info="none",trace="none",cexRow=0.5,dendrogram="row",Colv = F)


exprset22=t(data[,-c(1)])
pmeta=read.table("Annotation.txt", header=TRUE)
species=factor(pmeta$classLabel)

PC<-prcomp(exprset22)
PCi<-data.frame(PC$x,Species=pmeta$classLabel)
ggplot(PCi,aes(x=PC1,y=PC2,col=Species, fill=Species))+
  geom_point(size=3,alpha=0.5)+ #Size and alpha just for fun
  #geom_polygon( alpha = 0.5) +
  #geom_density2d(alpha=.5) +
  #stat_ellipse()+
  scale_color_manual(values = c("green","blue","red"))+ #your colors here
  #geom_path(data=PCi, aes(x=PC1, y=PC2,colour=Species), size=1, linetype=2)+
  theme_classic()

####Volcano Plot
plot(tT$logFC, 1-tT$adj.P.Val, xlim=c(-6, 6), 
     main="Effect Identification",
     xlab="log2Ratio", ylab="1-adj.P.Val")
abline(h=0.95, col="red")

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

######Convert gene symbols
source("http://bioconductor.org/biocLite.R")

biocLite("hgu133a.db")


# load NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL = TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))
# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by = "ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order
tT <- subset(tT, select = c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", 
                            "Gene.symbol", "Gene.title"))
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=FALSE)

ncbifd <- data.frame(attr(dataTable(platf), "table"))
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
library("annotate")
library("hgu133a.db")
OUT <- select(hgu133a.db, PROBES, c("SYMBOL", "ENTREZID", "GENENAME"))

####### change gene sumbol names
data=read.table("C:\\Users\\alokum\\Desktop\\FIO\\ex.txt", header=TRUE, sep="\t")
ID <- featureNames(tT$PROBE)

probe <- read.delim("ex.txt",stringsAsFactors=F, header = T, sep="\t")
probe$probeid<-tolower(probe$probeid)
names<-read.delim("GPL15207-14509.txt", as.is=T, stringsAsFactors=F, header=T)
##insted of dataframa we are sending out the vecotr
names<-names$probeid



biocLite('biomaRt')
library(biomaRt)
vignette('biomaRt')
m <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
map <- getBM(mart = m, attributes = c("ensembl_gene_id", "PROBEID"), filters = "PROBEID", 
             values = fData(tT)$PROBE)
head(map)

library(gdata)                   # load gdata package 
#help(read.xls)                   # documentation 
d1 = read.table("C:\\Users\\alokum\\Desktop\\GPL.txt", sep="\t", header=TRUE)

merged.data <- merge(d1, d2 , by.x="ID", by.y="PROBE", all=TRUE, )




# Get the annotation GPL id (see Annotation: GPL10558)

library("annotate")
PROBES<- as.data.frame(ex[,1],drop=false)
