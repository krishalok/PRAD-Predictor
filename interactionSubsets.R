require(ggplot2)
require(reshape2)
require(devtools)
require(igraph)
require(NetIndices)

# --- Bascompte et al. Quantitative network --------------

bmatrix <- read.csv("./Data/interactionsmarine.csv", header = TRUE, row.names = 1)
bmatrix <- as.matrix(bmatrix)
bgraph <- graph.adjacency(bmatrix, mode = "directed", weighted = TRUE)

plot.igraph(bgraph, layout = layout.circle, vertex.size = 1, edge.arrow.size = .25, vertex.label = NA)

bgraph

bel <- get.edgelist(bgraph)
bweights <- E(bgraph)$weight

quantile(bweights)
q75 <- which(bweights > .00239)
q50 <- which(bweights > .00030243 & bweights <= .00239)
q25 <- which(bweights > .000033 & bweights <= .00030243)
q0 <- which(bweights > 0 & bweights <= .000033)

b75 <- graph.edgelist(bel[q75,])
b50 <- graph.edgelist(bel[q50,])
b25 <- graph.edgelist(bel[q25,])
ball <- graph.edgelist(bel[q0,])

source_url("https://raw.github.com/jjborrelli/Ecological-Networks/master/FoodWebs/Rscripts/web_functions.R")

bsubs <- list(b75, b50, b25, ball)
motif.counts <- motif_counter(bsubs, webs = c("b75", "b50", "b25", "ball"))
null_counts <- null_motifs(bsubs, graph.names = c("b75", "b50", "b25", "ball"), sample = 100, iter = 5000)
null_counts <- split(null_counts, null_counts$web)
#null_counts <- melt(null_counts, id.vars = c("web", "s1", "s2", "s3", "s4", "s5", "d1", "d2", 
#                                             "d3", "d4", "d5", "d6", "d7", "d8"))

nulls <- list(null_counts[[3]][,2:14], null_counts[[2]][,2:14],
              null_counts[[1]][,2:14], null_counts[[4]][,2:14])

null.mean <- t(sapply(nulls, colMeans))
null.sd <- t(sapply(nulls, FUN = function(x){apply(x, 2, sd)}))
plot(t((motif.counts[,2:14] - null.mean) / null.sd)[,4])

bas.zscor <- t((motif.counts[,2:14] - null.mean) / null.sd)
bas.zscor[is.na(bas.zscor)] <- 0

#write.csv(bas.zscor, file = "basc_zscore.csv")
bas.zscor <- read.csv("~/Desktop/GitHub/Quantitative-Structure/Tables/basc_zscore.csv", row.names = 1)
colnames(bas.zscor) <- c("b75", "b50", "b25", "ball")
z.norm <- apply(bas.zscor, 2, FUN = function(x){x/sqrt(sum(x^2))})
plot(bas.zscor[,1], ylim = c(-10, 16), typ = "b")
points(bas.zscor[,2], col = "blue", typ = "b")
points(bas.zscor[,3], col = "red", typ = "b")
points(bas.zscor[,4], col = "orange", typ = "b")
abline(h = 0)


m <- melt(z.norm)
s <- c(9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 6, 7, 8)
s2 <- rep(s, 4)
m2 <- cbind(m, s2)
g <- ggplot(m, aes(x = Var1, y = value, col = Var2)) + geom_point(size = 3) 
g + geom_line(data = m2, aes(x = s2, y = value, col = Var2))


## ---- INDICES ------------------- 

bwebind <- get_fw_indices(adj.list = lapply(bsubs, get.adjacency, sparse = F), 
                          graphs = bsubs, web = c("b75", "b50", "b25", "ball"))

bnodeprops <- get_node_properties(adj.list = lapply(bsubs, get.adjacency, sparse = F), 
                                  web = c("b75", "b50", "b25", "ball"))


troplot <- ggplot(bnodeprops, aes(x = TL)) 
troplot <- troplot + geom_histogram(aes(y = ..density..), binwidth = .2, colour = "black", fill = "white") 
troplot + facet_grid(L1 ~ .) + geom_density(alpha = .2, fill = "#FF6666")


## --- Rezende Data

rdata <- read.csv("rezendeDATA.csv", header = TRUE)

redges <- rdata[rdata$Int == 1,]
redges <- redges[-1974,]

r.elist <- cbind(pred = as.character(redges$Pred), prey = as.character
                 (redges$Prey), weight = redges$Strength)

plot.igraph(graph.edgelist(r.elist[,1:2]), layout = layout.circle, vertex.size = 1, edge.arrow.size = .25, vertex.label = NA)

rweights <- as.numeric(r.elist[,3])
rq <- quantile(as.numeric(r.elist[,3]))
q75 <- which(rweights > rq[3])
q50 <- which(rweights > rq[2])
q25 <- which(rweights > rq[1])
q0 <- which(rweights > 0)

r75 <- graph.edgelist(r.elist[q75,1:2])
r50 <- graph.edgelist(r.elist[q50,1:2])
r25 <- graph.edgelist(r.elist[q25,1:2])
rall <- graph.edgelist(r.elist[q0,1:2])


rez <- list(r75, r50, r25, rall)
names(rez) <- c("r75", "r50", "r25", "rall")

motif_counter(rez, webs = names(rez))

motif_counter(list(graph.edgelist(r.elist[,1:2])), 1)

## --- Ulanowicz Quantitative Data ---- 
#setwd("~/Dropbox/Food Web Database/Ecosystem Flow/Ulan_Edges")

## Read in data
ULANwebnames<- c("baltic", "charca", "chesapeake", "chesapeakemeso", "crystala", "crystalb", "everglades",
                 "flbay", "lowerches", "middlechesa", "mondego", "narraganset", "rhode", "stmarks", "ythan")
networks.list <- list()
for(i in 1:15){
  networks.list[[i]] <- read.csv(paste(ULANwebnames[i], ".csv", sep = ""), row.names = 1)
}
names(networks.list) <- ULANwebnames

ulan.graphs <- list()
ulan.mat <- list()
for(i in 1:15){
  ulan.graphs[[i]] <- graph.edgelist(as.matrix(networks.list[[i]][1:2]))
  ulan.mat[[i]] <- get.adjacency(ulan.graphs[[i]], sparse = F)
}


ulan.motifs <- motif_counter(ulan.graphs, webs = ulan.names)
ulan.mot <- ulan.motifs[,2:14]

## --- Permutation ------------------------------
perm.ulan <- web_permutation(ulan.mat, fixedmar = "both", 100)
perm.ulan.mot <- lapply(perm.ulan, function(x){x[,2:14]})
perm.means <- t(sapply(perm.ulan.mot, colMeans))
perm.sd <- t(sapply(perm.ulan.mot, function(x){apply(x, 2, sd)}))

## --- Compare -----
z1 <- (ulan.mot - perm.means)/perm.sd
zN <- apply(z1, 2, function(x){x/sqrt(sum(x^2))})
rownames(zN) <- ULANwebnames

mUlan <- melt(zN)
s <- c(9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 6, 7, 8)
s2 <- rep(s, each=15)
mUlan2 <- cbind(mUlan, s2)
g <- ggplot(mUlan, aes(x = Var2, y = value, col = Var1)) + geom_point() 
g + geom_line(data = mUlan2, aes(x = s2, y = value, col = Var1))
