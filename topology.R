###############################
## Network Topology Analysis ##
###############################

library("igraph")

gene.coexpression.network <- read.graph(file="gene_coexpression_network_095.gml",format="gml")

## Photoperiodic flowering gene (hd17) neighbours at distance 3 

neighbour.2 <- gene.coexpression.network[,] %*% gene.coexpression.network[,]
neighbour.3 <- gene.coexpression.network[,] %*% neighbour.2

hd17.1 <- names(which(gene.coexpression.network["LOC_Os06g05060",] > 0))
hd17.2 <- names(which(neighbour.2["LOC_Os06g05060",]>0))
hd17.3 <- names(which(neighbour.3["LOC_Os06g05060",]>0))

neighbour.hd17.3 <- unique(c(hd17.1,hd17.2,hd17.3))
length(unique(neighbour.hd17.3))
length(hd17.3)

intersect(hd17.3,unique(neighbour.hd17.3))
neighbour.hd17.3[!(neighbour.hd17.3 %in% hd17.3)]

write.table(unique(neighbour.hd17.3),file="neighbours_hd17_distancia3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

## Free-scale network because its degree distribution follows a power law
network.degree.distribution <- degree.distribution(gene.coexpression.network)
help("power.law.fit")
fit.scale.free <- power.law.fit(network.degree.distribution)
fit.scale.free[["KS.p"]]

network.degrees <- degree(gene.coexpression.network)
network.degrees
degree.histogram <- hist(network.degrees,freq=FALSE,col="blue",xlab="Node degree", ylab="Probability",main="Degree distribution")

network.degrees <- degree(gene.coexpression.network)
degree.frequencies <- table(network.degrees)10
degree.frequencies.no.0 <- degree.frequencies[-1]

log10.degrees.frequencies <- log10(degree.frequencies.no.0)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies.no.0)))
plot(log10.node.degrees)

# Linear regression
help(lm)
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r)

## Clustering coefficient 
network.clustering.coefficient <- transitivity(gene.coexpression.network,type="global")
average.path.length(gene.coexpression.network, directed=FALSE)
## 12.74602 - nodes highly accessible because of the presence of highly connected nodes or hubs  


gene.coexpression.network
402442/38361
402442 - 10*38360
sum(c(rep(10,38360),18842))

## Generate random graphs and calculate medium clustering coefficient
number.of.added.edges <- c(rep(10,38360),18842)
random.scale.free.graph <- barabasi.game(n=38361,out.seq=number.of.added.edges,directed=FALSE)
transitivity(random.scale.free.graph)
help(transitivity)
## 0.002611862

clustering.coefficients <- vector(length=1000)

for(i in 1:1000)
{
  print(i)
  random.scale.free.graph <- barabasi.game(n=38361,out.seq=number.of.added.edges,directed=FALSE)
  clustering.coefficients[i] <- transitivity(random.scale.free.graph)
}

sum(clustering.coefficients > network.clustering.coefficient) / 1000
## 0
network.clustering.coefficient
## 0.7658085 - we can argue with certainty that it's a Free-scale network 


## Extract highly connected nodes or hubs 

network.hub.scores <- hub.score(gene.coexpression.network)
hub.score.attributes <-network.hub.scores[["vector"]]
write.table(hub.score.attributes,file="hub_score_attributes.txt")

sorted.hub.score.attributes <- sort(hub.score.attributes,decreasing=TRUE)
write.table(sorted.hub.score.attributes,file="hubs_sorted.txt")
network_hubs_GO <- names(sorted.hub.score.attributes[1:1000])
network_hubs_table <- matrix(c(network_hubs,rep("hub",length(network.hubs))),nrow=1000,ncol=2)
write.table(network_hubs_GO,file="network_hubs.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
