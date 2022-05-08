library("vegan")
library("cluster")
library("biomformat")
library("githubinstall")
# githubinstall("biom")
library("biom")
library("fpc")
# githubinstall("phyloseq")
library("phyloseq")
library("ggplot2")
library("ggbreak") 
library("patchwork")
library("grid")

biom_file <- read_biom("/Users/NEK/Downloads/v1v2abs.biom")
biom_file_cts <- as.matrix(biom_data(biom_file))
biom_file_cts <- t(biom_file_cts)
rownames(biom_file_cts)
map <- read.table("/Users/NEK/Downloads/UgandaMaternalV1V2.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)
common.ids <- intersect(rownames(map), rownames(biom_file_cts))
biom_file_cts <- biom_file_cts[common.ids,]
map <- map[common.ids,]
bc <- vegdist(biom_file_cts, method = "bray")
bc
# ps_bc <- prediction.strength(bc, Gmin = 2, Gmax = 10, M = 100)
# plot(2:10, ps_bc$mean.pred[2:10], type = "b", ylim = c(0,1), col = "blue", xlab = "Number of clusters", ylab = "Prediction strength")

p <- pam(bc, 2) # 2 clusters optimal 
pc <- cmdscale(bc, 2)
# plot(pc[,1], pc[,2], col = c("red", "orange", "blue")[p$clustering])

clusterz <- p$clustering
cluster1 <- names(clusterz[clusterz == 1])
length(cluster1) # 58
cluster2 <- names(clusterz[clusterz == 2])
length(cluster2) # 41

# relative abundance plots
choiceclust <- cluster2
level <- "Genus"

biom_file_cts <- read.csv("/Users/NEK/Downloads/splittaxa_relativeabundances_V1V2.csv", row.names = 1)
colnames(biom_file_cts) <- gsub("X", "", colnames(biom_file_cts))
OTU <- otu_table(biom_file_cts, taxa_are_rows = TRUE)

map <- read.table("/Users/NEK/Downloads/UgandaMaternalV1V2.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)
clusts <- p$clustering # add clusters from previous analysis
map$Cluster <- factor(clusts)

taxes <- read.csv("/Users/NEK/Downloads/taxonomy_tablev1v2.csv", row.names = 1)
taxes_table <- tax_table(as.matrix(taxes))
taxes_table[taxes_table == ""] <- NA

biom_file_cts_2 <- t(t(biom_file_cts)[choiceclust,]) # 58 for cluster 1
physeq2 <- merge_phyloseq(otu_table(biom_file_cts_2, taxa_are_rows = TRUE), sample_data(map), taxes_table)
p <- plot_bar(physeq2, fill = level) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  guides(fill = guide_legend(nrow = 30, byrow = TRUE))
  # theme(legend.position = "none")
# p

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p)
grid.draw(mylegend)

####################################################
# using vegan (http://metagenome.cs.umn.edu/microbiomecodebrowser/data/globalgut-66-adults/prediction.strength.html)
# read in biom file
biom_file <- read_biom("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom")

# extract OTU counts from biom table
# TODO: is it a problem that we haven't collapsed OTU counts into genus (L6) and phylum (L2)?
biom_file_cts <- as.matrix(biom_data(biom_file))

# transpose so that rows are samples and columns are genera
biom_file_cts <- t(biom_file_cts)
rownames(biom_file_cts)

# read in map
map <- read.table("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)

# find the overlapping samples
common.ids <- intersect(rownames(map), rownames(biom_file_cts))

# get just the overlapping samples
biom_file_cts <- biom_file_cts[common.ids,]
map <- map[common.ids,]

# calculate Bray-Curtis distances
bc <- vegdist(biom_file_cts, method = "bray")
bc

# run prediction strength analysis on Bray-Curtis table using 100 random splits
ps_bc <- prediction.strength(bc, Gmin = 2, Gmax = 10, M = 100)

# note that the prediction strength for 1 cluster is trivially 1, which is automatically included if GMin > 1

# plot prediction strength
plot(2:10, ps_bc$mean.pred[2:10], type = "b", ylim = c(0,1), col = "blue", xlab = "Number of clusters", ylab = "Prediction strength")

# run partitioning around medoids (PAM) clustering with 3 clusters
p <- pam(bc, 3)

# PCoA coordinates of Bray-Curtis distances
pc <- cmdscale(bc, 3)

# plot PCoA colored by cluster, with countries shown by shape (can toggle between different PCoAs)
plot(pc[,1], pc[,2], col = c("red", "orange", "blue")[p$clustering])

# TODO: label clusters with phenotypes

##########################
# using phyloseq (https://web.stanford.edu/class/bios221/labs/phyloseq/lab_phyloseq.html and https://joey711.github.io/phyloseq/import-data.html)

clusterz <- p$clustering

cluster1 <- names(clusterz[clusterz == 1])
length(cluster1) # 27

cluster2 <- names(clusterz[clusterz == 2])
length(cluster2) # 40

cluster3 <- names(clusterz[clusterz == 3])
length(cluster3) # 30


#####
biom_file_cts <- read.csv("/Users/NEK/Downloads/splittaxa_relativeabundances_V3V4.csv", row.names = 1)
colnames(biom_file_cts) <- gsub("X", "", colnames(biom_file_cts))
OTU <- otu_table(biom_file_cts, taxa_are_rows = TRUE)

map <- read.table("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)
clusts <- p$clustering # add clusters from previous analysis
map$Cluster <- factor(clusts)

taxes <- read.csv("/Users/NEK/Downloads/taxonomy_table.csv", row.names = 1)
taxes_table <- tax_table(as.matrix(taxes))
taxes_table[taxes_table == "0"] <- NA

choiceclust <- cluster3
level <- "Genus"

biom_file_cts_2 <- t(t(biom_file_cts)[choiceclust,])
physeq2 <- merge_phyloseq(otu_table(biom_file_cts_2, taxa_are_rows = TRUE), sample_data(map), taxes_table)
plot_bar(physeq2, fill = level) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  # theme(legend.position = "none")

#####

biom_file <- read_biom("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.biom")
biom_file_cts <- as.matrix(biom_data(biom_file))
OTU <- otu_table(biom_file_cts, taxa_are_rows = TRUE)

map <- read.table("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)
clusts <- p$clustering # add clusters from previous analysis
map$Cluster <- factor(clusts)

taxes <- read.csv("/Users/NEK/Downloads/taxonomy_table.csv", row.names = 1)
taxes_table <- tax_table(as.matrix(taxes))
taxes_table[taxes_table == "0"] <- NA

physeq1 <- merge_phyloseq(OTU, sample_data(map), taxes_table)

# get distribution of microbes in given cluster
choiceclust <- cluster3
level <- "Order"

biom_file_cts_2 <- t(t(biom_file_cts)[choiceclust,])
biom_file_cts_2[biom_file_cts_2 == 0] <- NA

map <- read.table("/Users/NEK/Downloads/UgandaMaternalV3V4.16s_DADA2.sample_details.tsv", sep = "\t", head = T, row.names = 1)
clusts <- p$clustering # add clusters from previous analysis
map$Cluster <- factor(clusts)

map <- map[choiceclust,]
physeq2 <- merge_phyloseq(otu_table(biom_file_cts_2, taxa_are_rows = TRUE), sample_data(map), taxes_table)
plot_bar(physeq2, fill = level) + theme(legend.position = "none") + scale_y_log10() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# plot clusters
plist <- vector("list", length(2))

Dist1 <- distance(physeq1, method = "bray")
MDS1  <- ordinate(physeq1, "MDS", distance = Dist1)
p1 <- plot_ordination(physeq1, MDS1, color = "Cluster", shape = "Cytmogealovirus") # can add color = "HIV", for example
p1 <- p1 + ggtitle("MDS using Bray-Curtis") + theme_minimal()
plist[[1]] <- p1

Dist2 <- distance(physeq1, method = "jsd")
MDS2 <- ordinate(physeq1, "MDS", distance = Dist2)
p2 <- plot_ordination(physeq1, MDS2, color = "Cluster", shape = "HIV")
p2 <- p2 + ggtitle("MDS using Jenssen-Shannon divergence")
plist[[2]] <- p2

plist[[1]]
plist[[2]]


