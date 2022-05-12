# Microbiome final project

# Environment ---------

  #libraries
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("MMUPHin")
  BiocManager::install("curatedMetagenomicData")
  
  library(MMUPHin)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(vegan)
  library(curated)

  um_demog_v34 <- read.table(file = '../data/UgandaMaternalV3V4.16s_DADA2.sample_details.tsv', sep = '\t', header = TRUE)
  um_taxon_v34 <- read.table(file = '../data/UgandaMaternalV3V4.16s_DADA2.taxon_abundance.tsv', sep = '\t', header = TRUE)
  
  um_demog_v12 <- read.table(file = '../data/UgandaMaternalV1V2.16s_DADA2.sample_details.tsv', sep = '\t', header = TRUE)
  um_taxon_v12 <- read.table(file = '../data/UgandaMaternalV1V2.16s_DADA2.taxon_abundance.tsv', sep = '\t', header = TRUE)
  
  View(um_demog_v12)
  View(um_demog_v34)
  View(um_taxon_v12) 
  View(um_taxon_v34) 
  
# Data Cleaning ---------------------
  
  #fix city/village situation
  um_demog_v12$City..village..or.region <- ifelse(um_demog_v12$City..village..or.region=="Mbale, Uganda", "Mbale", "Mbarara")
  
  #Construct a feature-abundance matrix with the taxon abundances
  
  #Make all NAs -> 0
  um_taxon_v12[is.na(um_taxon_v12)] <- 0
  um_taxon_v34[is.na(um_taxon_v34)] <- 0
  
  #Create matrix of abundances currently the columns all add up differently
  #sum(um_taxon$X5025.SRR12919592)
  #sum(um_taxon$X5031.SRR12919543)
  
  #Need to divide all cells in the column by the column total
  #Come back to this - unsure if the absolute abundance will work but going to try it 
  #If not will have to do relative abundance by hand
  
  #Distance/dissimilarity matrix for MMuphin - Bray-Curtis distances
  #flip so that species are columns
  um_taxon_t12 <- as.data.frame(t(um_taxon_v12[-1]))
  um_taxon_t34 <- as.data.frame(t(um_taxon_v34[-1]))
  colnames(um_taxon_t12) <- um_taxon_v12[,1]
  colnames(um_taxon_t34) <- um_taxon_v34[,1]
  View(um_taxon_t12)
  View(um_taxon_t34)
  
  #Set up batch argument
  um_demog_v12$batch <- "v1v2"
  um_demog_v34$batch <- "v3v4"
  
  #check which columns are the same
  library(janitor)
  compare_df_cols(um_taxon_t12, um_taxon_t34, return="match") #at least 2/3 or more (get.print cutoff are the same of the 34 options)
  
  #full join the two
  um_taxon_both <- full_join(um_taxon_t12, um_taxon_t34)
  View(um_taxon_both)
  rownames(um_taxon_both) <- c(rownames(um_taxon_t12), rownames(um_taxon_t34))
  #all NAs -> 0
  um_taxon_both[is.na(um_taxon_both)] <- 0
  
  #full join the two demographics data frames too
  um_demog_both <- full_join(um_demog_v12, um_demog_v34)
  
  #create duplicated matrix of v1v2
  #make the names possible
  um_taxon_t12_a <- um_taxon_t12
  rownames(um_taxon_t12_a) <- paste(rownames(um_taxon_t12), "a", sep="")
  um_taxon_t12_b <- um_taxon_t12
  rownames(um_taxon_t12_b) <- paste(rownames(um_taxon_t12), "b", sep="")
  
  v1v2_times2 <- rbind(um_taxon_t12_a, um_taxon_t12_b)
  
  #create duplicated matrix of v3v4
  um_taxon_t34_a <- um_taxon_t34
  rownames(um_taxon_t34_a) <- paste(rownames(um_taxon_t34), "a", sep="")
  um_taxon_t34_b <- um_taxon_t34
  rownames(um_taxon_t34_b) <- paste(rownames(um_taxon_t34), "b", sep="")
  
  v3v4_times2 <- rbind(um_taxon_t34_a, um_taxon_t34_b)
  
  #Just v3v4 matrix - compare Mbarara and Mbale
  umtt_mat34 <- as.matrix(um_taxon_t34)
  #Both matrix
  umtt_matall <- as.matrix(um_taxon_both)
  #Double v1v2 matrix
  umtt_mat12x2 <- as.matrix(v1v2_times2)
  #Double v3v4 matrix
  umtt_mat34x2 <- as.matrix(v3v4_times2)

  #Compute distance matrix
  #Distance matrix just v3v4
  um_dist34 <- vegdist(umtt_mat34, method="bray")
  #Distance matrix both
  um_dist_both <- vegdist(umtt_matall, method="bray")
  #Distance matrix v1v2
  um_dist12 <- vegdist(umtt_mat12x2, method="bray")
  #Distance matrix double v3v4
  um_dist34_double <- vegdist(umtt_mat34x2, method="bray")
  
  #need the rownames of the metadata to match the colnames of the taxon abundance data
  #for both
  um_demog_both$rowname <- paste("X", um_demog_both$X, sep="")
  rownames(um_demog_both) <- um_demog_both$rowname
  #for v3v4
  um_demog_v34$rowname <- paste("X", um_demog_v34$X, sep="")
  rownames(um_demog_v34) <- um_demog_v34$rowname
  #for v1v2
  #duplicate the data frame with different batch labels
  um_demog_v12_a <- um_demog_v12
  um_demog_v12_a$batch <- 1
  um_demog_v12_a$rowname <- paste("X", um_demog_v12_a$X, "a", sep="")
  rownames(um_demog_v12_a) <- um_demog_v12_a$rowname
  
  um_demog_v12_b <- um_demog_v12
  um_demog_v12_b$batch <- 2
  um_demog_v12_b$rowname <- paste("X", um_demog_v12_b$X, "b", sep="")
  rownames(um_demog_v12_b) <- um_demog_v12_b$rowname
  
  um_demog_v12_double <- rbind(um_demog_v12_a, um_demog_v12_b)
  
  #v3v4 duplicate data
  um_demog_v34_a <- um_demog_v34
  um_demog_v34_a$batch <- 1
  um_demog_v34_a$rowname <- paste("X", um_demog_v34_a$X, "a", sep="")
  rownames(um_demog_v34_a) <- um_demog_v34_a$rowname
  
  um_demog_v34_b <- um_demog_v34
  um_demog_v34_b$batch <- 2
  um_demog_v34_b$rowname <- paste("X", um_demog_v34_b$X, "b", sep="")
  rownames(um_demog_v34_b) <- um_demog_v34_b$rowname
  
  um_demog_v34_double <- rbind(um_demog_v34_a, um_demog_v34_b)
  

#MMuphin: Discrete community ----------------------
  
  #Primary results: v1v2 compared against itself
  
  fit_discrete_um12 <- discrete_discover(D=um_dist12,
                                         batch = "batch",
                                         data=um_demog_v12_double,
                                         control = list(k_max = 8))
  
  #create the graph 
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_um12$internal_mean[,1],
                         se = fit_discrete_um12$internal_se[,1],
                         type = "internal")
  
  internal %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (V1V2)")
  
  #Primary results: v3v4 compared against itself
  
  fit_discrete_um34d <- discrete_discover(D=um_dist34_double,
                                         batch = "batch",
                                         data=um_demog_v34_double,
                                         control = list(k_max = 8))
  
  #create the graph 
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_um34d$internal_mean[,1],
                         se = fit_discrete_um34d$internal_se[,1],
                         type = "internal")
  
  internal %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (V3V4)")
  
  #Secondary results: v1v2 vs v3v4
  fit_discrete_umall <- discrete_discover(D=um_dist_both,
                                          batch = "batch",
                                          data=um_demog_both,
                                          control = list(k_max = 8))
  
  #create comparison for v3v4 internal and external variation
  study_id = "v3v4"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_umall$internal_mean[, study_id],
                         se = fit_discrete_umall$internal_se[, study_id],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_umall$external_mean[, study_id],
                         se = fit_discrete_umall$external_se[, study_id],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure in vaginal microbiome (V3V4)")
  
  #Do the same with v1v2
  study_id = "v1v2"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_umall$internal_mean[, study_id],
                         se = fit_discrete_umall$internal_se[, study_id],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_umall$external_mean[, study_id],
                         se = fit_discrete_umall$external_se[, study_id],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure in vaginal microbiome (V1V2)")
  
  #Secondary results: city/village/region
  
  #Compare v1v2 to v3v4

  #compare city/village/region across v3 and v4
  #comparing v1v2 to v3v4 (gives internal for v3v4)
  fit_discrete_umallbycity <- discrete_discover(D=um_dist_both,
                                                batch = "City..village..or.region",
                                                data=um_demog_both,
                                                control = list(k_max = 8))
  
  #Compare cities
  
  #create comparison for Mbarara internal and external variation
  batch = "Mbarara"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_umallbycity$internal_mean[,2],
                         se = fit_discrete_umallbycity$internal_se[,2],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_umallbycity$external_mean[,2],
                         se = fit_discrete_umallbycity$external_se[,2],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (Mbarara, v1v2 compared to v3v4)")
  
  #create comparison for Mbale internal and external variation
  batch = "Mbale"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_umallbycity$internal_mean[,1],
                         se = fit_discrete_umall$internal_se[,1],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_umallbycity$external_mean[,1],
                         se = fit_discrete_umall$external_se[,1],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (Mbale, v1v2 compared to v3v4)")
  
  
  #Again with just v3v4 data
  
  #create comparison for Mbarara internal and external variation
  batch = "Mbarara"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_um34$internal_mean[,2],
                         se = fit_discrete_umall$internal_se[,2],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_um34$external_mean[,2],
                         se = fit_discrete_umall$external_se[,2],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (Mbarara, v3v4)")
  
  #create comparison for Mbale internal and external variation
  batch = "Mbale"
  
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_um34$internal_mean[,1],
                         se = fit_discrete_um34$internal_se[,1],
                         type = "internal")
  
  external <- data.frame(K = 2:8,
                         statistic = fit_discrete_um34$external_mean[,1],
                         se = fit_discrete_um34$external_se[,1],
                         type = "external")
  
  rbind(internal, external) %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (Mbale, v3v4)")
  
  #looking at a doubled up version of just v1v2
  fit_discrete_um12 <- discrete_discover(D=um_dist12,
                                         batch = "batch",
                                         data=um_demog_v12_double,
                                         control = list(k_max = 8))
  
  #create the graph 
  internal <- data.frame(K = 2:8,
                         statistic = fit_discrete_um12$internal_mean[,1],
                         se = fit_discrete_um12$internal_se[,1],
                         type = "internal")
  
  internal %>% 
    ggplot(aes(x = K, y = statistic, color = type)) +
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                  position = position_dodge(width = 0.5), width = 0.5) +
    ggtitle("Evaluation of discrete structure (V1V2)")
                                  
#MMuPhin: Continuous community -------  
  
  #Have to go back and set up feature by abundance tables
  #Just do a PCA in continuous - use stats princomp and prcomp
  
  #v3v4 - doing this first and checking it works
  rownames(um_taxon_v34) <- um_taxon_v34$X
  um_taxon_v34_fix <- um_taxon_v34[,-1]
  umtt_mat34_fix <- as.matrix(um_taxon_v34_fix)
  
  #v3v4 - goes back to feature/sample abundances
  fit_continuous_v3v4 <- continuous_discover(feature_abd = umtt_mat34_fix,
                                        batch = "City..village..or.region",
                                        data = um_demog_v34,
                                        control = list(var_perc_cutoff = 0.5,
                                                       verbose = FALSE))
  
  #all
  
  
  
  
  