require(tidyverse)
require(pbapply)
require(ggprism)

#Reading data and spliting motifs into residues

tt_rank <- readRDS("data/tt_kinaseRank_updated_clust.rds")

old_cluster <- readRDS("clusters_old.rds") %>% dplyr::rename(AssignedClusterOld = AssignedCluster)
tt_rank <- readRDS("data/tt_kinaseRank_updated_clust.rds") %>% mutate(AccMod = paste(Accession, MasterMod, sep = ":")) %>%
                left_join(old_cluster, by = c("GeneMod")) 

tt_rank %>% group_by(AssignedClusterOld) %>% count()

#Splitting the motif into single amino acids for pLogo generation 
#Setting new colnames for split animo acids
myCols <- paste("A", 0:21, sep = "")
tt_aa <- tt_rank %>% 
                    separate(motif, sep = "", into = myCols) %>% 
                    select(AssignedCluster, AssignedClusterOld, K1,K2,K3, FDR_90min, FDR_40min, CDK1, CDK2, CDK5, A1:A21, enzyme_genesymbol,AccMod) %>% 
                    #drop_na() %>%
                    filter(!duplicated(AccMod))


#Functions for calculation 

#Calcualtes aminoAcid frequency
calcCount <- function(df) {
#Count the amino acids at specific position 
aa_table <- apply(df,2,function(x)table(x))

#Generate a dataframe to merge to 
temp <- data.frame(AminoAcid =  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"))

#Loop the process
for (i in 1:length(aa_table)) {
    #covert the nested table to dataframe
    temp2 <- data.frame(aa_table[i])
    #Merge
    colnames(temp2) <- c("AminoAcid", names(aa_table[i]))
    temp <- left_join(temp, temp2, by = c("AminoAcid"))

}
#Fill NA with 0
temp[is.na(temp)] <- 0 
#Move aminoAcid column
temp <- temp[-1]

#Setting rownames
rownames(temp) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

temp

}

#Calculating comparison between fast (8) and slow (7) cluster. Using cluster 7 as background

#Calcualtes aminoAcid frequency
calcFreq <- function(df) {
#Count the amino acids at specific position 
aa_table <- apply(df,2,function(x)table(x))

#Generate a dataframe to merge to 
temp <- data.frame(AminoAcid =  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"))

#Loop the process
for (i in 1:length(aa_table)) {
    #covert the nested table to dataframe
    temp2 <- data.frame(aa_table[i])
    #Merge
    colnames(temp2) <- c("AminoAcid", names(aa_table[i]))
    temp <- left_join(temp, temp2, by = c("AminoAcid"))

}
#Fill NA with 0
temp[is.na(temp)] <- 0 
#Move aminoAcid column
temp <- temp[-1]

#Calculate frequency 
temp <- temp/dim(df)[1]

#Setting rownames
rownames(temp) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

temp

}


##Empirical classification
cdk_fast <- tt_aa %>% 
                dplyr::filter(CDK1 > .563 | CDK2 > .374 | CDK5 > .452 | grepl("CDK", enzyme_genesymbol), FDR_90min < 0.05, AssignedClusterOld == 7) %>% 
                dplyr::select(A1:A21) %>% distinct()
cdk_slow <- tt_aa %>% 
                dplyr::filter(CDK1 > .563 | CDK2 > .374 | CDK5 > .452 | grepl("CDK", enzyme_genesymbol), FDR_90min < 0.05, AssignedClusterOld == 8) %>% 
               dplyr::select(A1:A21) %>% distinct()
#
#cdk_fast <- tt_aa %>% 
#                dplyr::filter(grepl("CDK", K3) | grepl("CDK", K2) | grepl("CDK", K3), FDR_90min < 0.1, AssignedCluster == 3) %>% 
#                dplyr::select(A1:A21) 
#cdk_slow <- tt_aa %>% 
#                dplyr::filter(grepl("CDK", K3) | grepl("CDK", K2) | grepl("CDK", K3), FDR_90min < 0.1, AssignedCluster == 8) %>% 
#               dplyr::select(A1:A21) 

fast_freq <- calcFreq(cdk_fast)
slow_freq <- calcFreq(cdk_slow)
fast_count <- calcCount(cdk_fast)
slow_count <- calcCount(cdk_slow)


# Calculating the probability difference with boostraping 

#Step 1 Extract boostrapped samples

#Create list to store the boostrapped motifs
fast_list <- list()
slow_list <- list()

#Use for loop to repeat loop the process N times 
N <- 200

for(i in 1:N){
    #Generate integer for sampling 
    #First sampling vector 
    s1 <- sample.int(dim(cdk_fast)[1],dim(cdk_fast)[1], replace = TRUE)
    #Extract the motifs and add to list
    fast_list[[i]] <- cdk_fast[s1,]
    }

for(i in 1:N){
    #Generate integer for sampling 
    #First sampling vector 
    s1 <- sample.int(dim(cdk_slow)[1],dim(cdk_slow)[1], replace = TRUE)
    #Extract the motifs and add to list
    slow_list[[i]] <- cdk_slow[s1,]
    }

#Step 2 Calcylate the frequency for each sample 
fast_freq <- pblapply(fast_list, calcFreq)
slow_freq <- pblapply(slow_list, calcFreq)

#Step 3 substract the slow frequencies from the fast frequencies in each sample 
freq_list <- list()

for(i in 1:N){
    freq_list[[i]] <- fast_freq[[i]] - slow_freq[[i]]
}

#Step 4 calculate mean frequence and sd
m <- matrix(ncol = 21, nrow = 20, NA) #means
#Adding row and column names 
colnames(m) <- colnames(slow_freq[[1]]);rownames(m) <- rownames(slow_freq[[1]])
s <- m # standart deviations
#Looping the samples

#Looping the residues(rows)
for(j in 1:20){
    #Looping the positions(columns)
   for(k in 1:21){
       #Looping samples
       v <- c()
       for (i in 1:N) {
        #Concatenat all frequencies for the position into one vector
        v <- c(v, freq_list[[i]][j,k])
        }
    m[j,k] <- mean(v)
    s[j,k] <- sd(v)
     }
}

#Cuting the motif down to 6 flanking residues 

m <- m[,4:18];colnames(m) <- as.character(c((7:1)*(-1),0,1:7))
s <- s[,4:18]

#Calculating p-value matrix 
p <- 1 - pnorm(abs(m), mean = 0, sd = s)
p_sig <- p < 0.05

p2 <- p.adjust(p, method  = "fdr")
p_sig <- p2 < 0.05

#Making a nice tile plot 
mdf <- as.data.frame(m) %>% rownames_to_column("AA") %>% gather(pos, diff, `-7`:`7`) %>% 
mutate(sig = as.vector(p_sig), std = as.vector(s), p = as.vector(p), p.adj= as.vector(p2),
       ) 

ggplot(mdf, aes(y = AA, x = pos, fill = diff, color = p_sig)) +
  geom_tile(color = "white",
            lwd = 1.,
            linetype = 1) +
  viridis::scale_fill_viridis(option ="D", direction = -1, begin = 0.05) +
  #coord_fixed() +
  theme_prism(border = TRUE) +
  scale_color_manual(values = c("white", "tomoto3")) +
  xlab("") + ylab("")

pheatmap::pheatmap(m, display_numbers = p_sig, cluster_cols = F, cluster_rows = F)
corrplot::corrplot(m, p.mat = p, 
                   sig.level =  0.05,
                   pch.cex = 0.9,
                    is.corr = F)
#Step 5 reformat the data so it can be used for ggplot 

#m2 <- data.frame(freq = as.vector(m), 
#                  std = as.vector(s),
#                  residue = rownames(slow_freq[[1]]), 
#                  position = as.numeric(gsub("A", "",rep(colnames(slow_freq[[1]]), each = 20)))
#                  ) %>% 
#                  mutate(pval = pnorm(-abs(freq), mean = 0 , sd = std),
#                         p.adj = p.adjust(pval, method = "fdr"),
#                         sig = ifelse(pval < 0.05, "*", "n.s"))

m2 <- mdf %>% rename(freq = diff, 
                     position = pos,
                     residue = AA) %>% 
              mutate(std = as.vector(s),
                     er = std/sqrt(99),
                     pval = pnorm(-abs(freq), mean = 0 , sd = std),
                     p.adj = p.adjust(pval, method = "fdr"),
                     mycol = ifelse(freq < 0, "a", "b"),
                     position = factor(position, levels = as.character(c((7:1)*(-1),0,1:7))))

m2 %>% filter(residue %in% c("K")) %>% 
    ggplot(aes(x = position, y = freq, fill = mycol)) +
    geom_bar(stat = "identity") + 
    geom_linerange(aes(ymin = freq - std, ymax = freq + std), size = 1) +
    #scale_x_continuous(breaks=seq(1,21,1)) +
    ggprism::theme_prism(border = F) +
    scale_fill_manual(values = hcl.colors(10, "teal")[c(7,2)]) +
    scale_y_continuous(guide = "prism_offset_minor") +
    xlab("Position") + ylab("Δ Frequency")  +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    theme(legend.position = "none")

m2 %>% filter(residue %in% c("S")) %>% 
    ggplot(aes(x = position, y = freq, fill = mycol)) +
    geom_bar(stat = "identity") + 
    geom_linerange(aes(ymin = freq - std, ymax = freq + std), size = 1) +
    #scale_x_continuous(breaks=seq(1,21,1)) +
    ggprism::theme_prism(border = F) +
    scale_fill_manual(values = hcl.colors(10, "mint")[c(7,2)]) +
    scale_y_continuous(guide = "prism_offset_minor") +
    ylab("Position") + ylab("Δ Frequency")  +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    theme(legend.position = "none")

m2 %>% filter(residue %in% c("T")) %>% 
    ggplot(aes(x = position, y = freq, fill = mycol)) +
    geom_bar(stat = "identity") + 
    geom_linerange(aes(ymin = freq - std, ymax = freq + std), size = 1) +
    #scale_x_continuous(breaks=seq(1,21,1)) +
    ggprism::theme_prism(border = F) +
    scale_fill_manual(values = hcl.colors(10, "darkmint")[c(7,2)]) +
    scale_y_continuous(guide = "prism_offset_minor") +
    ylab("Position") + ylab("Δ Frequency")  +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    theme(legend.position = "none")
