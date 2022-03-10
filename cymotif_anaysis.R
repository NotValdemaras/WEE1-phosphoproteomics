library(tidyverse)
library(ggpubr)
library(ggprism)

old_cluster <- readRDS("clusters_old.rds") %>% dplyr::rename(AssignedClusterOld = AssignedCluster)
tt_rank <- readRDS("data/tt_kinaseRank_updated_clust.rds") %>% mutate(AccMod = paste(Accession, MasterMod, sep = ":")) %>%
                left_join(old_cluster %>% dplyr::select(-MasterMod), by = c("GeneMod"))


locateCyMotifs <- function(df){
    #Function that extracts cymotifs locations
    #Amino acid seqeunce column has to be named "AAseq"

    #Function for caluclating means in the sublists row-wise
    calcMean <- function(x){
        apply(x, 1, mean)} 

    #set input to temporary variable
    temp <- df
    #Extract the locations 
    l1 <- str_locate_all(pattern = "[R].[L]", temp$AAseq) # Edit this part to change to classical RxL
    l2 <- rapply(l1, calcMean, how = "list")
    l3 <- rapply(l2, function(x) paste(x, collapse = ";"))

    #Extract the mathced Cy motifs
    c1 <- str_extract_all(pattern = "[R].[L]", temp$AAseq) # Edit this part to change to classical RxL
    c2 <- rapply(c1, function(x) paste(x, collapse = ";"))

    #Add back to dataframe
    temp$cy_location <- l3
    temp$cy_motif <- c2

    #return dataframe
    as_tibble(temp)
}

#find Cy motif locations 
cy <- tt_rank %>% select(GeneMod, MasterMod, AAseq, location) %>%
        locateCyMotifs()

closesCyMotif <- function(df){


    #Some functiosn for index extraction 
    myfunc <- function(a,b){
        which(a == b)
    }

    myfunc2 <- function(a,b){
            a[b]}

    temp <- df

    #merge site location to the cy_locations 
    c1 <- paste(temp$cy_location, temp$location, sep = ";")

    #Split the locations
    c2 <- strsplit(c1,split=';', fixed=TRUE)

    #convert to integers and sort
    c3 <- rapply(c2, function(x) sort(as.integer(x)), how = "list")

    #collapse back into one string
    c4 <- rapply(c3, function(x) paste(x, collapse=";"))
    temp$merged_location <- c4


    #Some times the motifs that contain S,T,Y in the middle of the Cy motif are on top of the detected phosphoSites
    c5 <- mapply(myfunc, strsplit(c4,split=';', fixed=TRUE),temp$location)
    is.na(c5) <- lengths(c5) == 0
    temp$siteIndexList <- c5
    #Marking the "on top" sites
    temp$CyOnSite <- ifelse(lengths(temp$siteIndexList) == 2, "YES", "NO")

    #Filtering out the the merged_location to the remove the cy motif and site duplication 

    #Split the merged_locations
    c6 <- strsplit(temp$merged_location,split=';', fixed=TRUE)
    #Remove the non-unique elements 
    c7 <- rapply(c6, function(x)unique(x), how = "list")
    #Collapsing back into one string
    c8 <- rapply(c7,function(x) paste(x, collapse = ";"))
    temp$merged_location_filtered <- c8
    #Extract the location indexes 
    c9 <- mapply(myfunc, strsplit(c8,split=';', fixed=TRUE),temp$location)
    is.na(c9) <- lengths(c9) == 0
    temp$siteIndex <- unlist(c9)


    temp$leftIndex <- temp$siteIndex - 1
    temp$rightIndex <- temp$siteIndex + 1


    #Extracting the flanking cy motif locations and replacing the out of bound locations with NA
    temp$location_left <- mapply(myfunc2, strsplit(temp$merged_location_filtered,split=';', fixed=TRUE), temp$leftIndex)
    temp$location_left[lengths(temp$location_left) == 0] <- NA
    temp$location_left <- as.numeric(unlist(temp$location_left))

    temp$location_right <- mapply(myfunc2, strsplit(temp$merged_location_filtered,split=';', fixed=TRUE), temp$rightIndex)
    temp$location_right[lengths(temp$location_right) == 0] <- NA
    temp$location_right <- as.numeric(unlist(temp$location_right))

    #Extracting the cy motifs 
    #Extracting the flanking cy motif locations and replacing the out of bound locations with NA
    temp$cy_left <- mapply(myfunc2, strsplit(temp$cy_motif,split=';', fixed=TRUE), temp$leftIndex)
    temp$cy_left[lengths(temp$cy_left) == 0] <- NA
    temp$cy_left <- unlist(temp$cy_left)

    temp$cy_right <- mapply(myfunc2, strsplit(temp$cy_motif,split=';', fixed=TRUE), temp$rightIndex)
    temp$cy_right[lengths(temp$cy_right) == 0] <- NA
    temp$cy_right <- unlist(temp$cy_right)

    temp
}

#Find the closes Cy motif to the phosphorylation site 
cy2 <- closesCyMotif(cy) %>% distinct()


#Extracting known CDK sites
cdk <- tt_rank %>% 
        filter(grepl("CDK1|CDK2", enzyme_genesymbol)) %>% 
        select(AssignedCluster, GeneMod, MasterMod, AAseq, location) %>% 
        distinct() %>%
        left_join(cy2, by = c("GeneMod", "MasterMod", "AAseq", "location")) 


#Extracting predicted CDK sites 
cdk <- tt_rank %>% 
        filter(CDK1 > .563 | CDK2 > .374 | CDK5 > .452, FDR_90min < 0.05, AssignedClusterOld %in% c(7,8))  %>% 
        select(AssignedClusterOld, GeneMod, MasterMod, AAseq, location) %>% 
        distinct() %>%
        left_join(cy2, by = c("GeneMod", "MasterMod", "AAseq", "location")) 
#Trasforming the data for plotting 


new_cdk <- cdk %>% 
            mutate(loc_left = abs(location_left - location),
                   loc_right = abs(location_right - location)) %>% 
            filter(AssignedClusterOld %in% c(7,8)) %>% 
            mutate(AssignedCluster = ifelse(AssignedClusterOld == 7, "Fast", "Slow")) %>% 
            select(-AAseq)

#Calculate some summary statistic
cy_stats <- new_cdk %>% group_by(AssignedClusterOld) %>% 
                summarise(notFound = sum(loc_left > 100, na.rm =T),
                          Found = sum(loc_left < 100, na.rm =T),
                          med = median(loc_left, na.rm = T),
                          er = mad(loc_left, na.rm = T))
cy_stats

ggplot(new_cdk %>% filter(loc_right< 80 & loc_right > 15), aes(x = loc_right, fill = AssignedCluster)) +
        #geom_density(size = 1.5) +
        geom_histogram(aes(y = stat(count) / sum(count)), bins = 10, position = "dodge") + 
        theme_prism() +
        xlab("Cy motif distance") +ylab("Frequency") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              legend.position = "top") +
        ylim(c(0,.3)) +
        scale_y_continuous(guide = "prism_offset", position = "left") +
         scale_x_continuous(guide = "prism_offset") + 
        scale_fill_manual(values = hcl.colors(n=3, palette = "viridis")) +
        xlim(c(100,0)) 


ggplot(new_cdk %>% filter(loc_left < 100), aes(x = loc_left, color = AssignedCluster)) +
  stat_ecdf(size =2) +
  scale_color_manual(values = hcl.colors(n=3, "viridis")) +
  theme_prism() +
        xlab("Cy motif distance") +ylab("Probability") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              legend.position = "top") +
        ylim(c(0,.3)) +
        scale_y_continuous(guide = "prism_offset", position = "left") +
         scale_x_continuous(guide = "prism_offset") + 
        scale_fill_manual(values = hcl.colors(n=3, palette = "viridis")) +
        xlim(c(0,100)) 

## Testing out the new approach

#        theme(axis.text.x = element_text(size = 12),
#              axis.text.y = element_text(size = 12),
#              legend.position = "top") +
#        scale_color_manual(values = hcl.colors(n=3, palette = "viridis")) +
#        ggpubr::labs_pubr() + 
#        xlim(c(100,0)) #+
#        scale_y_continuous(position = "right") 

#Adjusted function that calculates the distance of the Cy motifs phosphorylation site 

closesCyMotif2 <- function(df){

    #Assign temporary variable
    temp <- df

    #Split Cy location string into single numbers "24|100|150" -> ["24", "100", "1500"]
    c1 <- strsplit(temp$cy_location ,split=';', fixed=TRUE)

    #Convert the string list to integers 
    c2 <- lapply(c1, function(x) as.integer(x))

    #Replacing the integer(0) for rows without citations with 0, to avoid errors in the substration bellow
    c2[which(lengths(c2) == 0)] <- 0

    #Substrat the location of the phosphorylation from the cy motif location 

        #Simple substraction function to use in a for loop
        substractionFunction <- function(l, y) {
            rapply(l, function(x) x - y)}

        #Creating empty list to save the data 
        c3 <- list()

        #Looping for all the rows
        for (i in 1:length(temp$location)){
        c3[[i]] <- substractionFunction(c2[i], temp$location[i])}    

    #Collapse all the integers into one string so they can be added as a row to the dataframe 
    c4 <- rapply(c3,function(x) paste(x, collapse = ";"))

    #Add to the temporary dataframe and return
    temp$fromSite <- c4

    temp

}

cy3 <- closesCyMotif2(cy)


#Function that extract the motifs that are a specific distance from the modification 

specificDistance <- function(df, d) { 
    #df - dataframe from closesCyMotif2
    #d - mas distance to factor in motifs

    temp <- df 

    #Split the the string into integers 
    c1 <- strsplit(temp$fromSite ,split=';', fixed=TRUE)

    #Convert the string list to integers 
    c2 <- lapply(c1, function(x) as.integer(x))

    #Keep only the Cy motifs that are and a certain distance (d) from the from the phosphorylation site 
    c3 <- lapply(c2, function(x) x[abs(x) <= d])

    #Collapse all the integers into one string so they can be added as a row to the dataframe 
    c4 <- rapply(c3,function(x) paste(x, collapse = ";"))

    #Add to the temporary dataframe and return
    temp$fromSiteFiltered <- c4

    temp
}

cy4 <- specificDistance(cy3, 100)

#Merge with datase data 
temp_tt <- tt_rank %>% select(AssignedCluster,AssignedClusterOld, K1:K3, CDK1, CDK2, CDK5, FDR_90min, FDR_20min)

data <- tibble(cbind(cy4, temp_tt))

df_1 <- data %>% 
                distinct() %>%
                filter(!is.na(fromSiteFiltered)) %>% 
                separate_rows(fromSiteFiltered, sep = ";", convert = TRUE) %>% 
                filter(abs(fromSiteFiltered) > 10) %>% 
                #filter(grepl("CDK1|CDK2", enzyme_genesymbol), AssignedCluster %in% c(3,8))
                #filter(grepl("CDK", K1) | grepl("CDK", K2) | grepl("CDK", K3),   FDR_90min < 0.05, AssignedCluster %in% c(3,8))
                filter(CDK1 > .563 | CDK2 > .374 | CDK5 > .452, FDR_90min < 0.05, AssignedClusterOld %in% c(7,8)) %>%
                mutate(newgroups = case_when(fromSiteFiltered < 0 & AssignedClusterOld == 7 ~ "A",
                                             fromSiteFiltered > 0 & AssignedClusterOld == 7 ~ "B",
                                             fromSiteFiltered < 0 & AssignedClusterOld == 8 ~ "C",
                                             fromSiteFiltered > 0 & AssignedClusterOld == 8 ~ "D"))
mycol = c(rep(hcl.colors(n = 10, "darkmint")[1],2),
          rep(hcl.colors(n = 10, "darkmint")[8],2))

ggplot(df_1, aes(x = fromSiteFiltered, color = as.character(newgroups))) + 
        geom_density(size = 1.5) +
        theme_prism(border = F) +
        scale_y_continuous(guide = "prism_offset") +
        xlab("Distance, aa") + ylab("Density") +
        scale_color_manual(values = mycol) +
        theme(legend.position = "none")

df_1

ggplot(df_1, aes(x = fromSiteFiltered, fill = as.character(AssignedClusterOld))) + 
        geom_histogram(bins = 20, position = "dodge") +
        xlab("Distance, aa") + ylab("Density")

df_2 <- df_1  %>% group_by(GeneMod, AssignedClusterOld) %>% count(fromSiteFiltered) 


ggplot(df_2, aes(y = n, x = fromSiteFiltered, group = AssignedClusterOld, color = as.character(AssignedClusterOld))) + 
        geom_point(size = 3)


ggplot(df_1, aes(x = fromSiteFiltered, fill = AssignedClusterOld)) +
        #geom_density(size = 1.5) +
        geom_histogram(aes(y = stat(count) / sum(count)), bins = 20, position ="dodge") + 
        theme_prism() +
        xlab("Cy motif distance") +ylab("Frequency") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              legend.position = "top") +
        ylim(c(0,.3)) +
        scale_y_continuous(guide = "prism_offset", position = "left") +
         scale_x_continuous(guide = "prism_offset") + 
        scale_fill_manual(values = hcl.colors(n=3, palette = "viridis")) +
        xlim(c(100,0)) 

test <- ggplot_build(hist1)$data[[1]]

#
test2 <- ggplot_build(myhist)$data[[1]]


totals <- df_1 %>% group_by(AssignedClusterOld) %>% count() %>% pull(n)
test3 <- df_1 %>% group_by(AssignedClusterOld) %>% 
        mutate(bins = cut_interval(fromSiteFiltered, n = 20, width = 10)) %>% 
        count(bins) %>% mutate(freq = case_when(AssignedClusterOld == 7 ~ n/totals[1],
                                                AssignedClusterOld == 8 ~ n/totals[2]))

#Boostraping for frequency uncertainty                                                
N = 1000
temp_a <- matrix(nrow = 20, ncol = 1000)
temp_b <- matrix(nrow = 20, ncol = 1000)
a <- df_1 %>% filter(AssignedClusterOld == 7) %>% pull(fromSiteFiltered)
b <- df_1 %>% filter(AssignedClusterOld == 8) %>% pull(fromSiteFiltered)
#Taking two thirs as the sample size 

for(i in 1:N){
    #Generate integer for sampling 
    #First sampling vector 
    s1 <- sample.int(length(a),length(a), replace = TRUE)
    #Extract the motifs and add to list
    temp <- as.numeric(table(cut_interval(a[s1], n = 20, width = 10)))/totals[1]
    temp_a[,i] <- temp
    }

for(i in 1:N){
    #Generate integer for sampling 
    #First sampling vector 
    s1 <- sample.int(length(b),length(b), replace = TRUE)
    #Extract the motifs and add to list
    temp <- as.numeric(table(cut_interval(b[s1], n = 20, width = 10)))/totals[2]
    temp_b[,i] <- temp
    }

#Calculating mean and std 
freq_df <- data.frame( freq = c( apply(temp_a, 1, mean), apply(temp_b, 1, mean)),
                       er = c(apply(temp_a, 1, sd), apply(temp_b, 1, sd)),
                       AssignedClusterOld = c(rep("Fast", 20), rep("Slow", 20)),
                       bins = rep(c(-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80,90,100),2)) %>% 
                       mutate(mygroup = case_when(AssignedClusterOld == "Fast" & bins > 0 ~ "A",
                                                  AssignedClusterOld == "Fast" & bins < 0 ~ "B",
                                                  AssignedClusterOld == "Slow" & bins > 0 ~ "C",
                                                  AssignedClusterOld == "Slow" & bins < 0 ~ "D"))


ggplot(subset(freq_df, abs(bins) != 10), aes( x= bins, y = freq, color = AssignedClusterOld, group = mygroup)) + 
        geom_point(size = 2) +
        geom_linerange(aes(ymin = freq -er, ymax = freq + er), size = 1) +
        geom_line(size = 1) + 
        theme_prism(border = F) +
        scale_y_continuous(guide = "prism_offset") +
        scale_color_manual(values = hcl.colors(n=6, "viridis")[c(1,4)])+
        xlab("Distance, aa") + ylab("Frequency") +
        theme(legend.position = "top")

ggplot(test3, aes(x = bins, y = n, color = AssignedClusterOld, group = AssignedClusterOld)) + 
    geom_point() +
    geom_line()
sum(is.na(df_1$fromSiteFiltered))

test3


### Adding disorcer and conservation prediction to select only one Cy motif per site (run the CyDistance.R script first to set up the dataframe)
cy4 <- cy4 %>% mutate(substrate = tt_rank$Accession)
c1 <- cy4 %>% separate_rows(cy_location, sep = ";")  %>% pull(cy_location)
cy5 <- cy4 %>% separate_rows(fromSite, sep = ";") %>% mutate(cy_loc = c1, cyAcc = paste(substrate, c1, sep = "_"))


#Add the conservation and disorder prediction for each Cy motif

test <- list()

for (i in unique(cy5$substrate)){
    prot <- dis_cons[[i]] # extract substrate information for specific substrate 
    #loop over all potential 
    temp_dis <- c(); temp_cons <- c(); temp_aa <- c()
    for (j in cy5$cy_loc[cy5$substrate == i]){ 
        j <- as.numeric(j)
        temp_dis <- c(temp_dis, mean(as.numeric(prot$property[c(j-1,j,j+1)])))
        temp_cons <- c(temp_cons,mean(as.numeric(prot$conservation_level[c(j-1,j,j+1)])))
        temp_aa <- c(temp_aa, paste(prot$AminoAcid[c(j-1,j,j+1)], collapse = ""))          
    }
    temp_info <- data.frame(substrate = i, cy_loc = cy5$cy_loc[cy5$substrate == i],dis = temp_dis, cons = temp_cons, aa = temp_aa) 
    test[[i]] <- temp_info
}

test2 <- bind_rows(test)


cy6 <- cy5 %>% 
                left_join(test2, by = c("substrate", "cy_loc")) %>%
                mutate(cy_id = paste(substrate, cy_loc, sep = "_"),
                       fromSite = as.numeric(fromSite)) %>% 
                distinct()

#Filtering out motif in the range between 10 and 100 AA 

cdkSites <- tt_rank %>% filter(CDK1 > .563 | CDK2 > .374 | CDK5 > .452, FDR_90min < 0.05, AssignedClusterOld %in% c(7,8)) %>% 
            select(GeneMod,MasterMod, AssignedClusterOld) %>% distinct()
cy7 <- cy6 %>% filter(fromSite < 100 & fromSite > 10) %>% 
        left_join(cdkSites, by = c("GeneMod", "MasterMod")) %>% filter(!is.na(AssignedClusterOld))

table(cy7$AssignedClusterOld)

cy_cons <- cy7 %>% group_by(GeneMod, MasterMod) %>% filter(cons == max(cons)) 

cy_dis <- cy7 %>% group_by(GeneMod, MasterMod) %>% filter(dis == max(dis)) 


ggplot(cy_dis, aes(x = fromSite, fill = as.character(AssignedClusterOld))) + 
        geom_histogram(bins = 20, position = "dodge") 
        #geom_density()

