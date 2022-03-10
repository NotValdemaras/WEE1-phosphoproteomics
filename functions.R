require(ellipse)
require(tidyverse)
require(stringr)
#helper functions 

#Extract Accession and abundance data from PD export for a specific sample
PDextract <- function(df, sample, remove, myOrder) {
    #df <- data exporter from PD
    #sample <- the sample number, e.g. 'F5', 'F6' and etc FRACTION LABEL, NOT CONDITION LABEL
    #remove <- columns to be removed, e.g. 'EMPTY'
    l <- paste('Abundance\\.*', sample, sep ='') # generate string to search
    i <- grepl(l, names(df)) # search the columns for the string
    temp1  <- df[,i] # index the dataframe
    temp2 <- temp1[,!(grepl(remove, names(temp1)))] # remove empty samples, the
    temp3 <- df$Accession
    rownames(temp2) <- temp3
    first <- str_extract_all(names(temp2), "Sample\\.\\..*")
    second <- gsub("Sample\\.\\.", "", first)
    names(temp2) <- second
    temp2 <- temp2[myOrder]
    temp2
}

#Extract Peptide ID, annotated sequence and abundance data from PD export for a specific sample
PeptideExtract <- function(df, sample, remove, myOrder) {
    #df <- data exporter from PD
    #sample <- the sample number, e.g. 'F5', 'F6' and etc. 
    #remove <- columns to be removed, e.g. 'EMPTY'
    l <- paste('Abundance\\.*', sample, sep ='') # generate string to search
    i <- grepl(l, names(df)) # search the columns for the string
    temp1  <- df[,i] # index the dataframe
    temp2 <- temp1[,!(grepl(remove, names(temp1)))] # remove empty samples, the
    temp3 <- df$Peptide.Groups.Peptide.Group.ID # unique ID for peptide
    temp4 <- df$Annotated.Sequence # the annotated sequence (non-unique)
    temp4 <- gsub("\\[.\\]", "", temp4);temp4 <- gsub("\\.", "", temp4)
    temp5 <- df$Modifications  # modifications of the peptides
    temp6 <- df$Master.Protein.Accessions # Extract protein the peptide belongs to accession
    rownames(temp2) <- temp3
    temp2 <- cbind(temp4,temp2);temp2 <- cbind(temp5,temp2);temp2 <- cbind(temp6,temp2)

    first <- str_extract_all(names(temp2), "Sample\\.\\..*")
    second <- gsub("Sample\\.\\.", "", first)[-c(1,2,3)] # removes the names for the accession, modification and sequence
    names(temp2) <- c("Accession", "Modifications", "Sequence", second)
    temp2 <- temp2[c("Accession", "Modifications", "Sequence", myOrder)]

    temp2
}

#Sample loading correction, calculates a mean for all the samples and then scales the columns to it. 
SL_normalization <- function(df){
  
  target <- mean(colSums(df, na.rm = TRUE)) 
  norm_facs <- target / colSums(df, na.rm = TRUE)
  df_sl <- sweep(df, 2, norm_facs, FUN = "*")
  df_sl
}


# Normalizes the samples to an anritifcially generated reference standart.
IRS_normalization <- function(df, rep1=c(1:5), rep2=c(6:10), rep3=c(11:15)){
  #IRS normalization for TMT experiment with 3 replicates
  IRS <- tibble(rowSums(df[rep1], na.rm = TRUE), rowSums(df[rep2], na.rm = TRUE), rowSums(df[rep3], na.rm = TRUE)) # genrates a tibble that contains IRS for each experiment. IRS is the mean of all rows in the experiment.
  colnames(IRS) <- c("sum1", "sum2", "sum3")
  
  IRS$average <- apply(IRS, 1, function(x) exp(mean(log(x))))  #Get the geometric average intensity of each protein. Why is the exp needed????? Because log calculates the natural logarithm.
  
  # Calculate the scaling factor for each experiment by dividing the total geometric average by the sum of the experiment
  IRS$fac1 <- IRS$average / IRS$sum1
  IRS$fac2 <- IRS$average / IRS$sum2
  IRS$fac3 <- IRS$average / IRS$sum3
  # Apply the factor to the data
  df_IRS <- df 
  df_IRS[rep1] <- df_IRS[rep1] * IRS$fac1
  df_IRS[rep2] <- df_IRS[rep2] * IRS$fac1
  df_IRS[rep3] <- df_IRS[rep3] * IRS$fac1
  return(df_IRS)
}


#My PCA function
myPCA <- function(df, mylabs, myTitle){
distance.matrix <- dist(scale(t(df), center = TRUE, scale = TRUE), method = "euclidean")
mds.stuff <- cmdscale(distance.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values <- mds.stuff$points
mds.data <- data.frame(X = mds.values[,1], Y = mds.values[,2], Sample = mylabs)
mds.data$Sample <- factor(mds.data$Sample, levels = unique(mylabs))


centroids <- aggregate(cbind(X,Y)~Sample,mds.data,mean)
conf.rgn  <- do.call(rbind,lapply(unique(mds.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(mds.data[mds.data$Sample==t,1:2]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))


ggplot(data= mds.data, aes(x=X, y=Y, label = Sample, color = Sample)) +
  geom_point(size = 6, shape=18) +
  ggsci::scale_color_jama() +
  theme_bw() +
  xlab(paste("PC1 - ", mds.var.per[1], "%", sep = " ")) +
  ylab(paste("PC2 - ", mds.var.per[2], "%", sep = " ")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  geom_path(data=conf.rgn, size = 1.1, linetype = 2) +ggtitle(myTitle)
}

#calculate the mean and std of  a single extracted protein (only works for ExpG)
#TODO make the function more general
mleFit <- function(df){
    temp.fit <- apply(df, 1, function(x) coefficients(summary(lm(x ~ group + 0)))[1:8])
    temp.df <- data.frame(avg = as.numeric(temp.fit[1:4]), std = as.numeric(temp.fit[5:8]), samples = c("A", "B", "C", "D"))
    temp.df
}
#PEPTIDE FUNCTIONS. MIGHT MOVE THEN TO A SEPARE SCRIPT LATER ON 



#Check if the mapped(nonphosphorylated) peptide is detected
peptideMapped <- function(df){
    if (sum(is.na(df$Modifications)) > 0){
        TRUE
    }
    else {FALSE}
}

#Check if only one peptide is detected 
onePeptide <-function(df){
    if (dim(df)[1] < 2){
        TRUE
    }
    else {FALSE}
}

#merges peptides that have different TMT labelling
mergePeptides <- function(df){
    #Modifications and peptide sequnce must be in "Modifications" and "Sequence" column respectively.
    temp <- df
    p1 <- "\\dxPhospho.*"
    p2 <- "(\\[.*\\])" 
    p3 <- "[STY]\\d+|S/Y|S/T|T/S|T/Y|Y/S|Y/T|S/T/Y|S/Y/T|T/S/Y|T/Y/S|Y/S/T|Y/T/S|S|T|Y" # very nice pattern 

    m1 <-  str_extract(df$Modifications,p1) # extract the phospho modification part
    m2 <- str_extract(m1, p2) # extract the modification residue in the brackets
    m3 <- str_extract_all(m2, p3) # extract the residue with number or just the residues for not localized modifications
    m4 <- rapply(m3, function(x) paste(x, collapse = ";")) # merge the modifications into one strings

    temp$mod <- m4 # adding the extracting modification to df
    temp <- temp %>% group_by(Sequence, mod) %>% mutate_at(vars(colnames(temp[4:18])), sum) # sum the rows that have the same peptide sequence and modification
    temp <- temp[!duplicated(temp[,c("Sequence", "mod")]), ] # remove the duplicate rows left behind by the previous line
    temp
}


#convert the peptide abundance to occupancy
calcOccupancy <- function(df) {
    df <- mergePeptides(df) #merges peptides with different TMT modifications
    total <- apply(df[4:11],2,sum); total <- t(replicate(dim(df)[1], total)) #calculte sum and replicate in into a matrix for further calculation
    occ <- 1 - (total - df[4:11])/total # Occupancy calculation
    occ <- cbind(df[,1:3], occ)
    occ$mapped <- peptideMapped(df)
    occ$one <- onePeptide(df)
    occ
}

#Extract the  phospho modification on the peptide
extractSite <- function(df){
    #Modifications and peptide sequnce must be in "Modifications" and "Sequence" column respectively.
    temp <- df
    p1 <- "(\\[.*?\\])" # pattern one for mathcing modification. Extract all the string between []
    p2 <- "[STY]\\d*" # pattern two for matching modification. Extracts the aminoacid and number from the first match e.g. S18

    m1 <- str_extract(df$Modifications, p1) # first match
    m2 <- str_extract_all(m1, p2) # second match
    temp$modSite <- paste(temp$Accession, m2, sep = "-")
    temp
}


# Functions for changing the modification text from PD into a nicer text format
mergeSites <- function(a,b,c){
  #Works for mergin three modifications max
  if (!(is.na(a))){
    if (!(is.na(c))){ site <- paste(a,b,c, sep = ";")}
    else if ((!is.na(b))){ site <- paste(a,b, sep = ";")}
    else {site = a}
  }
  else { site <- "n.m."}
  site
}


extMasterMod <- function(df){
  #Modification column has to be called Mod
  modLoc1 <- str_extract(df$Modification.Master, pattern = "(\\[.*\\])") #modifications with brackets
  modLoc2 <- gsub("\\].*","", modLoc1) # remove first bracket
  modLoc3 <- gsub("\\[","", modLoc2) # remove second bracket
  df$AllMods <- modLoc3 # Adding the modification sites

  temp <- df %>% separate(AllMods, c("a", "b", "c"), sep = ";") # separating sites into columns by ;
  a <- str_extract(temp$a, ".[0-9]+") # Extracting the modification from each column, not mapped values will be NA
  b <- str_extract(temp$b, ".[0-9]+")
  c <- str_extract(temp$c, ".[0-9]+")

  df$MasterMod <- mapply(mergeSites, a, b,c)
  df
}


extMasterMod2 <- function(df){
  #Modification column has to be called Mod
  modLoc1 <- str_extract(df$ModificationMaster, pattern = "(\\[.*\\])") #modifications with brackets
  modLoc2 <- gsub("\\].*","", modLoc1) # remove first bracket
  modLoc3 <- gsub("\\[","", modLoc2) # remove second bracket
  df$AllMods <- modLoc3 # Adding the modification sites

  temp <- df %>% separate(AllMods, c("a", "b", "c"), sep = ";") # separating sites into columns by ;
  tofind <- ".[0-9]+|T/S/Y|S/T|S/Y|Y/S|S|Y|T|T/Y|Y/T"
  a <- str_extract(temp$a, tofind) # Extracting the modification from each column, not mapped values will be NA
  b <- str_extract(temp$b, tofind)
  c <- str_extract(temp$c, tofind)

  df$MasterMod <- mapply(mergeSites, a, b,c)
  df
}


#Function for writing fasta files

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"ID"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Motifs"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


#Extract motif function. Used with mapply
extractAA <- function(modLoc,protSeq,l){
  #Check if position is mapped
  if (is.integer(modLoc)){
    loc <- as.integer(modLoc)
    motif <- substr(protSeq, loc-l, loc+l)
  }
  else {
    motif <- "NotMapped"
  }
  motif
}
makePad <- function(n,s){ 
  #Function for making the padding strings
  if (n == TRUE) {
    paste(rep("X",s), collapse = "")
  } else  {""} # adds empty string
}
padMotif <- function(modLoc, protSeq, motif,l){
    #Function for applying the pading to motifs
    temp <- data.frame(loc = modLoc, protSeq = protSeq, motif = motif)
    temp$x1 <- nchar(temp$protSeq) - (temp$loc+l) #calculate missing residues on left or right side
    temp$x2 <- temp$loc - (l+1)
    
    #Checking if padding is needed
    temp$p1 <- temp$x1 < 0
    temp$p2 <- temp$x2 < 0
    #Making padding strings
    temp$pad1 <- mapply(makePad, temp$p1, abs(temp$x1))
    temp$pad2 <- mapply(makePad, temp$p2, abs(temp$x2))
    #Add the padding
    temp$motif <- paste(temp$motif, temp$pad1, sep = "")
    temp$motif <- paste(temp$pad2, temp$motif, sep = "")
    temp$motif
}
extractMotif <- function(df, l){
  #motif extraction function 
  #The modification location (integer) has to be in column modLoc
  #The protein sequence in the AAseq column
  # l - the motif size
  #extracting non-padded motifs
  #Modification have to be in MasterMod colunmn
  temp <- df
  temp$modLoc <- as.integer(str_extract(df$MasterMod, pattern = "[:digit:]+"))
  temp$motif <- extractAA(temp$modLoc, temp$AAseq, l)
  #padding the motif
  motifs <- padMotif(temp$modLoc, temp$AAseq, temp$motif, l)
  motifs
}


mleFit2 <- function(df){
    temp.fit1 <- apply(df, 1, function(x) coefficients(summary(lm(x ~ group + 0)))[1:2])
    rownames(temp.fit1) <- paste(unique(group), "avg", sep = "")

    temp.fit2 <- apply(df, 1, function(x) coefficients(summary(lm(x ~ group + 0)))[3:4])
    rownames(temp.fit2) <- paste(unique(group), "sd", sep = "")

    temp.fit <- cbind(t(temp.fit1), t(temp.fit2))
    temp.fit
}