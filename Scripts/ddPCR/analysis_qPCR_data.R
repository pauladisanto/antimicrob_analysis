process_df <- function(data) {
  # Calculate the mean and standard dispersion of "data" and return them as a list with two arrays

  # Numbers useful to build the new array (the last two are to be discarded)
  step = 3
  n_rep = ( length(data) - 2 )

  # Create an array with the new averages. Initialize as a subarray of the original one to have the right names
  # Calculate the indices for subarray
  idx <- c(seq(1, n_rep, by = 2*step), 
           seq(2, n_rep, by = 2*step), 
           seq(3, n_rep, by = 2*step))

  idx <- sort(idx) # Sort idx in increasing order
  means <- data[idx]

  for(i in 1:n_rep/2) {
    means[i] = mean( c(data[idx[i]],data[idx[i]+3]), na.rm = TRUE )
  }

  # Calculate the third mean, among the technical replicates
  idx2 <- c(seq(1, n_rep/2, by = 3))
  mean_2 <- means[idx2]
  sigma_2 <- means[idx2]
  for(i in 1:n_rep/6) {
    subarray = c(means[idx2[i]], means[idx2[i]+1], means[idx2[i]+2])
    mean_2[i] = mean(subarray)
    sigma_2[i] = sd(subarray)
  }
  
  # Create a list to store the output and return it
  output <- list(mean_2, sigma_2)
  
  # Return the list containing the arrays
  return(output)
}

#=============================================================================

# Analysis of the qPCR data
library(ggplot2)
library(tidyverse)
library("agricolae")


#filename = "DATA_1.csv"
filename = "DATA_2.csv"

df = read.table(filename, header=TRUE, sep="\t")
df[df == "0"] <- "NA"
df$Concentration=  as.numeric(df$Concentration)

# Define a new dataframe in which we apply log_2 to the Concentration column
#df2 <- df
df$Concentration <- log(df$Concentration, base = 2)

# Names of all the genes to loop over
gene_names = c("gyr_reference", "agrA", "hld", "spa", "icaA")

# Calculate mean and sd for each gene
mean = list() 
sd = list()
raw_mean = list()
for(k in 1:length(gene_names)) {

  print(paste("Results for", gene_names[k]))

  # Change df for df2????
  table <- subset(df, Gene==gene_names[k], select=c(Gene, Sample, Concentration))
  raw_mean[[k]] = tapply(table$Concentration, factor(table$Sample), mean, na.rm=TRUE)

  result <- process_df(raw_mean[[k]])

  # Access the output arrays from the returned list
  mean[[k]] <- result[[1]]
  sd[[k]] <- result[[2]]
  # Print the output arrays
  print(mean[[k]])
  print(sd[[k]])
  print(" ")
}

#------------------------------------------------------------
#Multiple comparisons
#raw_mean[[]]
test_agr = na.omit(tibble::rownames_to_column(as.data.frame(raw_mean[[2]])))  
test_hld = na.omit(tibble::rownames_to_column(as.data.frame(raw_mean[[3]])))  
test_spa = na.omit(tibble::rownames_to_column(as.data.frame(raw_mean[[4]])))  
test_icaA = na.omit(tibble::rownames_to_column(as.data.frame(raw_mean[[5]]))) 

tests = list(test_agr, test_hld, test_spa, test_icaA)
for(test in tests){
  colnames(test) = c("Treatment","fold_change")
  test = filter(test, Treatment != 'ZZPC')
  test$Treatment <- gsub("^AIP.*", "AIP", test$Treatment)
  test$Treatment <- gsub("^SUP.*", "SUP", test$Treatment)
  test$Treatment <- gsub("^ACN.*", "ACN", test$Treatment)
  test$Treatment <- gsub("^MED.*", "MED", test$Treatment)
  test$Treatment <- gsub("^N10.*", "N10", test$Treatment)
  test$Treatment <- gsub("^N1_.*", "N1", test$Treatment)
  test$Treatment <- gsub("^N0.1.*", "N0.1", test$Treatment)
  test$Treatment <- gsub("^N0.01.*", "N0.01", test$Treatment)

  aovScreening<-aov(fold_change~Treatment,data=test)
  summary(aovScreening)

  postestScreening<-LSD.test(aovScreening, "Treatment")
  print(postestScreening)
}


#===============================================

#Delta calculation of the means

agrA_gyrA=mean[[2]]-mean[[1]]
hld_gyrA=mean[[3]]-mean[[1]]
spa_gyrA=mean[[4]]-mean[[1]]
icaA_gyrA=mean[[5]]-mean[[1]]

#Substraction of the  delta media value and reorder of the rows
target = c("AIP_A_BR_1", "SUP_A_BR_1", "ACN_A_BR_1", "MED_A_BR_1", "N10_A_BR_1",  "N1_A_BR_1", "N0.1_A_BR_1", "N0.01_A_BR_1")

N_agrA_gyrA = tibble::rownames_to_column(as.data.frame(agrA_gyrA - agrA_gyrA[3])) %>% arrange(factor(rowname, levels = target))
N_hld_gyrA = tibble::rownames_to_column(as.data.frame(hld_gyrA - hld_gyrA[3])) %>% arrange(factor(rowname, levels = target))
N_spa_gyrA = tibble::rownames_to_column(as.data.frame(spa_gyrA - spa_gyrA[3])) %>% arrange(factor(rowname, levels = target))
N_icaA_gyrA = tibble::rownames_to_column(as.data.frame(icaA_gyrA - icaA_gyrA[3])) %>% arrange(factor(rowname, levels = target))

#Calculation of the SD for the deltas reorder of the rows
#SD = sqrt(a^2 + a^2)
#apply potency two to SD list
sd2 = lapply(sd,function(x) x^2)
#Sqrt of the sum
SD_agrA = tibble::rownames_to_column(as.data.frame(sqrt(sd2[[1]]+sd2[[2]]))) %>% arrange(factor(rowname, levels = target))
SD_hld = tibble::rownames_to_column(as.data.frame(sqrt(sd2[[1]]+sd2[[3]]))) %>% arrange(factor(rowname, levels = target))
SD_spa = tibble::rownames_to_column(as.data.frame(sqrt(sd2[[1]]+sd2[[4]]))) %>% arrange(factor(rowname, levels = target))
SD_icaA = tibble::rownames_to_column(as.data.frame(sqrt(sd2[[1]]+sd2[[5]])))%>% arrange(factor(rowname, levels = target))


#merge means and SD values

agrA_dataset = N_agrA_gyrA %>% 
  left_join(SD_agrA, by=c("rowname"))

hld_dataset = N_hld_gyrA %>% 
  left_join(SD_hld, by=c("rowname"))

spa_dataset = N_spa_gyrA %>% 
  left_join(SD_spa, by=c("rowname"))

icaA_dataset = N_icaA_gyrA %>% 
  left_join(SD_icaA, by=c("rowname"))


colnames(agrA_dataset) = c("Treatment","fold_change", "sd")
colnames(hld_dataset) = c("Treatment","fold_change", "sd")
colnames(spa_dataset) = c("Treatment","fold_change", "sd")
colnames(icaA_dataset) = c("Treatment","fold_change", "sd")


#===================================================================
#Plots
vector_name = c("agrA", "hld", "spa", "icaA")
vector = list(agrA_dataset, hld_dataset, spa_dataset, icaA_dataset)

#significances_agr = c("c", "a", "ab", "b", "ab", "ab", "ab", "b")
#significances_hld = c("d", "a", "b", "ab", "c", "ab", "bc", "bc")
#significances_spa = c("a", "b", "b", "b", "a", "b", "b", "b")
#significances_icaA = c("ab", "ab", "c", "abc", "a", "bc", "abc", "ab")
significances_agr = c("a", "d", "c", "b", "e", "bc", "a", "a")
significances_hld = c("a", "d", "b", "b", "c", "b", "a", "a")
significances_spa = c("ab", "ab", "ab", "ab", "b", "a", "a", "ab")
significances_icaA = c("abc", "bc", "bc", "bc", "c", "a", "bc", "ab")

significances = list(significances_agr, significances_hld, significances_spa, significances_icaA)

distance = c(0.6, 0.6, 0.6, 2.4, 0.6, 0.6, 0.6, 0.6)
xc = c(1,2,3,4,5,6,7,8)
textSize=12

for(i in 1:length(vector)) {
  Mean_i = vector[[i]]$fold_change
  SD_i = vector[[i]]$sd

  x_order = c("AIP_A_BR_1" , "SUP_A_BR_1" , "ACN_A_BR_1" , "MED_A_BR_1" , "N10_A_BR_1", "N1_A_BR_1" , "N0.1_A_BR_1", "N0.01_A_BR_1") 
  vector[[i]]$Treatment <- factor(vector[[i]]$Treatment, levels = x_order)
  p<-ggplot(data= vector[[i]], aes(x=Treatment, y=fold_change,fill = Treatment)) +
    geom_bar(stat="identity") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    scale_fill_discrete(labels=c("AIP-I", "agrIII SUP", "ACN", "TSB", "10mM", "1mM", "0.1mM","0.01mM")) +
    theme(legend.position="top", legend.title=element_text(size=textSize), legend.text=element_text(size=textSize),
      axis.text=element_text(size=textSize)) +
    #labs(title = paste0("S. aureus ATTC 25923-", vector_name[i])) +
    labs(title = paste0("S. aureus 8325-4-" , vector_name[i])) +
    ylim(-10, 9) +
    #annotate("text")+
    geom_errorbar(aes(ymin=fold_change-sd, ymax=fold_change+sd), width=0.2)+
    annotate("text", x = xc, y = Mean_i+sign(Mean_i)*(SD_i+distance)+(1-abs(sign(Mean_i)))*distance, label = significances[[i]], size = 7, family = "Times", fontface = "bold.italic", colour = "black")
    #+ggtitle(paste0("Plot ", i)) # Add title
  print(p)
  ggsave(filename = paste0("plot_new_", i, ".png"), plot = p)
}  

#, axis.text.x=element_text(size=textSize), axis.text.y=element_text(size=textSize)

#target = c("AIP-I", "agrIII SUP", "ACN", "TSB", "10mM", "1mM", "0.1mM","0.01mM")
#labels=c("ACN", "AIP-I", "TSB", "0.01mM", "0.1mM", "1mM", "10mM", "agrIII SUP")


#===========================


#Comparisons among the samples t-test

n = nrow(agrA_dataset)
tagr <- matrix(nrow=n, ncol=n)  
for(i in 1:(n-1)) {
  for(j in (i+1):n) {
    tagr[i,j] = abs(agrA_dataset[i,2] - agrA_dataset[j,2]) / sqrt(agrA_dataset[i,3]^2 + agrA_dataset[j,3]^2)
  }
}
