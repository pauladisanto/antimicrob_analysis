library(ggplot2)
library(reshape2)
library(data.table)
library("agricolae")
library(DescTools)


#generates a list with the three dataframes           
filenames <- c("CFU_counting_experimento_1_S_epidermididis_only_Control_BT.csv", "CFU_counting_experimento_2_S_epidermididis_only_Control_BT.csv", "CFU_counting_experimento_3_S_epidermididis_only_Control_BT.csv")

# Create an empty list to store the data frames
data_frames <- list()

# Loop through the filenames and read the data into data frames
for (filename in filenames) {
  data <- read.table(filename, header = TRUE, sep = ",")
  data_frames[[filename]] <- data
}

#======================================================================================================================

# List of dataframes
dataframe_list <- data_frames  

# Create an empty list to store the results
result_list <- list()

# Loop through each dataframe in the list
for (df in dataframe_list) {
  new_data <- data.frame()

  # Calculate the means for each pair of rows
  for (i in 1:(nrow(df) / 2)) {
    pair <- df[(2 * i - 1):(2 * i), ]  # Get a pair of consecutive rows
    pair_mean <- mean(pair$CFU.ml)  # Calculate the mean for the pair
    new_row <- data.frame(
      Concentration_number = pair$Concentration_number[1],
      Concentration = pair$Concentration[1],
      Mean_CFU_ml = pair_mean
    )
    new_data <- rbind(new_data, new_row)
  }

  # Append the results for this dataframe to the result_list
  result_list <- append(result_list, list(new_data))
}

# Print the list of dataframes
print(result_list)

#======================================================================================
#replace zeros as NA values 

# Loop through each dataframe in the list
for (i in 1:length(result_list)) {
  # Access the dataframe within the list using double brackets
  df <- result_list[[i]]
  
  # Replace 0 with NA in the "ValueColumn"
  df$Mean_CFU_ml[df$Mean_CFU_ml == 0] <- NA
  
  # Update the dataframe within the list
  result_list[[i]] <- df
}

# Print the updated list of dataframes
for (i in 1:length(result_list)) {
  print(paste("Dataframe", i))
  print(result_list[[i]])
}

#======================================================================================
#Apply log10 to the means in a list

df_list <- result_list
# Column name to apply log10
column_name <- "Mean_CFU_ml"

# Apply log10 to the specified column in each dataframe
df_list <- lapply(df_list, function(df) {
  df[[column_name]] <- log10(df[[column_name]])
  return(df)
})

# Print the updated list of dataframes
for (i in 1:length(df_list)) {
  print(paste("Dataframe", i))
  print(df_list[[i]])
}

#======================================================================================
#replace NA as zero values 

# Loop through each dataframe in the list
for (i in 1:length(df_list)) {
  # Access the dataframe within the list using double brackets
  df <- df_list[[i]]
  
  # Replace NA with 0 in the "Mean_CFU_ml"
  df$Mean_CFU_ml[is.na(df$Mean_CFU_ml)] <- 0
  
  # Update the dataframe within the list
  df_list[[i]] <- df
}

# Print the updated list of dataframes
for (i in 1:length(df_list)) {
  print(paste("Dataframe", i))
  print(df_list[[i]])
}

#======================================================================================

#mean calculations

mean_CFU <- list()

# Loop through the data frames
for (i in 1:length(df_list)) {
  mean_CFU[[i]] <- tapply(df_list[[i]]$Mean_CFU_ml, factor(df_list[[i]]$Concentration_number), na.rm = TRUE, mean)
}

print(mean_CFU)
#====================================================================

df_list <- mean_CFU

# Combine the dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)

# Print the combined dataframe
print(combined_df)

#traspose the dataframe with the means
transposed_mean_CFU <- as.matrix(t(combined_df))


# Convert the matrix into a dataframe
df <- as.data.frame(transposed_mean_CFU)

#I include the concentration number as a column 

setDT(df, keep.rownames = TRUE)[]

# Set the column names
colnames(df) <- c("Concentration", "Column1", "Column2", "Column3")


# Melt the dataframe
melted_df <- melt(df, id.vars = "Concentration")

# Print the melted dataframe
print(melted_df)

melted_df = melted_df[,c(1,3)]

colnames(melted_df) <- c("Concentration", "Log10_mean")

melted_df$Concentration= as.numeric(melted_df$Concentration)

#add a column with treatment and control

treatment =  rev(c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C0"))

melted_df = data.frame(melted_df, Treatment = treatment)

#====================================================================
aovScreening<-aov(Log10_mean~Treatment,data=melted_df)
summary(aovScreening)

#perform Dunnett's Test
DunnettTest(x=melted_df$Log10_mean, g=melted_df$Treatment)

#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#====================================================================
#agricolae test

#aovScreening<-aov(Log10_mean~Concentration,data=melted_df)
#summary(aovScreening)

#postestScreening<-LSD.test(aovScreening, "Concentration")
#print(postestScreening)

#====================================================================
#Means
  mean_CFU <- tapply(melted_df$Log10_mean, factor(melted_df$Concentration), na.rm = TRUE, mean)

#SDs 
  sd_CFU <- tapply(melted_df$Log10_mean, factor(melted_df$Concentration), na.rm = TRUE, sd)

#====================================================================
#plot

Concentration = c(960) 
for(i in 2:12) {Concentration[i] = Concentration[i-1]/2}

Concentration = c(Concentration , 0)



df_plot = data.frame(
    "Conc" = rev(Concentration),
    "Log10_mean" = mean_CFU, 
    "SD" = sd_CFU
    )


#obliga a ordenar a que el eje x se ordene de manera equidistante
df_plot$Conc <- factor(df_plot$Conc, levels = df_plot$Conc)

title = "CFU S. epidermidis C1 Replicas 1, 2 and 3"
#significances = c("abc", "a", "a", "a", "a", "ab", "bcd", "d", "cd", "de", "d", "de", "e")
significances = c("", "", "", "", "", "*", "**", "*", "**", "***", "***", "***", "***")
xc = c(1,2,3,4,5,6,7,8,9,10,11,12,13)
distance = c(0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16)
textSize=6
Mean_i = df_plot$Log10_mean
SD_i = df_plot$SD

p <- ggplot(data=df_plot, aes(x=Conc, y=Log10_mean)) + 
          theme_bw() +
          geom_bar(stat="identity",  width=0.9, fill = "coral2")+ 
          ylim(0, 7)+
          geom_errorbar(aes(ymin=Log10_mean-SD, ymax=Log10_mean+SD), width=0.2, position=position_dodge(0.05)) +  
              theme(axis.text.x = element_text(color="black", size=10, angle=10)) +
              theme(legend.position="top",legend.title = element_text(size = textSize), legend.text  = element_text(size = textSize)) + 
          labs(title=element_blank(), x=expression(paste("Compound concentration (", mu, "M)")), y = "Log 10 (Mean)")+
          annotate("text", x = xc, y = Mean_i+SD_i+distance, label = significances, size = 3, family = "Times", fontface = "bold.italic", colour = "black")+
                                                                             ggtitle(title)

ggsave(filename = "CFU_S_epidermidis.pdf", plot = p, width = 7, height = 4, units = "in", dpi = 300)

ggsave(filename = "CFU_S_epidermidis.png", plot = p, width = 7, height = 4, units = "in", dpi = 300)