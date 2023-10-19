
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(agricolae)
read_table <- function(filename) {  
	# Read the raw data from the file, skipping the first 14 rows which contain irrelevant information
	raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")
}

# define the files to analyse

filename1 = "MIC_b_S_epidermidis_experimento_1_compuesto_C1.csv"
filename2 = "MIC_b_S_epidermidis_experimento_1_compuesto_C1_replica_2.csv"
filename3 = "MIC_b_S_epidermidis_experimento_1_compuesto_C1_replica_3.csv"
 


# read the tables
raw_data1 = read_table(filename1)
raw_data2 = read_table(filename2)
raw_data3 = read_table(filename3)

#=============================================================================================================
# List of dataframes (replace with your actual list of dataframes)
file_list <- list(raw_data1, raw_data2, raw_data3)

# Define a function to process each dataframe and return a dataframe
process_dataframe <- function(df) {
  # Set column names
  column_names <- c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
  colnames(df) <- column_names

  # Remove the 'X' column
  data <- subset(df, select = -c(X))

  # Create a new dataframe for compounds
  compounds <- data.frame()
  sum <- 0

  for (i in 1:n_row) {
    for (j in 1:n_col) {
      sum <- sum + 1

      # Obtain the index "i_c" for each compound
      i_c <- as.integer((sum - 1) / n_m) + 1
      j_c <- (j - 1) %% 12 + 1
      name_c <- paste('Condition', sprintf("%02d", i_c), sep = "_")
      compounds[n_m * (i_c - 1) + j_c, 1] <- name_c
      compounds[n_m * (i_c - 1) + j_c, 2] <- data[i, j]

      # Obtain the index "i_b" for each bacteria
      i_b <- as.integer((sum - 1) / n_b) + 1
      name_b <- paste('bacteria', as.character(i_b), sep = "_")
      compounds[n_m * (i_c - 1) + j_c, 3] <- name_b
    }
  }

  # Add a header to compounds dataframe
  column_names <- c('Condition', 'Absorbance', 'Bacteria')
  colnames(compounds) <- column_names

  # Return the processed dataframe
  return(compounds)
}

# Initialize a list to store processed dataframes
processed_data_list <- list()

# Define the number of compounds, rows, and columns
n_m <- 12
n_b <- 96
n_row <- 8
n_col <- 12

# Loop through each dataframe and process it, adding the result to the list
for (i in 1:length(file_list)) {
  processed_data <- process_dataframe(file_list[[i]])
  processed_data_list[[i]] <- processed_data
}

# Print or return the list of processed dataframes
processed_data_list

#=============================================================================================================


bacteria_condition =  c("Replicate_1","Replicate_2","Replicate_3","Negative_control","Replicate_4","Replicate_5","Replicate_6","Positive_control")

processed_data_list_1 = processed_data_list[[1]]
processed_data_list_2 = processed_data_list[[2]]
processed_data_list_3 = processed_data_list[[3]]

# Define a function to replace values in a dataframe
replace_values <- function(dataframe, conditions, replacements) {
  for (i in 1:length(conditions)) {
    condition <- conditions[i]
    replacement <- replacements[i]
    dataframe[dataframe == condition] <- replacement
  }
  return(dataframe)
}

# Define conditions and replacements
conditions <- paste('Condition', sprintf("%02d", 1:8), sep = "_")
replacements <- bacteria_condition[1:8]

# Process and replace values in each dataframe
processed_data_list_1 <- replace_values(processed_data_list_1, conditions, replacements)
processed_data_list_2 <- replace_values(processed_data_list_2, conditions, replacements)
processed_data_list_3 <- replace_values(processed_data_list_3, conditions, replacements)


concentration = c(960) 
for(i in 2:12) {concentration[i] = concentration[i-1]/2}
#================================================================================
processed_data_list_1$Bacteria <- concentration
processed_data_list_2$Bacteria <- concentration
processed_data_list_3$Bacteria <- concentration

colnames(processed_data_list_1)= c("Condition", "Absorbance",  "concentration")
colnames(processed_data_list_2)= c("Condition", "Absorbance",  "concentration")
colnames(processed_data_list_3)= c("Condition", "Absorbance",  "concentration")

#================================================================================
# Define a function to calculate the mean for the negative control condition
calculate_mean_negative_control <- function(dataframe) {
  # Subset the dataframe for the "Negative_control" condition
  Tabla_Ensayo_control_negativo <- subset(dataframe, Condition == "Negative_control", select = c(Condition, concentration, Absorbance))
  
  # Calculate the mean for the "Negative_control" condition
  mean_negative_control <- tapply(Tabla_Ensayo_control_negativo$Absorbance, factor(Tabla_Ensayo_control_negativo$Condition), mean,  na.rm = TRUE)
  
  return(mean_negative_control)
}

# Apply the function to each dataframe
mean_negative_control_1 <- calculate_mean_negative_control(processed_data_list_1)
mean_negative_control_2 <- calculate_mean_negative_control(processed_data_list_2)
mean_negative_control_3 <- calculate_mean_negative_control(processed_data_list_3)

# Print or access the means for each dataframe
print(mean_negative_control_1)
print(mean_negative_control_2)
print(mean_negative_control_3)

#================================================================================

# Create a list of your processed_data_list dataframes
dataframe_list <- list(processed_data_list_1, processed_data_list_2, processed_data_list_3)

# Define a function to process a dataframe
process_dataframe <- function(dataframe) {
  # Subset the dataframe for the specified conditions
  Tabla_Ensayo <- subset(dataframe, Condition %in% c("Replicate_1", "Replicate_2", "Replicate_3",  "Replicate_4",  "Replicate_5",  "Replicate_6" ,"Positive_control"), select = c(concentration, Absorbance, Condition))
  
  # Replace the values of the concentration for the positive control
  Tabla_Ensayo$concentration[Tabla_Ensayo$Condition == "Positive_control"] <- 0

  # Filter rows
  Tabla_Ensayo <- Tabla_Ensayo[c(1:72, 77:82), ]
  
  return(Tabla_Ensayo)
}

# Use lapply to apply the process_dataframe function to all dataframes
processed_data_list_processed <- lapply(dataframe_list, process_dataframe)

# Access or use the processed dataframes as needed

#============================================================================
#Extract the negative OD value

processed_data_list_processed[[1]]$Absorbance = processed_data_list_processed[[1]]$Absorbance-mean_negative_control_1
processed_data_list_processed[[2]]$Absorbance = processed_data_list_processed[[2]]$Absorbance-mean_negative_control_2
processed_data_list_processed[[3]]$Absorbance = processed_data_list_processed[[3]]$Absorbance-mean_negative_control_3

#============================================================================
#Calculate the means and SDs
mean_absorbance_1 <- tapply(processed_data_list_processed[[1]]$Absorbance,factor(processed_data_list_processed[[1]]$concentration),mean)
mean_absorbance_2 <- tapply(processed_data_list_processed[[2]]$Absorbance,factor(processed_data_list_processed[[2]]$concentration),mean)
mean_absorbance_3 <- tapply(processed_data_list_processed[[3]]$Absorbance,factor(processed_data_list_processed[[3]]$concentration),mean)

#============================================================================
#Generate a new dataframe with the means
mean_MIC = list(mean_absorbance_1, mean_absorbance_2, mean_absorbance_3)

df_list <- mean_MIC

# Combine the dataframes into a single dataframe
combined_df <- do.call(rbind, df_list)

# Print the combined dataframe
print(combined_df)

#traspose the dataframe with the means
transposed_mean_MBEC <- as.matrix(t(combined_df))


# Convert the matrix into a dataframe
df <- as.data.frame(transposed_mean_MBEC)

#I include the concentration number as a column 

setDT(df, keep.rownames = TRUE)[]

# Set the column names
colnames(df) <- c("Concentration", "Column1", "Column2", "Column3")


# Melt the dataframe
melted_df <- melt(df, id.vars = "Concentration")

melted_df = melted_df[,c(1,3)]

# Print the melted dataframe
print(melted_df)

colnames(melted_df) <- c("Concentration", "Mean")

melted_df$Concentration= as.numeric(melted_df$Concentration)

#============================================================================
#Test estadistico
aovScreening<-aov(Mean~Concentration,data=melted_df)
summary(aovScreening)

postestScreening<-LSD.test(aovScreening, "Concentration")
print(postestScreening)


#=====================================================================================
#Means
  mean_MIC_b <- tapply(melted_df$Mean, factor(melted_df$Concentration), mean)

#SDs 
  sd_MIC_b <- tapply(melted_df$Mean, factor(melted_df$Concentration), sd)

#====================================================================
#plot

Concentration = c(960) 
for(i in 2:12) {Concentration[i] = Concentration[i-1]/2}

Concentration = c(Concentration , 0)



df_plot = data.frame(
    "Conc" = rev(Concentration),
    "mean" = mean_MIC_b, 
    "SD" = sd_MIC_b
    )

#====================================================================
#obliga a ordenar a que el eje x se ordene de manera equidistante
df_plot$Conc <- factor(df_plot$Conc, levels = df_plot$Conc)

title = "MIC (pre MBEC) S. epidermidis C1 3 Replicates"
#significances = c("a", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "bc", "ab", "ab", "bc", "c")
significances = c("", "", "", "", "", "", "", "", "", "", "", "", "")

xc = c(1,2,3,4,5,6,7,8,9,10,11,12,13)
distance = c(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)
textSize=8
Mean_i = df_plot$mean
SD_i = df_plot$SD

p <- ggplot(data=df_plot, aes(x=Conc, y=mean)) + 
          theme_bw() +
          geom_bar(stat="identity",  width=0.9, fill = "red")+ 
          ylim(-0.08, 1)+
          geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.2, position=position_dodge(0.05)) +  
              theme(axis.text.x = element_text(color="black", size=12, angle=10)) +
              theme(legend.position="top",legend.title = element_text(size = textSize), legend.text  = element_text(size = textSize)) + 
          labs(title=element_blank(), x=expression(paste("Compound concentration (", mu, "M)")), y =  "Planctonik growth after 24 hs (OD_600)")+
          annotate("text", x = xc, y = Mean_i+SD_i+distance, label = significances, size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
                                                                             ggtitle(title)


