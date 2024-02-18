
library(dplyr)
library(ggplot2)
library(reshape2)
library(data.table)
library(agricolae)
library(DescTools)
library(gridExtra)

read_table <- function(filename) {  
	# Read the raw data from the file, skipping the first 14 rows which contain irrelevant information
	raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")
}

# define the files to analyse

filename1 = "Interaction_S.epidermidis_C1_7_5_plus_ATB_AFTER_recovery_media_replica_1_MBEC.csv"
filename2 = "Interaction_S.epidermidis_C1_7_5_plus_ATB_AFTER_recovery_media_replica_2_MBEC.csv"


# read the tables
raw_data1 = read_table(filename1)
raw_data2 = read_table(filename2)

#============================================================================================================

  # Dummy column that is later removed
  column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')

  # Define the column names
  colnames(raw_data1) = column_names
  colnames(raw_data2) = column_names
  raw_data1 <- subset(raw_data1, select = -c(X))
  raw_data2 <- subset(raw_data2, select = -c(X))

 
#Controls
concentration0 = c(0) 
for(i in 2:8) {concentration0[i] = concentration0[i-1]/2}

#RIF and LEVO
concentration1 = c(128) 
for(i in 2:14) {concentration1[i] = concentration1[i-1]/2}

#FUS
concentration2 = c(256) 
for(i in 2:10) {concentration2[i] = concentration2[i-1]/2}

#CLI
concentration3 = c(256) 
for(i in 2:12) {concentration3[i] = concentration3[i-1]/2}

#LZD
concentration4 = c(128) 
for(i in 2:7) {concentration4[i] = concentration4[i-1]/2}

#SXT
concentration5 = c(64) 
for(i in 2:7) {concentration5[i] = concentration5[i-1]/2}

#OXA
concentration6 = c(64) 
for(i in 2:8) {concentration6[i] = concentration6[i-1]/2}

#VAN
concentration7 = c(128) 
for(i in 2:9) {concentration7[i] = concentration7[i-1]/2}

#FOX
concentration8 = c(64) 
for(i in 2:7) {concentration8[i] = concentration8[i-1]/2}
concentration = c(rev(concentration0),rev(concentration1), rev(concentration1), rev(concentration2), rev(concentration3), rev(concentration4), rev(concentration5), rev(concentration6), rev(concentration7),rev(concentration8))


CN = c("Negative_control")
CN <- rep(CN, each = 2)
CP = c("Positive_control")
CP <- rep(CP, each = 6)
RIF = c("RIF")
RIF = rep(RIF, each = 14)
LEVO = c("LEVO")
LEVO = rep(LEVO, each = 14)
FUS = c("FUS")
FUS = rep(FUS, each = 10)
CLI= c("CLI")
CLI = rep(CLI, each = 12)
LZD= c("LZD")
LZD = rep(LZD, each = 7)
SXT= c("SXT")
SXT = rep(SXT, each = 7)
OXA = c("OXA")
OXA = rep(OXA, each = 8)
VAN = c("VAN")
VAN = rep(VAN, each = 9)
FOX = c("FOX")
FOX = rep(FOX, each = 7)


antibiotic = c(CN, CP, RIF, LEVO, FUS, CLI, LZD,  SXT, OXA, VAN,  FOX)


create_compounds_dataframe <- function(raw_data) {
  n_row <- nrow(raw_data)
  n_col <- ncol(raw_data)
  n_m <- 3
  
  compounds <- data.frame()
  
  for (j in 1:n_col) {
    for (i in 1:n_row) {
      index <- (j - 1) * n_row + i
      compounds[index, 1] <- raw_data[i, j]
    }
  }
  
  return(compounds)
}

compounds1 <- create_compounds_dataframe(raw_data1)

compounds2 <- create_compounds_dataframe(raw_data2)

 column_names = c("Absorbance_600")

colnames(compounds1) = column_names
colnames(compounds2) = column_names


 #Add the column with the concentrations
 compounds1$Concentration <- concentration
 compounds2$Concentration <- concentration

 #Add the column with the antibiotics name
 compounds1$antibiotic <- antibiotic
 compounds2$antibiotic <- antibiotic

#=============================================================================================================
  
# Define a function to calculate the mean for the negative control condition
calculate_mean_negative_control <- function(dataframe) {
  # Subset the dataframe for the "Negative_control" condition
  Tabla_Ensayo_control_negativo <- subset(dataframe, antibiotic == "Negative_control", select = c(Absorbance_600, Concentration, antibiotic))
  
  # Calculate the mean for the "Negative_control" condition
  mean_negative_control <- tapply(Tabla_Ensayo_control_negativo$Absorbance_600, factor(Tabla_Ensayo_control_negativo$antibiotic), mean,  na.rm = TRUE)
  
  return(mean_negative_control)
}


# Apply the function to each dataframe
mean_negative_control_1 <- calculate_mean_negative_control(compounds1)
mean_negative_control_2 <- calculate_mean_negative_control(compounds2)

#Substract baseline 
compounds1$Absorbance_600 <-compounds1$Absorbance_600 - mean_negative_control_1
compounds2$Absorbance_600 <-compounds2$Absorbance_600 - mean_negative_control_2

#=============================================================================================================
antibiotics = c("Positive_control" , "RIF", "LEVO", "FUS", "CLI", "LZD",  "SXT", "OXA", "VAN",  "FOX")

create_Tabla_Ensayo_control <- function(raw_data, antibiotics) {
  # Create an empty list to store the subsets
  Tabla_Ensayo_control <- vector("list", length(antibiotics))

  # Loop through the antibiotics
  for (i in 1:length(antibiotics)) {
    Tabla_Ensayo_control[[i]] <- subset(raw_data, antibiotic == antibiotics[i], select = c(Absorbance_600, Concentration, antibiotic))
  }

  return(Tabla_Ensayo_control)
}

# Example usage with compounds1
Tabla_Ensayo_control1 <- create_Tabla_Ensayo_control(compounds1, antibiotics)

# Example usage with compounds2
Tabla_Ensayo_control2 <- create_Tabla_Ensayo_control(compounds2, antibiotics)

#=============================================================================================================
# Define a condition (e.g., remove rows where Value is less than 0.1) for the first dataframe POSITIVE CONTROL
condition <- 0.1

# Access the second dataframe and remove rows based on the condition
Tabla_Ensayo_control1[[1]] <- Tabla_Ensayo_control1[[1]][Tabla_Ensayo_control1[[1]]$Absorbance_600 >= condition, ]
Tabla_Ensayo_control2[[1]] <- Tabla_Ensayo_control2[[1]][Tabla_Ensayo_control2[[1]]$Absorbance_600 >= condition, ]

#=============================================================================================================

create_mean_list <- function(Tabla_Ensayo_control) {
  mean_list <- list()
  
  # Loop through the elements in Tabla_Ensayo_control
  for (i in 1:length(Tabla_Ensayo_control)) {
    mean_list[[i]] <- tapply(
      Tabla_Ensayo_control[[i]]$Absorbance_600,
      factor(Tabla_Ensayo_control[[i]]$Concentration),
      mean
    )
  }
  
  return(mean_list)
}

# Example usage with Tabla_Ensayo_control1
mean_list1 <- create_mean_list(Tabla_Ensayo_control1)

# Example usage with Tabla_Ensayo_control2
mean_list2 <- create_mean_list(Tabla_Ensayo_control2)


L11=mean_list1[[1]]
L12=mean_list2[[1]]
L21=mean_list1[[2]]
L22=mean_list2[[2]]
L31=mean_list1[[3]]
L32=mean_list2[[3]]
L41=mean_list1[[4]]
L42=mean_list2[[4]]
L51=mean_list1[[5]]
L52=mean_list2[[5]]
L61=mean_list1[[6]]
L62=mean_list2[[6]]
L71=mean_list1[[7]]
L72=mean_list2[[7]]
L81=mean_list1[[8]]
L82=mean_list2[[8]]
L91=mean_list1[[9]]
L92=mean_list2[[9]]
L101=mean_list1[[10]]
L102=mean_list2[[10]]


#=============================================================================================================
#Make dataframes
# Create some example lists
concentration_0 = c(0) 

df1 <- as.data.frame(t(rbind(L11, L12)))
df2 <- as.data.frame(t(rbind(L21, L22)))
df3 <- as.data.frame(t(rbind(L31, L32)))
df4 <- as.data.frame(t(rbind(L41, L42)))
df5 <- as.data.frame(t(rbind(L51, L52)))
df6 <- as.data.frame(t(rbind(L61, L62)))
df7 <- as.data.frame(t(rbind(L71, L72)))
df8 <- as.data.frame(t(rbind(L81, L82)))
df9 <- as.data.frame(t(rbind(L91, L92)))
df10 <- as.data.frame(t(rbind(L101, L102)))


# Create a list of your dataframes
df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)

# List to store the results
result_list <- list()

# Loop through the dataframes
for (df in df_list) {
  # Use the stack function to reshape the dataframe
  reshaped_data <- stack(df)
  
  # Rename the columns
  colnames(reshaped_data) <- c("OD", "Name")
  
  # Append the reshaped dataframe to the result list
  result_list[[length(result_list) + 1]] <- reshaped_data
}

# Access the individual results  

df1_f = result_list[[1]]
df2_f = result_list[[2]]
df3_f = result_list[[3]]
df4_f = result_list[[4]]
df5_f = result_list[[5]]
df6_f = result_list[[6]]
df7_f = result_list[[7]]
df8_f = result_list[[8]]
df9_f = result_list[[9]]
df10_f = result_list[[10]]

# Add a concentration column

df1_f$Conc = rev(c(concentration_0,concentration_0))
df2_f$Conc = rev(c(concentration1,concentration1))
df3_f$Conc = rev(c(concentration1,concentration1))
df4_f$Conc = rev(c(concentration2,concentration2))
df5_f$Conc = rev(c(concentration3,concentration3))
df6_f$Conc = rev(c(concentration4,concentration4))
df7_f$Conc = rev(c(concentration5,concentration5))
df8_f$Conc = rev(c(concentration6,concentration6))
df9_f$Conc = rev(c(concentration7,concentration7))
df10_f$Conc = rev(c(concentration8,concentration8))


df2_f = rbind(df1_f, df2_f)
df3_f = rbind(df1_f, df3_f)
df4_f = rbind(df1_f, df4_f)
df5_f = rbind(df1_f, df5_f)
df6_f = rbind(df1_f, df6_f)
df7_f = rbind(df1_f, df7_f)
df8_f = rbind(df1_f, df8_f)
df9_f = rbind(df1_f, df9_f)
df10_f = rbind(df1_f, df10_f)


#================================================================================================================
#DunnettTest
# List of dataframes
df_list <- list(df2_f, df3_f, df4_f, df5_f, df6_f, df7_f, df8_f, df9_f, df10_f)

# Lists to store the results
aov_results <- list()
Dunnett_results <- list()

# Loop through the dataframes
for (i in 1:9) {
  # Perform ANOVA
  aov_result <- aov(OD ~ Conc, data = df_list[[i]])
  summary_aov <- summary(aov_result)
  aov_results[[i]] <- summary_aov

  # Perform Dunnett's Test
  Dunnett_result <- DunnettTest(x = df_list[[i]]$OD, g = df_list[[i]]$Conc)
  Dunnett_results[[i]] <- Dunnett_result
}

# Access individual results 

Dunnett_df2_f = Dunnett_results[[1]]
Dunnett_df3_f =  Dunnett_results[[2]]
Dunnett_df4_f =  Dunnett_results[[3]]
Dunnett_df5_f =  Dunnett_results[[4]]
Dunnett_df6_f =  Dunnett_results[[5]]
Dunnett_df7_f =  Dunnett_results[[6]]
Dunnett_df8_f =  Dunnett_results[[7]]
Dunnett_df9_f =  Dunnett_results[[8]]
Dunnett_df10_f =  Dunnett_results[[9]]


#=============================================================================================================

# List of dataframes
df_list <- list(df2_f, df3_f, df4_f, df5_f, df6_f, df7_f, df8_f, df9_f, df10_f)

# Lists to store the results
mean_results <- list()
sd_results <- list()

# Loop through the dataframes
for (i in 1:9) {
  # Means
  mean_result <- tapply(df_list[[i]]$OD, factor(df_list[[i]]$Conc), mean)
  mean_results[[i]] <- mean_result

  # SDs
  sd_result <- tapply(df_list[[i]]$OD, factor(df_list[[i]]$Conc), sd)
  sd_results[[i]] <- sd_result
}

# Access individual results 
 mean_df2_f =mean_results[[1]] 
 sd_df2_f = sd_results[[1]] 
 mean_df3_f =mean_results[[2]] 
 sd_df3_f = sd_results[[2]] 
 mean_df4_f =mean_results[[3]] 
 sd_df4_f = sd_results[[3]] 
 mean_df5_f =mean_results[[4]] 
 sd_df5_f = sd_results[[4]] 
 mean_df6_f =mean_results[[5]] 
 sd_df6_f = sd_results[[5]] 
 mean_df7_f =mean_results[[6]] 
 sd_df7_f = sd_results[[6]] 
 mean_df8_f =mean_results[[7]] 
 sd_df8_f = sd_results[[7]] 
 mean_df9_f =mean_results[[8]] 
 sd_df9_f = sd_results[[8]] 
 mean_df10_f =mean_results[[9]] 
 sd_df10_f = sd_results[[9]] 
 
#=============================================================================================================
#Plots

generate_plot_df <- function(concentration, mean, sd) {
  df_plot = data.frame(
    "Conc" = rev(c(concentration, concentration_0)),
    "mean" = mean,
    "SD" = sd
  )
  return(df_plot)
}

df_plot2 = generate_plot_df(concentration1, mean_df2_f, sd_df2_f)
df_plot3 = generate_plot_df(concentration1, mean_df3_f, sd_df3_f)
df_plot4 = generate_plot_df(concentration2, mean_df4_f, sd_df4_f)
df_plot5 = generate_plot_df(concentration3, mean_df5_f, sd_df5_f)
df_plot6 = generate_plot_df(concentration4, mean_df6_f, sd_df6_f)
df_plot7 = generate_plot_df(concentration5, mean_df7_f, sd_df7_f)
df_plot8 = generate_plot_df(concentration6, mean_df8_f, sd_df8_f)
df_plot9 = generate_plot_df(concentration7, mean_df9_f, sd_df9_f)
df_plot10 = generate_plot_df(concentration8, mean_df10_f, sd_df10_f)

#======================================================================================================================


titles = c("Rifampin (RIF)", "Levofloxacin (LEVO)", "Fusidate (FUS)", "Clindamycin (CLI)", "Linezolid (LZD)", "Trimethoprim/Sulfamethoxazole (SXT)",
 "Oxacilin + 2% NaCl (OXA+)", "Vancomycin (VAN)", "Cefoxitin (FOX)")

# These are obtained from the postestScreening (TO DO: AUTOMATE THIS)
significances = list()
significances[[1]] = c("", "", "", "", "", "", "", "", "", "", ".", "", ".", ".", "")
significances[[2]] = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
significances[[3]] = c("", "", "", "", "", "", "", "", "", "", "*")
significances[[4]] = c("", "", "", "", "", "", "", "", "", "", "", "", "")
significances[[5]] = c("", "", "", "", "", "", "", "") 
significances[[6]] = c("", "", "", "", "", "", "", "") 
significances[[7]] = c("", ".", "", "", "", "", "", "", "") 
significances[[8]] = c("", "", "", "", "", "", "", "", "","") 
significances[[9]] = c("", "", "", "", "", "", "", "")  

# Define a distance for plotting and store in a list of the correct size
dist <- 0.08
for (i in 1:9) {
  distance[[i]] <- rep(dist, length(significances[[i]]))
}

dfs_plot = list(df_plot2, df_plot3, df_plot4, df_plot5, df_plot6, df_plot7, df_plot8, df_plot9, df_plot10)

#obliga a ordenar a que el eje x se ordene de manera equidistante
for (i in 1:9){
dfs_plot[[i]]$Conc <- factor(dfs_plot[[i]]$Conc, levels = dfs_plot[[i]]$Conc)
}

for (i in 1:9){
  Mean_i = dfs_plot[[i]]$mean
  SD_i = dfs_plot[[i]]$SD
  label = significances[[i]]
  distances = distance[[i]]

  # First check if lengths of vectors are consistent
  if (length(Mean_i) == length(SD_i) && length(Mean_i) == length(distances) && length(Mean_i) == length(label)) {
  p <- ggplot(data=dfs_plot[[i]], aes(x=Conc, y=mean)) + 
          theme_bw() + 
          geom_bar(stat="identity", fill="deeppink4")+ 
          geom_errorbar(aes(ymin=Mean_i-SD_i, ymax=Mean_i+SD_i), width=.2, position=position_dodge(0.05)) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.8)) +
          theme(legend.position="top") + 
          labs(title=element_blank(), x=expression(paste("Compound concentration (", mu, "M)")), y = "OD_600")+
          ggtitle(titles[i])

  # Add annotations using a loop
  for (j in 1:length(Mean_i)){
    p <- p + annotate("text", x = j+0.01, y = Mean_i[j] + SD_i[j] + distances[j], label = label[j],
                      size = 9, family = "Times", fontface = "bold.italic", colour = "blue")
  }

  p_list[[i]] <- p} else {
    warning(paste("Skipping plot_ATB", i, "due to inconsistent vector lengths"))
  }

  name_plot = paste( paste("plot_ATB", i, sep="_"), ".jpg", sep="")
  ggsave(name_plot)
}


#---- Fuera de R

montage -geometry +2+2 -tile 3x3 plot_ATB_{1..9}.jpg output.jpg
