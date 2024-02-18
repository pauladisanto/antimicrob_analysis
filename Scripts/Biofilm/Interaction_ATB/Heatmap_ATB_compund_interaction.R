#heatmap
library(reshape)                                                
library(ggplot2)

read_table <- function(filename) {  
	# Read the raw data from the file, skipping the first 14 rows which contain irrelevant information
	raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")
}

# define the files to analyse

filename1 = "Interaction_P_putida_C1_30_plus_ATB_AFTER_recovery_media_replica_1.csv"
filename2 = "Interaction_P_putida_C1_30_plus_ATB_AFTER_recovery_media_replica_2.csv"

# read the tables
raw_data1_table = read_table(filename1)
raw_data2_table = read_table(filename2)

#==================================================================================

# Dummy column that is later removed
column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')

# Define the column names
colnames(raw_data1_table) = column_names
colnames(raw_data2_table) = column_names
raw_data1 <- subset(raw_data1_table, select = -c(X))
raw_data2 <- subset(raw_data2_table, select = -c(X))

average_raw_data <- (raw_data1 + raw_data2) / 2

#Normalization by positive control
positive_control = average_raw_data[3:5,1] 
positive_control = unlist(positive_control,recursive = F, use.names = F)    
mean_positive_control = mean(positive_control)

average_raw_data = average_raw_data / mean_positive_control

#Add column and row names
colnames(average_raw_data)<- paste0(1:12)                             
rownames(average_raw_data)<- paste0(1:8) 


average_raw_data = as.matrix(average_raw_data)

data_melt1 <- melt(average_raw_data)                                           
ggp1 <- ggplot(data_melt1, aes(X2, X1)) +                           
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "darkred")+
  scale_y_reverse()  # Flip the plot horizontally 

print(ggp1)
#ggsave("C1_ATB_biofilm_30_S_Pseudomonas.png")
ggsave("C1_ATB_biofilm_30_S_Pseudomonas.png", ggp1, width = 8, height = 4, units = "in")









