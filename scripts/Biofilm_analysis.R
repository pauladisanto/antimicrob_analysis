
read_table <- function(filename) {  
	# Read the raw data from the file, skipping the first 14 rows which contain irrelevant information
	raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")
}

# define the files to analyse
filename1 = "Absorbancia550nm_ensayoBF_MH_24h_02_y_2_glucosa.csv"
filename2 = "Absorbancia550nm_ensayoBF_MH_48h_02_y_2_glucosa.csv"

# read the tables
raw_data1 = read_table(filename1)
raw_data2 = read_table(filename2)

# merge the two data frames into one
compounds_merged <- rbind(raw_data1, raw_data2)


#===================================================================================================

# define the column names
column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
colnames(compounds_merged) = column_names

# remove the first column
data <- subset( compounds_merged, select = -c(X))

# We have n_m measurements for each of the n_c conditions
n_m = 3 
n_b = 12    # amount of wells for each bacteria
n_row = 16
n_col = 12
n_c = n_row * (n_col/n_m) # number of conditions
n_bac = 7   # number of bacteria (including the negative control)

# Create a new dataframe separating the measurements for each compound
# The first column has the compound name, and the second column is the measurement
compounds = data.frame()
sum = 0
for(i in 1:n_row) {
	for(j in 1:n_col) {
		sum = sum + 1

		# Obtain the index "i_c" for each compound
		i_c = as.integer( (sum-1)/n_m ) + 1
		j_c = (j-1)%%3 + 1
		name_c = paste('Condition', sprintf("%02d", i_c), sep = "_")
		compounds[n_m*(i_c-1)+j_c,1] = name_c
		compounds[n_m*(i_c-1)+j_c,2] = data[i,j]

		# Obtain the index "i_b" for each bacteria
		i_b = as.integer( (sum-1)/n_b ) + 1
		name_b = paste('bacteria', as.character(i_b), sep = "_")
		compounds[n_m*(i_c-1)+j_c,3] = name_b
	}
}


#===============================================================================================================

#bacteria 24 h incubation
#bacteria_condition =  c("EntC 24h 30C 0.2GLU", "EntC 24h 37C 0.2GLU", "EntC 24h 30C 2GLU", "EntC 24h 37C 2GLU",
#											"StaphE 24h 30C 0.2GLU", "StaphE 24h 37C 0.2GLU", "StaphE 24h 30C 2GLU", "StaphE 24h 37C 2GLU",
#											"EColi 24h 30C 0.2GLU", "EColi 24h 37C 0.2GLU", "EColi 24h 30C 2GLU", "EColi 24h 37C 2GLU",
#											"Bsub 24h 30C 0.2GLU", "Bsub 24h 37C 0.2GLU", "Bsub 24h 30C 2GLU", "Bsub 24h 37C 2GLU",
#											"Pput 24h 30C 0.2GLU", "Pput 24h 37C 0.2GLU", "Pput 24h 30C 2GLU", "Pput 24h 37C 2GLU",
#										  "PC 24h 30C 0.2GLU", "PC 24h 37C 0.2GLU", "PC 24h 30C 2GLU", "PC 24h 37C 2GLU",
#										  "PB 24h 30C 0.2GLU", "PB 24h 37C 0.2GLU", "PB 24h 30C 2GLU", "PB 24h 37C 2GLU",
#										  "NEG", "NEG", "NEG", "NEG")


#bacteria 48 h incubation
#bacteria_condition =  c("EntC 48h 30C 0.2GLU", "EntC 48h 37C 0.2GLU", "EntC 48h 30C 2GLU", "EntC 48h 37C 2GLU",
#						"StaphE 48h 30C 0.2GLU", "StaphE 48h 37C 0.2GLU", "StaphE 48h 30C 2GLU", "StaphE 48h 37C 2GLU",
#						"EColi 48h 30C 0.2GLU", "EColi 48h 37C 0.2GLU", "EColi 48h 30C 2GLU", "EColi 48h 37C 2GLU",
#						"Bsub 48h 30C 0.2GLU", "Bsub 48h 37C 0.2GLU", "Bsub 48h 30C 2GLU", "Bsub 48h 37C 2GLU",
#						"Pput 48h 30C 0.2GLU", "Pput 48h 37C 0.2GLU", "Pput 48h 30C 2GLU", "Pput 48h 37C 2GLU",
#						"PC 48h 30C 0.2GLU", "PC 48h 37C 0.2GLU", "PC 48h 30C 2GLU", "PC 48h 37C 2GLU",
#						"PB 48h 30C 0.2GLU", "PB 48h 37C 0.2GLU", "PB 48h 30C 2GLU", "PB 48h 37C 2GLU",
#						"NEG", "NEG", "NEG", "NEG")


bacteria_condition =  c("EntC 24h 30C 0.2GLU", "EntC 24h 37C 0.2GLU", "EntC 24h 30C 2GLU", "EntC 24h 37C 2GLU",
						"StaphE 24h 30C 0.2GLU", "StaphE 24h 37C 0.2GLU", "StaphE 24h 30C 2GLU", "StaphE 24h 37C 2GLU",
						"EColi 24h 30C 0.2GLU", "EColi 24h 37C 0.2GLU", "EColi 24h 30C 2GLU", "EColi 24h 37C 2GLU",
						"Bsub 24h 30C 0.2GLU", "Bsub 24h 37C 0.2GLU", "Bsub 24h 30C 2GLU", "Bsub 24h 37C 2GLU",
						"Pput 24h 30C 0.2GLU", "Pput 24h 37C 0.2GLU", "Pput 24h 30C 2GLU", "Pput 24h 37C 2GLU",
					    "PC 24h 30C 0.2GLU", "PC 24h 37C 0.2GLU", "PC 24h 30C 2GLU", "PC 24h 37C 2GLU",
						"PB 24h 30C 0.2GLU", "PB 24h 37C 0.2GLU", "PB 24h 30C 2GLU", "PB 24h 37C 2GLU",
						"NEG", "NEG", "NEG", "NEG","EntC 48h 30C 0.2GLU", "EntC 48h 37C 0.2GLU", "EntC 48h 30C 2GLU", "EntC 48h 37C 2GLU",
						"StaphE 48h 30C 0.2GLU", "StaphE 48h 37C 0.2GLU", "StaphE 48h 30C 2GLU", "StaphE 48h 37C 2GLU",
						"EColi 48h 30C 0.2GLU", "EColi 48h 37C 0.2GLU", "EColi 48h 30C 2GLU", "EColi 48h 37C 2GLU",
						"Bsub 48h 30C 0.2GLU", "Bsub 48h 37C 0.2GLU", "Bsub 48h 30C 2GLU", "Bsub 48h 37C 2GLU",
						"Pput 48h 30C 0.2GLU", "Pput 48h 37C 0.2GLU", "Pput 48h 30C 2GLU", "Pput 48h 37C 2GLU",
						"PC 48h 30C 0.2GLU", "PC 48h 37C 0.2GLU", "PC 48h 30C 2GLU", "PC 48h 37C 2GLU",
						"PB 48h 30C 0.2GLU", "PB 48h 37C 0.2GLU", "PB 48h 30C 2GLU", "PB 48h 37C 2GLU",
						"NEG", "NEG", "NEG", "NEG")


for(i in 1:n_c) {
	# Retrieve a string as "Condition_i" for each "i"
	name_c = paste('Condition', sprintf("%02d", i), sep = "_")

	# Replace when the compounds has this name
	compounds[compounds == name_c] <- bacteria_condition[i]
}


# Add a header
column_names = c('Condition', 'Absorbance', 'Bacteria')
colnames(compounds) = column_names

# select blanks and calculate the media to rest the value to the crystal violect data
negative_control <- subset(compounds, Bacteria == "bacteria_8" | Bacteria == "bacteria_16" , select=c(Condition, Absorbance, Bacteria))

mean_negative_control <- tapply(negative_control$Absorbance,factor(negative_control$Condition),mean)

compounds<-transform(compounds,Absorbance=Absorbance-mean_negative_control)


#=========================================subselect tables====================
#df["Column Name"][df["Column Name"] == "Old Value"] <- "New Value"

#DEFINIR LAS BACTERIAS DE 9 A 16 COMO SUS PARES DE 1 A 8
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_9"]  <- "bacteria_1"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_10"]  <- "bacteria_2"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_11"]  <- "bacteria_3"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_12"]  <- "bacteria_4"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_13"]  <- "bacteria_5"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_14"]  <- "bacteria_6"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_15"]  <- "bacteria_7"
compounds["Bacteria"][compounds["Bacteria"] == "bacteria_16"]  <- "bacteria_8"


# Make a list of all the dataframes. The element "i" of the list is the df for the "i" bacteria.
dfs = list()
for(i in 1:n_row) {
	dfs[[i]] <- subset(compounds, Bacteria == paste("bacteria", i, sep = "_"), select=c(Condition, Absorbance, Bacteria))
}


#============================================post-test==============================================

library(agricolae)

aovScreening1<-aov(Absorbance~Condition,data=dfs[[1]])
summary(aovScreening1)
postestScreening1<-LSD.test(aovScreening1,"Condition")
print(postestScreening1)

aovScreening2<-aov(Absorbance~Condition,data=dfs[[2]])
summary(aovScreening2)
postestScreening2<-LSD.test(aovScreening2,"Condition")
print(postestScreening2)

aovScreening3<-aov(Absorbance~Condition,data=dfs[[3]])
summary(aovScreening3)
postestScreening3<-LSD.test(aovScreening3,"Condition")
print(postestScreening3)

aovScreening4<-aov(Absorbance~Condition,data=dfs[[4]])
summary(aovScreening4)
postestScreening4<-LSD.test(aovScreening4,"Condition")
print(postestScreening4)

aovScreening5<-aov(Absorbance~Condition,data=dfs[[5]])
summary(aovScreening5)
postestScreening5<-LSD.test(aovScreening5,"Condition")
print(postestScreening5)

aovScreening6<-aov(Absorbance~Condition,data=dfs[[6]])
summary(aovScreening6)
postestScreening6<-LSD.test(aovScreening6,"Condition")
print(postestScreening6)

aovScreening7<-aov(Absorbance~Condition,data=dfs[[7]])
summary(aovScreening7)
postestScreening7<-LSD.test(aovScreening7,"Condition")
print(postestScreening7)

#=========================================================================

# CALCULATIONS TO MAKE THE PLOT

#=========================================================================

mean_conditions = list()
sd_conditions = list()
conditions = list()
dfs_plot = list() 
for(i in 1:n_bac){
	mean_conditions[[i]] <- tapply(dfs[[i]]$Absorbance,factor(dfs[[i]]$Condition),mean)
	sd_conditions[[i]] <- tapply(dfs[[i]]$Absorbance,factor(dfs[[i]]$Condition),sd)
	conditions[[i]] = names(mean_conditions[[i]])
	dfs_plot[[i]] = data.frame(
		"Cond" = conditions[[i]],
		"Mean" = mean_conditions[[i]], 
		"SD" = sd_conditions[[i]]
		)
}


#=========================================================================

library("ggplot2")

colors = c(1,2,3,4,5,6,7) #TO DO: DEFINE NICE COLOURS FOR FILLING!
titles = c("Enterococcus raffinosus", "Staphylococcus epidermidis", "Escherichia coli", 
		"Bacillus subtilis",  "Pseudomonas Putida",  "Pectobacterium carotovorum", 
		"Paraburkholderia caledonica")

# These are obtained from the postestScreening (TO DO: AUTOMATE THIS)
significances = list()
significances[[1]] = c("c", "bc", "bc", "a", "c", "bc", "c", "b")
significances[[2]] = c("d", "cd","e", "c", "a", "a", "b", "a")
significances[[3]] = c("d", "d", "b", "d", "c", "c", "a", "c")
significances[[4]] = c("g", "c", "fg", "de", "a", "b", "cd", "ef")
significances[[5]] = c("bc", "abc", "c", "a", "abc", "bc", "ab", "ab")
significances[[6]] = c("c", "ab", "bc", "c", "a", "bc", "bc", "bc")
significances[[7]] = c("b", "b", "b", "a", "b", "b", "b", "b")

distance = list()
distance[[1]] = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
distance[[2]] = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
distance[[3]] = c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
distance[[4]] = c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
distance[[5]] = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
distance[[6]] = c(0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005)
distance[[7]] = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
distance[[8]] = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)

label = c()
for (i in 1:n_bac){
	Mean_i = dfs_plot[[i]]$Mean
	SD_i = dfs_plot[[i]]$SD
  	label = significances[[i]]
  	distances = distance[[i]]
	p <- ggplot(data=dfs_plot[[i]], aes(x=Cond, y=Mean)) + 
	  			geom_bar(stat="identity", fill=i)+ 
	  			geom_errorbar(aes(ymin=Mean_i-SD_i, ymax=Mean_i+SD_i), width=.2,
	      		position=position_dodge(0.05)) + 
	      		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.8)) +
	      		theme(legend.position="top") + 
	            labs(title=element_blank(), x="Condition", y = "Absorbance")+
	            annotate("text", x = 1, y = Mean_i[1]+SD_i[1]+distances[1], label = label[1],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 2, y = Mean_i[2]+SD_i[2]+distances[2], label = label[2],
	    		   		size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 3, y = Mean_i[3]+SD_i[3]+distances[3], label = label[3],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 4, y = Mean_i[4]+SD_i[4]+distances[4], label = label[4],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 5, y = Mean_i[5]+SD_i[5]+distances[5], label = label[5],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	   			annotate("text", x = 6, y = Mean_i[6]+SD_i[6]+distances[6], label = label[6],
	   					size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 7, y = Mean_i[7]+SD_i[7]+distances[7], label = label[7],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	    		annotate("text", x = 8, y = Mean_i[8]+SD_i[8]+distances[8], label = label[8],
	    				size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	            ggtitle(titles[i])

	print(p)
	name_plot = paste( paste("plot_biofilm", i, sep="_"), ".jpg", sep="")
	ggsave(name_plot)
}




