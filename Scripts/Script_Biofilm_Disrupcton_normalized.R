
 #========PART 1=================================================================
 #PREPARE THE DATA SET TO BE USEFUL FOR R
 filename = "Biofilm_Pseudomona_compuesto_34_CIMb_Cristal_violeta.csv"
 raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")


column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
colnames(raw_data) = column_names

data <- subset( raw_data, select = -c(X))

n_m = 12 
n_b = 96    # amount of wells for each bacteria
n_row = 6
n_col = 12
#n_c = n_row * (n_col/n_m)

# Create a new dataframe separating the measurements for each compound
# The first column has the compound name, and the second column is the measurement
compounds = data.frame()
sum = 0
for(i in 1:n_row) {
	for(j in 1:n_col) {
		sum = sum + 1

		# Obtain the index "i_c" for each compound
		i_c = as.integer( (sum-1)/n_m ) + 1
		j_c = (j-1)%%12 + 1
		name_c = paste('Condition', sprintf("%02d", i_c), sep = "_")
		compounds[n_m*(i_c-1)+j_c,1] = name_c
		compounds[n_m*(i_c-1)+j_c,2] = data[i,j]

		# Obtain the index "i_b" for each bacteria
		i_b = as.integer( (sum-1)/n_b ) + 1
		name_b = paste('bacteria', as.character(i_b), sep = "_")
		compounds[n_m*(i_c-1)+j_c,3] = name_b
	}
}


# Add a header
column_names = c('Condition', 'Absorbance', 'Bacteria')
colnames(compounds) = column_names


bacteria_condition =  c("Replicate_1","Replicate_2","Replicate_3","Replicate_4","Positive_control","Negative_control")


#if the plate is inverted
#bacteria_condition =  c("Positive_control","Replicate_4","Replicate_3","Negative_control",
#	"Replicate_2","Negative_control","Replicate_1","Negative_control")


	for(i in 1:6) {
	# Retrieve a string as "Condition_i" for each "i"
	name_c = paste('Condition', sprintf("%02d", i), sep = "_")

	# Replace when the compounds has this name
	compounds[compounds == name_c] <- bacteria_condition[i]
    }


# Use this to define names and concentrations of the compounds.
# The first value is the MAXIMUM concentration used.
concentration = c(7.2) 
for(i in 2:12) {concentration[i] = concentration[i-1]/2}

#add column called 'concentration'

concentration = round(concentration,  digits = 3)

compounds$concentration <- concentration

#if the plate is inverted
#compounds$concentration <- rev(concentration)


compounds["Bacteria"][compounds["Bacteria"] == "bacteria_1"]  <- "Staphylococcus aureus"

#==NORMALIZE THE DATA BY THE POSITIVE AND NEGATIVE CONTROL============

Tabla_Ensayo_control_negativo = subset(compounds, Condition =="Negative_control", select=c(Bacteria, concentration, Absorbance))

mean_negative_control=tapply(Tabla_Ensayo_control_negativo$Absorbance,factor(Tabla_Ensayo_control_negativo$Bacteria),mean)

Tabla_Ensayo_control_positivo = subset(compounds, Condition =="Positive_control", select=c(Bacteria, concentration, Absorbance))

mean_positive_control=tapply(Tabla_Ensayo_control_positivo$Absorbance,factor(Tabla_Ensayo_control_positivo$Bacteria),mean)


Tabla_Ensayo<-subset(compounds, (Condition =="Replicate_1" | Condition =="Replicate_2" | Condition =="Replicate_3"| Condition =="Positive_control"), select=c(Bacteria, concentration, Absorbance, Condition))

#replace the values of the concentration of the positive control 

Tabla_Ensayo$concentration[Tabla_Ensayo$Condition == "Positive_control"] <- 0

Tabla_Ensayo$Absorbance<-Tabla_Ensayo$Absorbance/mean_positive_control

#Tabla_Ensayo$Absorbance <-Tabla_Ensayo$Absorbance-mean_negative_control




#Tabla_Ensayo$AbsorbanceTabla_Ensayo$Absorbance/mean_positive_control

#ANALISIS DE LOS PARAMETROS ESTADISTICOS=============================================

library("agricolae")
library("ggplot2")

aovScreening<-aov(Absorbance~concentration,data=Tabla_Ensayo)
summary(aovScreening)

postestScreening<-LSD.test(aovScreening, "concentration")
print(postestScreening)

#PLOTS=====================================================================================

mean_absorbance <- tapply(Tabla_Ensayo$Absorbance,factor(Tabla_Ensayo$concentration),mean)

sd_absorbance <- tapply(Tabla_Ensayo$Absorbance,factor(Tabla_Ensayo$concentration),sd)

#add the value 0 to the vector concentration
concentration <- append(concentration, 0)

df_plot = data.frame(
		"Conc" = rev(concentration),
		"Mean" = mean_absorbance, 
		"SD" = sd_absorbance
		)

#obliga a ordenar a que el eje x se ordene de manera equidistante
df_plot$Conc <- factor(df_plot$Conc, levels = df_plot$Conc)


#================================================================================

# Make the barplot

#title = "Pseudomonas putida JM MIS"
#title = "Escherichia coli JM MIS"
#title = "Escherichia coli JM MIR"
#title = "Staphylococcus epidermis JM MIS"
#title = "Staphylococcus epidermis JM MIR"
#title = "Pseudomonas putida JM MIR"
title = "Pseudomonas putida JM MIS"



#significances = c("a", "b", "bc", "ab", "ab", "ab", "b", "bc", "b", "ab", "bc", "cd", "d")
#significances = c("ab", "a", "a", "a", "ab", "abc", "bcd", "cd", "cd", "d", "d", "d", "d")
#significances = c("ab", "abc", "bcd", "ab", "abc", "bcd", "d", "cd", "cd", "cd", "cd", "a", "abc")
#significances = c("", "", "", "", "", "", "", "", "", "", "", "", "")
#significances = c("a", "a", "a", "ab", "ab", "abc", "c", "c", "c", "c", "c", "bc", "abc")
#significances = c("a", "abc", "abc", "ab", "ab", "bcd", "d", "d", "d", "d", "cd", "abcd", "a")
#significances = c("a", "b", "bc", "ab", "ab", "ab", "b", "bc", "b", "ab", "bc", "cd", "d")
significances = c("a", "abcd", "bcde", "abcde", "abcd", "abc", "abcde", "ab", "ab", "abc", "cde", "e", "de")



distance = c(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)
xc = c(1,2,3,4,5,6,7,8,9,10,11,12,13)
textSize=8

for (i in 1:12){
	Mean_i = df_plot$Mean
	SD_i = df_plot$SD
    distances = distance
	p <- ggplot(data=df_plot, aes(x=Conc, y=Mean)) + 
	            theme_bw() +
	  			geom_bar(stat="identity",  width=0.9, fill = "darkmagenta")+ 
	  			ylim(0, 1.5)+
	  			geom_errorbar(aes(ymin=Mean_i-SD_i, ymax=Mean_i+SD_i), width=0.2,
	      			position=position_dodge(0.05)) +  
	  			    theme(axis.text.x = element_text(color="black", size=12, angle=10)) +
	      			theme(legend.position="top",legend.title = element_text(size = textSize), legend.text  = element_text(size = textSize)) + 
	                     labs(title=element_blank(), x=expression(paste("Compound concentration (", mu, "g/",mu, "l)")), y = "Normalized biofilm growth")+
	                      annotate("text", x = xc, y = Mean_i+SD_i+distances, label = significances, size = 7, family = "Times", fontface = "bold.italic", colour = "black")+
	            	                                                   ggtitle(title)
 
	#p + theme_minimal()
}

	    		    		

#ggsave("Pseudomonas_JMMIR_normalizado1.pdf")
ggsave("Psudo_JMMIS_normalizado.pdf")