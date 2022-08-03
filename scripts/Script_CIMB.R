
 #========PART 1=================================================================
 #PREPARE THE DATA SET TO BE USEFUL FOR R
 filename = "Biofilm_E.coli_compuesto_35_CIMb.csv"
 raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")


column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
colnames(raw_data) = column_names

data <- subset( raw_data, select = -c(X))

n_m = 12 
n_b = 96    # amount of wells for each bacteria
n_row = 8
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
column_names = c('Condition', 'Fluorescence', 'Bacteria')
colnames(compounds) = column_names


bacteria_condition =  c("Negative_control","Replicate_1","Negative_control","Replicate_2",
	"Negative_control","Replicate_3","Replicate_4","Positive_control")


#if the plate is inverted
#bacteria_condition =  c("Positive_control","Replicate_4","Replicate_3","Negative_control",
#	"Replicate_2","Negative_control","Replicate_1","Negative_control")

	for(i in 1:8) {
	# Retrieve a string as "Condition_i" for each "i"
	name_c = paste('Condition', sprintf("%02d", i), sep = "_")

	# Replace when the compounds has this name
	compounds[compounds == name_c] <- bacteria_condition[i]
}


# Use this to define names and concentrations of the compounds.
# The first value is the MAXIMUM concentration used.
concentration = c(9.6) 
for(i in 2:12) {concentration[i] = concentration[i-1]/2}

#add column called 'concentration'
compounds$concentration <- concentration

#if the plate is inverted
#compounds$concentration <- rev(concentration)


compounds["Bacteria"][compounds["Bacteria"] == "bacteria_1"]  <- "Staphylococcus aureus"

#==NORMALIZE THE DATA BY THE POSITIVE AND NEGATIVE CONTROL============

Tabla_Ensayo_control_negativo = subset(compounds, Condition =="Negative_control", select=c(Bacteria, concentration, Fluorescence))

mean_negative_control=tapply(Tabla_Ensayo_control_negativo$Fluorescence,factor(Tabla_Ensayo_control_negativo$Bacteria),mean)

Tabla_Ensayo_control_positivo = subset(compounds, Condition =="Positive_control", select=c(Bacteria, concentration, Fluorescence))

mean_positive_control=tapply(Tabla_Ensayo_control_positivo$Fluorescence,factor(Tabla_Ensayo_control_positivo$Bacteria),mean)

#ANALISIS DE LOS PARAMETROS ESTADISTICOS=============================================

library("drc")

Tabla_Ensayo<-subset(compounds, (Condition =="Replicate_1" | Condition =="Replicate_2" | Condition =="Replicate_3"), select=c(Bacteria, concentration, Fluorescence))

Tabla_Ensayo$Fluorescence <-Tabla_Ensayo$Fluorescence-mean_negative_control

Tabla_Ensayo$Fluorescence <-Tabla_Ensayo$Fluorescence/mean_positive_control


Modelo_Ensayo<- drm(Fluorescence ~ concentration,
                   data = Tabla_Ensayo, 
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, 0, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("CIMb")
print(summary(Modelo_Ensayo))
print(confint(Modelo_Ensayo, level = 0.90))
EC90Comp<-ED(Modelo_Ensayo,c(10,50,90),interval="delta")


pdf("Titration_curve_generic_CIMb.pdf")

plot(Modelo_Ensayo,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText= "Escherichia coli - JM MIR", xlab = expression(paste("Compound concentration (", mu, "g/ ", mu, "l)")), ylab = "Viability", main = "CIMb")


plot(Modelo_Ensayo, col= "blue", add=TRUE, type='confidence')

dev.off()