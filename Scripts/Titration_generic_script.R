
#script titration

# We read the raw data from the file, skipping the first 14 rows which contain irrelevant information

#source("Titration_generic_script.R")
filename = "tritration_curve_first_assay_Ole.csv"

raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")

# Now we construct a new data frame with the good data
# For this, we remove the odd columns and rows which contain dummy data

data = data.frame()

# Define the dimension of the cured dataset
n_row = nrow(raw_data) / 2 
n_col = ( ncol(raw_data) - 1 ) / 2

for(i in 1:n_row) {
	for(j in 1:n_col) { 
		data[i,j] = raw_data[2*i-1, 2*j]
	}
}

print(data)

# Use this to define names and concentrations of the compounds.
# The first value is the MAXIMUM concentration used.
concentration1 = c(128) 
for(i in 2:n_row) {concentration1[i] = concentration1[i-1]/2}

concentration2 = c(128)
for(i in 2:n_row) {concentration2[i] = concentration2[i-1]/2}

concentration3 = c(128)
for(i in 2:n_row) {concentration3[i] = concentration3[i-1]/2}

concentration4 = c(128)
for(i in 2:n_row) {concentration4[i] = concentration4[i-1]/2}

concentration = rbind(concentration1, concentration2, concentration3, concentration4)

comp_short_name = c("Comp1", "Comp2", "Comp3", "Control")
comp_full_name = c("Compound1", "Compound2", "Compound3", "Control")

# We have n_m measurements for each of the n_c compounds
n_m = 3 

# Create a new dataframe separating the measurements for each compound
# The first column has the compound name, and the second column is the measurement
compounds = data.frame()
for(j in 1:n_col) {		
	# Obtain the index "j_c" for each compound
	j_c = as.integer((j-1)/n_m) + 1
	for(i in 1:n_row) {
		name_c = paste(comp_short_name[j_c], as.character(concentration[j_c,i]), sep = "_")
		index = (j-1)*n_row + i
		compounds[index, 1] = name_c
		compounds[index, 2] = data[i,j]
		compounds[index, 3] = comp_full_name[j_c]
		#compounds[index, 4] = as.double(concentration[j_c,i])
		if (comp_full_name[j_c]=="Control" & concentration[j_c,i]==64) {
		    compounds[index, 4] = as.numeric(2^(-10)) # positive control
		} else if (comp_full_name[j_c]=="Control" & concentration[j_c,i]==128) {
		    compounds[index, 4] = as.numeric(2^10) # negative control
		} else {
		    compounds[index, 4] = as.numeric(concentration[j_c,i])
		}
	}
}


#Define the column names
column_names = c('Name', 'Fluorescence', 'Compound', 'Concentration')
colnames(compounds) = column_names

#Filter the empty wells
min_fluo = 1000
df <- subset(compounds, Fluorescence >= min_fluo, select=c(Name, Fluorescence, Compound, Concentration))

# normalize the fluorescense to the mean of the positive control
pos_con = subset(compounds, Concentration == 0.5, select=c(Fluorescence))
df$Fluorescence <- df$Fluorescence / mean(pos_con$Fluorescence)

# obtain the value of the negative control mean for fittings
neg_con = subset(compounds, Concentration == 256, select=c(Fluorescence))
min = mean(neg_con$Fluorescence) / mean(pos_con$Fluorescence)

#===============EC50===================================


#estimation of the EC50 and plot for the titration curve
library("drc")

Tabla_Compound1<-subset(df, Compound=="Compound1", select=c(Name, Fluorescence, Concentration, Compound))
Modelo_Comp1<- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound1, 
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, min, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 1")
print(summary(Modelo_Comp1))
print(confint(Modelo_Comp1, level = 0.90))
EC90Comp1<-ED(Modelo_Comp1,c(10,50,90),interval="delta")

plot(Modelo_Comp1,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText= "Compound 1", xlab = "Compound concetration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp1, col= "blue", add=TRUE, type='confidence')


Tabla_Compound2<-subset(df, Compound=="Compound2", select=c(Name, Fluorescence, Concentration, Compound))
#remove outlier
Tabla_Compound2[12,2] = NA
Modelo_Comp2 <- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound2,  
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, min, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 2")
print(summary(Modelo_Comp2))
print(confint(Modelo_Comp2, level = 0.90))
EC90Comp2<-ED(Modelo_Comp2,c(10,50,90),interval="delta")

plot(Modelo_Comp2,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText = "Compound 2", xlab = "Compound concentration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp2, col= "blue", add=TRUE, type='confidence')


Tabla_Compound3<-subset(df, Compound=="Compound3", select=c(Name, Fluorescence, Concentration, Compound))
Modelo_Comp3 <- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound3, 
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, min, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 3")
print(summary(Modelo_Comp3))
print(confint(Modelo_Comp3, level = 0.90))
EC90Comp3<-ED(Modelo_Comp3,c(10,50,90),interval="delta")

plot(Modelo_Comp3,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText = "Compound 3", xlab = "Compound concentration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp3, col= "blue", add=TRUE, type='confidence')


#--------------------------------------------------------------------------------------------------------
# Make a plot of all curves and save it

# Open file
pdf("Titration_curve_generic.pdf")
xlabel = expression(paste("Compound concentration (", mu, "g/ml)"))
# Plot
plot(Modelo_Comp1,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8,
       xlab = xlabel, ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp1, col= "blue", add=TRUE, type='confidence')

plot(Modelo_Comp2, add=TRUE, type='all', col="red", pch=4, lwd=1, cex.axis=0.8)
plot(Modelo_Comp2, col= "red", add=TRUE, type='confidence')

plot(Modelo_Comp3, add=TRUE, type='all', col="green", pch=2, lwd=1, cex.axis=0.8)
plot(Modelo_Comp3, col= "green", add=TRUE, type='confidence')

legend("topright", legend=c("Compound 1", "Compound 2 ", "Compound 3"),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)
# Close file
dev.off()




