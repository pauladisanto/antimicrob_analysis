
filename = "MIC_P_putida_comp_catarina_1_3_8_17_fluorescencia.csv"
# We read the raw data from the file, skipping the first 14 rows which contain irrelevant information
single_data_raw = read.table(filename, skip=14, header=FALSE, sep="\t")

# Dummy column that is later removed
column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
# Define the column names
colnames(single_data_raw) = column_names

single_data_raw <- subset(single_data_raw, select = -c(X))


#select the wells of the data frame and convert it into a vector and estimate the mean
negative = mean(as.numeric(as.vector(single_data_raw[8,1:3])))
positive = mean(as.numeric(as.vector(single_data_raw[8,4:6])))

# Use this to define names and concentrations of the compounds.
# The first value is the MAXIMUM concentration used.

n_row=7
concentration1 = c(160) 
for(i in 2:n_row) {concentration1[i] = concentration1[i-1]/2}

concentration2 = c(160)
for(i in 2:n_row) {concentration2[i] = concentration2[i-1]/2}

concentration3 = c(160)
for(i in 2:n_row) {concentration3[i] = concentration3[i-1]/2}

concentration4 = c(160)
for(i in 2:n_row) {concentration4[i] = concentration4[i-1]/2}

concentration = rbind(concentration1, concentration2, concentration3, concentration4)


comp_short_name = c("Comp1", "Comp2", "Comp3", "Comp4")
comp_full_name = c("Compound1", "Compound2", "Compound3", "Compound4")


n_col=12
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
      compounds[index, 2] = single_data_raw[i,j]
      compounds[index, 3] = comp_full_name[j_c]
      compounds[index, 4] = as.numeric(concentration[j_c,i])
      }
  }

#Define the column names
column_names = c('Name', 'Fluorescence', 'Compound', 'Concentration')
colnames(compounds) = column_names


#normalize the data by the positive value
compounds$Fluorescence = ((compounds$Fluorescence / positive) - (negative/ positive))



#================================================

library("drc")


#Selection of the subtables according to the Compound
compound=c("Compound1", "Compound2", "Compound3", "Compound4")
Tabla_Compound = list()
Modelo_Comp = list()
EC_10_5_90 = list()

for(i in 1:4) {

Tabla_Compound[[i]]<-subset(compounds, Compound==compound[i], select=c(Name, Fluorescence, Concentration, Compound))



Modelo_Comp[[i]] <- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound[[i]],  
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, 0, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))



#print("Compound ")
print(summary(Modelo_Comp[[i]]))
print(confint(Modelo_Comp[[i]], level = 0.90))
EC_10_5_90[[i]]<-ED(Modelo_Comp[[i]],c(10,50,90),interval="delta")



# Open file
pdf(paste("Titration_curve_planktonic", i, ".pdf"))
xlabel = expression(paste("Compound concentration (", mu, "g/ml)"))


plot(Modelo_Comp[[i]],  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText= compound[i], xlab = "Compound concetration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp[[i]], col= "blue", add=TRUE, type='confidence')


#legend("topright", legend=c("Compound 1", "Compound 2 ", "Compound 3"),
       #col=c("blue", "red", "green"), lty=1:2, cex=0.8)
# Close file
dev.off()

}


#================================================



