
#===================================================================================================

  filename = "Citotoxicidad_Compuestos_35_32_34_HepG2.csv"
  # We read the raw data from the file, skipping the first 14 rows which contain irrelevant information
  single_data_raw = read.table(filename, skip=14, header=FALSE, sep="\t")


  # Dummy column that is later removed
  column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')

  # Define the column names
  colnames(single_data_raw) = column_names

  single_data_raw <- subset(single_data_raw, select = -c(X))

  neg_con = 4200

  single_data_raw <-single_data_raw - neg_con

  #normalize forn the row number 7 (lowest concentration of compound)
  lowest_concec_c1 = single_data_raw[7,1:3]
  lowest_concec_c2 = single_data_raw[7,4:6]  
  lowest_concec_c3 = single_data_raw[7,7:9] 
  lowest_concec_c1 = unlist(lowest_concec_c1,recursive = F, use.names = F)    
  lowest_concec_c2 = unlist(lowest_concec_c2,recursive = F, use.names = F)    
  lowest_concec_c3 = unlist(lowest_concec_c3,recursive = F, use.names = F)    
  mean_lowest_concec_c1 = mean(lowest_concec_c1)
  mean_lowest_concec_c2 = mean(lowest_concec_c2)
  mean_lowest_concec_c3 = mean(lowest_concec_c3)

  single_data_raw[1:7,1:3] = single_data_raw[1:7,1:3] / mean_lowest_concec_c1
  single_data_raw[1:7,4:6] = single_data_raw[1:7,4:6] / mean_lowest_concec_c2
  single_data_raw[1:7,7:9] = single_data_raw[1:7,7:9] / mean_lowest_concec_c3

n_row=7
n_col=12

# Use this to define names and concentrations of the compounds.
# The first value is the MAXIMUM concentration used.
concentration1 = c(0.12) 
for(i in 2:n_row) {concentration1[i] = concentration1[i-1]/2}

concentration2 = c(0.054)
for(i in 2:n_row) {concentration2[i] = concentration2[i-1]/2}

concentration3 = c(0.36)
for(i in 2:n_row) {concentration3[i] = concentration3[i-1]/2}

concentration4 = c(12000)
for(i in 2:n_row) {concentration4[i] = concentration4[i-1]/2}

concentration = rbind(concentration1, concentration2, concentration3, concentration4)

comp_short_name = c("Comp1", "Comp2", "Comp3", "Control")
comp_full_name = c("Compound1", "Compound2", "Compound3", "Control")



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




# normalize the fluorescense to the mean of the positive control
#pos_con = subset(compounds, Concentration == concentration4[7], select=c(Fluorescence))
#compounds$Fluorescence <- compounds$Fluorescence / mean(pos_con$Fluorescence)


# obtain the value of the negative control mean for fittings
#neg_con = 4441/mean(pos_con$Fluorescence)

#compounds$Fluorescence <- compounds$Fluorescence - neg_con


#================================================================

library("drc")

Tabla_Compound1<-subset(compounds, Compound=="Compound1", select=c(Name, Fluorescence, Concentration, Compound))

Modelo_Comp1<- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound1, 
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, 0, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 1")
print(summary(Modelo_Comp1))
print(confint(Modelo_Comp1, level = 0.90))
EC90Comp1<-ED(Modelo_Comp1,c(10,50,90),interval="delta")

plot(Modelo_Comp1,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText= "Compound 1", xlab = "Compound concetration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp1, col= "blue", add=TRUE, type='confidence')


Tabla_Compound2<-subset(compounds, Compound=="Compound2", select=c(Name, Fluorescence, Concentration, Compound))
#remove outlier
#Tabla_Compound2[12,2] = NA
Modelo_Comp2 <- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound2,  
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, 0, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 2")
print(summary(Modelo_Comp2))
print(confint(Modelo_Comp2, level = 0.90))
EC90Comp2<-ED(Modelo_Comp2,c(10,50,90),interval="delta")

plot(Modelo_Comp2,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText = "Compound 2", xlab = "Compound concentration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp2, col= "blue", add=TRUE, type='confidence')


Tabla_Compound3<-subset(compounds, Compound=="Compound3", select=c(Name, Fluorescence, Concentration, Compound))
Modelo_Comp3 <- drm(Fluorescence ~ Concentration,
                   data = Tabla_Compound3, 
                   robust = 'mean', 
                   fct = LL.4 (fixed = c(NA, 0, 1, NA), names = c("Slope", "Lower limit", "Upper limit", "EC50")))

print("Compound 3")
print(summary(Modelo_Comp3))
print(confint(Modelo_Comp3, level = 0.90))
EC90Comp3<-ED(Modelo_Comp3,c(10,50,90),interval="delta")

plot(Modelo_Comp3,  type='all', col="blue", pch=16, lwd=1, cex.axis=0.8, legend = TRUE,
       legendText = "Compound 3", xlab = "Compound concentration", ylab = "Mortality", main = "Estimation EC50")
plot(Modelo_Comp3, col= "blue", add=TRUE, type='confidence')


##################################################################
#plot junto


# Open file
pdf("Cytotoxicity_generic_poster2.pdf")
xlabel = expression(paste("Compound concentration (", mu, "g/",mu, "l)"))
# Plot
plot(Modelo_Comp1,  type='all', col="#003300", pch=4, lwd=1, cex.axis=0.8,
       xlab = xlabel, ylab = "Normalized viability",  xlim=c(0.005, 1), ylim=c(0.005, 1.1), main = "Cytotoxicity HepG2")
plot(Modelo_Comp1, col= "#003300", add=TRUE, type='confidence')

plot(Modelo_Comp3, add=TRUE, type='all', col="#660066", pch=16, lwd=1, cex.axis=0.8)
plot(Modelo_Comp3, col= "#660066", add=TRUE, type='confidence')

#plot(Modelo_Comp3, add=TRUE, type='all', col="green", pch=2, lwd=1, cex.axis=0.8)
#plot(Modelo_Comp3, col= "green", add=TRUE, type='confidence')

legend("topright", legend=c("JM MIS", "JM MIR"),
       col=c("#660066", "#003300"), lty=1:2, cex=0.8)
# Close file
dev.off()
















#files[1] = "Citotoxicidad_Compuestos_35_32_34_HepG2.csv"
#files[2] = "Citotoxicidad_Compuestos_75_69_64_HepG2.csv"
#files[3] = "Citotoxicidad_Compuestos_76_66_41_HepG2.csv"
#files[4] = "Citotoxicidad_Compuestos_38_39_33_HepG2.csv"
#files[5] = "Citotoxicidad_Compuestos_36_40_15_HepG2_Plato1.csv"
#files[6] = "Citotoxicidad_Compuestos_108_104_105_HepG2_Plato2.csv"
#files[7] = "Citotoxicidad_Compuestos_10_16_21_HepG2_Plato1.csv"
#files[8] = "Citotoxicidad_Compuestos_43_48_49_HepG2_Plato2.csv"




#concentrations[[1]] = c(0.12, 0.054, 0.36,450)
#concentrations[[2]] = c(0.03, 0.324, 0.276,450)
#concentrations[[3]] = c(0.21, 0.348, 0.3,450)
#concentrations[[4]] = c(0.378, 0.348, 0.552,450)
#concentrations[[5]] = c(0.264, 0.384, 0.864,450)
#concentrations[[6]] = c(0.138, 0.198, 0.138,450)
#concentrations[[7]] = c(0.12, 0.06, 0.156,450)
#concentrations[[8]] = c(0.066, 0.066, 0.048,450)






















