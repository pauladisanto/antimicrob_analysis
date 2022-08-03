#levantar planilla sinergia
filename = "Primera_prueba_Sinergia_E_coli_Kanamicina_Gentamicina_resazurina.csv"
 raw_data = read.table(filename, skip=14, header=FALSE, sep="\t")


column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')
colnames(raw_data) = column_names

data <- subset( raw_data, select = -c(X))

n_m = 12 
n_b = 1    # amount of wells for each bacteria
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
		name_b = paste('Concentration', as.character(i_b), sep = "_")
		compounds[n_m*(i_c-1)+j_c,3] = name_b
	}


}

# Add a header
column_names = c('Condition', 'response', 'Condicion')
colnames(compounds) = column_names


concentration1 = c(16)
for(i in 2:8) {concentration1[i] = concentration1[i-1]/2}

concentration2 = c(32) 
for(i in 2:12) {concentration2[i] = concentration2[i-1]/2}



concentration2 = matrix(concentration2)

concentrationB = rbind(concentration2, concentration2, concentration2, concentration2,concentration2, concentration2, concentration2, concentration2)


concentrationA = rbind(concentration1,concentration1,concentration1,concentration1,concentration1,concentration1,concentration1,concentration1,concentration1,concentration1,concentration1, concentration1)

new_vector <- c(concentrationA)

 concentrationA = matrix(concentrationA)


#add column called 'concentration'
compounds$conc_c <- concentrationA
compounds$conc_r <- rev(concentrationB)


#reemplazo los valores de las concentraciones minimas que debería ser cero tengo que hacerlo de una forma mas general

compounds["conc_r"][compounds["conc_r"] == 0.015625]  <- 0
compounds["conc_c"][compounds["conc_c"] == 0.125]  <- 0



#add columns called block_id drug_row  drug_col y conc_r_unit conc_c_unit
#compounds$drug_col <- "Ampicilin"
compounds$drug_col <- "Kanamycin"
compounds$drug_row <- "Gentamicin"
compounds$conc_r_unit<-"ug/ml"
compounds$conc_c_unit<-"ug/ml"
compounds$block_id<-1
compounds$Seleccion<-"Seleccion"




#expression(paste(("\U00B5"g/"\U00B5"l))


#block_id  drug_row  drug_col conc_r  conc_c  response conc_r_unit conc_c_unit

#Tabla_Analisis<-subset(compounds, Seleccion=="Seleccion", select=c(block_id , drug_row,  drug_col, conc_r,  conc_c,  response))

Tabla_Analisis<-subset(compounds, Seleccion=="Seleccion", select=c(block_id , drug_row,  drug_col, conc_r,  conc_c,  response, conc_r_unit, conc_c_unit))

#Tabla_Analisis$response = Tabla_Analisis$response / 1.897
#Tabla_Analisis$response = Tabla_Analisis$response / 2.167
#Tabla_Analisis$response = Tabla_Analisis$response /2.59
#Tabla_Analisis$response = Tabla_Analisis$response /2.307


#Tabla_Analisis$response = Tabla_Analisis$response / 74955
#Tabla_Analisis$response = Tabla_Analisis$response / 112075
#Tabla_Analisis$response = Tabla_Analisis$response /69495
Tabla_Analisis$response = Tabla_Analisis$response /70613


Tabla_Analisis$response = Tabla_Analisis$response *100





####################analisis_interaccion#################################

library(synergyfinder)
res <- ReshapeData(
  data = Tabla_Analisis,
  data_type = "viability",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  seed = 1)

  str(res)


  res <- CalculateSynergy(
  data = res,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  Emin = NA,
  Emax = NA,
  correct_baseline = "non")


  res$drug_pairs

  str(res$synergy_scores)


  #The output adds one data frame “synergy_scores” to the input R object. It contains:

  #“block_id” - The identifier for drug combination blocks;
   # “concX” - The concentration for combined drugs;
   # “ZIP_ref”, “HSA_ref”, “Bliss_ref”, and “Loewe_ref” - The reference additive effects calculated from corresponding models;
   # “ZIP_synergy”, “HSA_synergy”, “Bliss_synergy” and "Loewe_synergy - The synergy scores calculated from corresponding models;
   # “ZIP_fit” - The fitted %inhibition for combinations which is used to calculate the ZIP synergy score.

   #The mean of synergy scores for the combination matrix (excluding the monotherapy observations) and the p-values for them are 
   #added to the “drug_pairs” table in the input R object.
#==============================================================================

  res <- CalculateSensitivity(
  data = res,
  correct_baseline = "non"
  )



   sensitive_columns <- c(
  "block_id", "drug1", "drug2",
  "ic50_1", "ic50_2",
  "ri_1", "ri_2",
  "css1_ic502", "css2_ic501", "css")
res$drug_pairs[, sensitive_columns]


#Visualization

# Dose-response curve


for (i in 1){
  PlotDoseResponseCurve(
    data = res,
    plot_block = 1,
    drug_index = i,
    plot_new = FALSE,
    record_plot = FALSE
  )
}




for (i in 2){
  PlotDoseResponseCurve(
    data = res,
    plot_block = 1,
    drug_index = i,
    plot_new = FALSE,
    record_plot = FALSE
  )
}


#####grafica la titulacion de los compuestos solos y el porcentaje de inhibicion


PlotDoseResponse(
res,
block_ids = c(1),
drugs = c(1, 2),
adjusted = TRUE,
statistic = NULL,
summary_statistic = "mean",
high_value_color = "#A90217",
low_value_color = "#2166AC",
point_color = "#C24B40",
curve_color = "black",
curve_ylim = NULL,
curve_grid = TRUE,
text_size_scale = 1,
heatmap_text_label_size_scale = 1,
heatmap_text_label_color = "#000000",
heatmap_color_range = NULL,
curve_plot_title = NULL,
heatmap_plot_title = NULL,
Emin = NA,
Emax = NA,
save_file = TRUE,
file_type = "pdf",
file_name = NULL,
width = 12,
height = 6
)





#The PlotDoseResponseCurve function will plot the dose-response curve fitted by 4-parameter log logistic function for selected drug. 
#Important parameters are:

    #data: The R object generated by function ReshapeData;
    #plot_block: Select the block from which to extract the monotherapy dose-response data;
    #drug_indes: Select the drug from the combination to draw the dose-response curve. Available values are 1, 2, 3, … 
    #to indicates the drug 1, drug 2, drug 3, … in the combination.

#This function will return an object of class recordedplot. User could use replayPlot to plot it. 
#User could modify the parameter plot_setting to control the themes of some items in the plot.


#Two-drug combination visualization
# Heatmap

#If dynamic is set as FALSE, this function will return a ggplot object. User could use + theme() to adjust the plot theme.




A = c("ZIP", "HSA", "Bliss", "Loewe")


for (i in c(1, 2,3,4)){
PlotSynergy(
res,
type = "3D",
method = A[i],
block_ids = c(1),
drugs = c(1, 2),
row_range = NULL,
col_range = NULL,
color_range = NULL,
z_range = NULL,
axis_line = FALSE,
statistic = NULL,
summary_statistic = "mean",
plot_title = NULL,
interpolate_len = 3,
high_value_color = "#A90217",
low_value_color = "#2166AC",
text_size_scale = 1,
heatmap_text_label_size_scale = 1,
heatmap_text_label_color = "#000000",
grid = TRUE,
dynamic = FALSE,
display = TRUE,
save_file = TRUE,
file_type = "pdf",
file_name = NULL,
file_path = NULL,
height = 6,
width = 6,
units = "in"
)
}

#c("ZIP", "HSA", "Bliss", "Loewe")
#type = "3D", or 2D

#####################otra forma de hacer los plots################################
Plot2DrugHeatmap(
    data = res,
    plot_block = 1,
    drugs = c(1, 2),
    plot_value = "response",
    dynamic = FALSE,
    summary_statistic = c("mean",  "median"),
    )

Plot2DrugHeatmap(
    data = res,
    plot_block = 1,
    drugs = c(1, 2),
    plot_value = "ZIP_synergy",
    dynamic = FALSE,
    summary_statistic = c( "quantile_25", "quantile_75")
    )


#2D contour plot

#If dynamic is set as FALSE, this function will return a ggplot object. User could use + theme() to adjust the plot theme.

Plot2DrugContour(
    data = res,
    plot_block = 1,
    drugs = c(1, 2),
    plot_value = "response",
    dynamic = FALSE,
    summary_statistic = c("mean", "median")
  )

Plot2DrugContour(
    data = res,
    plot_block = 1,
    drugs = c(1, 2),
    plot_value = "ZIP_synergy",
    dynamic = FALSE,
    summary_statistic = c("quantile_25", "quantile_75")
  )


Plot2DrugSurface(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "response",
  dynamic = FALSE,
  summary_statistic = c("mean", "quantile_25", "median", "quantile_75")
)


###############################################

################################################
Plot2DrugSurface(
  data = res,
  plot_block = 1,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean", "quantile_25", "median", "quantile_75")
)




#############################plot de todos los modelos juntos##########################


PlotMultiDrugBar(
  data = res,
  plot_block = 1,
  plot_value = c("response", "ZIP_synergy", "Loewe_synergy", "HSA_synergy", "Bliss_synergy"),
  sort_by = "response",
  #highlight_row = c(9.7656, 50),
  highlight_label_size = 8
)



####################plot de un modelos dos combinaciones de drogas######################

PlotSensitivitySynergy(
  data = res,
  plot_synergy = "ZIP",
  show_labels = TRUE,
  dynamic = FALSE
)