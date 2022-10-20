normalize_table <- function(df) {
  # This function normalizes a dataframe using the minimum value from the last row
  n_c = 4
  lowest_concec_c = data.frame()
  data_norm = data.frame()
  for(j_c in 1:n_c) {
    j_min = 3*j_c - 2
    j_max = 3*j_c
    # Use the average of the last row as the minimum value
    lowest_concec_c = df[7,j_min:j_max]
    lowest_concec_c = unlist(lowest_concec_c,recursive = F, use.names = F)    
    mean_lowest_concec_c = mean(lowest_concec_c)
    for(i in 1:8) {
      data_norm[i,j_min:j_max] = df[i,j_min:j_max] / mean_lowest_concec_c
    }
  }
  return(data_norm)
}

#======================================================================================

process_table <- function(filename) {  
  # This function processes a table to obtain a curated dataframe

  # We read the raw data from the file, skipping the first 14 rows which contain irrelevant information
  single_data_raw = read.table(filename, skip=14, header=FALSE, sep="\t")


  # Dummy column that is later removed
  column_names = c('X','A', 'B', 'C', 'D','E', 'F','G','H','I','J','K','L')

  # Define the column names
  colnames(single_data_raw) = column_names

  single_data_raw <- subset( single_data_raw, select = -c(X))

  # Now we normalize the data
  single_data = normalize_table(single_data_raw)
  return(single_data)
}

#======================================================================================

df_for_plot <- function(df, max_concentrations, i_table) {  
  # This function processes a dataframe "df" to make it useful for plotting. 
  # "max_concentrations" is a vector with the maximum concentrations used 
  # "i_table" is the number of the table (to keep track of the compound names)

  #row 8 has the negative control
  n_row = 7
  n_col = 12

  concentration1 = c(max_concentrations[1]) 
  for(i in 2:n_row) {concentration1[i] = concentration1[i-1]/2}

  concentration2 = c(max_concentrations[2]) 
  for(i in 2:n_row) {concentration2[i] = concentration2[i-1]/2}

  concentration3 = c(max_concentrations[3]) 
  for(i in 2:n_row) {concentration3[i] = concentration3[i-1]/2}

  concentration4 = c(max_concentrations[4]) 
  for(i in 2:n_row) {concentration4[i] = concentration4[i-1]/2}

  concentration = rbind(concentration1, concentration2, concentration3, concentration4)

  comp_short_name = c()
  comp_full_name = c()
  for(i in 1:4) {
    idx = i + 4*(i_table-1) # index
    comp_short_name[i] = paste("Comp",idx,sep="")
    comp_full_name[i] = paste("Compound",idx,sep="")
  }

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
      compounds[index, 2] = df[i,j]
      compounds[index, 3] = comp_full_name[j_c]
      compounds[index, 4] = as.numeric(concentration[j_c,i])
      }
  }

  #Define the column names
  column_names = c('Name', 'Fluorescence', 'Compound', 'Concentration')
  colnames(compounds) = column_names

  #normalize for the lowest concentration of each compound
  return(compounds)
}

#======================================================================================



plot_data <- function(df, i_table, nombre) {  
  # This function takes a dataframe "df" with 4 compounds from table number "i_table"
  # and makes a plot and saves it as a pdf file
  # PASARLE UN VECTOR DE TITULOS COMO ARGUMENTO!!!

  #vector de nombres


  for(i in 1:4) {
    idx = i + 4*(i_table-1) # index
    comp_full_name = paste("Compound",idx,sep="")
    
    nombre= c("FRI-040-1",  "FRI-050-1", "FRI-039-1", "CONTROL1" ,"FRI-057-1", "FRI-105", "FRI-034-1", "FRI-070-1", "FRI-0.71-1",  "FRI-086-1", "FRI-087-1", "CONTROL2","FRI-016-2", "FRI-094-1", "FRI-053-1","CONTROL3", "FRI-089-1", "FRI-049-1", "FRI-055-1","CONTROL4", "FRI-091-1", "FRI-090-1", "FRI-0.67-1", "CONTROL5" ,"FRI-016-1", "FRI-051-1", "FRI-054-1", "CONTROL6","KAT-C-006-H1F1",  "FRI-058-1", "REPLICAS1", "CONTROLCONSDS")
 

    Tabla_Compound <- subset(df, Compound==comp_full_name, select=c(Name, Fluorescence, Concentration, Compound))
    
    print(Tabla_Compound)


    mean_fluorescence <- tapply(Tabla_Compound$Fluorescence, factor(Tabla_Compound$Concentration),mean,na.rm=TRUE)

    sd_fluorescence <- tapply(Tabla_Compound$Fluorescence, factor(Tabla_Compound$Concentration),sd, na.rm=TRUE)

    concentration <- Tabla_Compound$Concentration # CHECK!!!!!!!

    Tabla <- data.frame(mean_fluo=mean_fluorescence, sd_fluo=sd_fluorescence, concentration=rev(concentration))

    # Make a plot of all curves and save it

    # Default line plot
    p <- ggplot(Tabla, aes(x=concentration, y=mean_fluo)) + 
      geom_line(linetype = "dashed") +
      geom_point()+
      geom_pointrange(aes(ymin=mean_fluo-sd_fluo, ymax=mean_fluo+sd_fluo), width=0.1)+
      scale_x_log10()+
      ylim(0, 2)
    #  print(p)

    # TO DO > PASARLE UN VECTOR DE TITULOS!!!
    p + labs(title=paste("Cytotoxicity MCF-7", nombre[[idx]], sep="_"), x= "concentration (\u00b5g/\u00b5l)", y = "Relative Fluorescence") + theme_classic() + theme(text = element_text(size = 20))

    file_name = nombre[[idx]]
    file_name = paste(file_name, ".pdf", sep="")
    ggsave(file_name)
  }

}



  
#======================================================================================


library("ggplot2")


files = list()

files[1] = "Citotoxicidad_plato1_compuestos1a3_MCF7.csv"
files[2] = "Citotoxicidad_plato2_compuestos4a7_MCF7.csv"
files[3] = "Citotoxicidad_plato1_compuestos8a10_MCF7.csv"
files[4] = "Citotoxicidad_plato2_compuestos11a13_MCF7.csv"
files[5] = "Citotoxicidad_plato3_compuestos14a16_MCF7.csv"
files[6] = "Citotoxicidad_plato1_compuestos17a19_MCF7.csv"
files[7] = "Citotoxicidad_plato2_compuestos20a22_MCF7.csv"
files[8] = "Citotoxicidad_plato1_compuestos23_24_y-repeticiones_MCF7.csv"


concentrations = list()
concentrations[[1]] = c(0.048, 0.0912, 0.0444, 450)
concentrations[[2]] = c(0.084, 0.1308, 0.1488, 0.0648)
concentrations[[3]] = c(0.09,0.036,0.1428,450)
concentrations[[4]] = c(0.1236,0.1212,0.0756,450)
concentrations[[5]] = c(0.0552,0.06,0.0288,450)
concentrations[[6]] = c(0.2052,0.1296,0.09,450)
concentrations[[7]] = c(0.1128,0.054,0.0648,450)
concentrations[[8]] = c(0.1044,0.0756,700,450)

# process the tables and merge all individual data frames into a complete one
for(i in 1:8) {
  single_data = process_table(files[[i]])
  transformed_data = df_for_plot(single_data, concentrations[[i]], i)
  plot_data(transformed_data, i)
}
                                          
