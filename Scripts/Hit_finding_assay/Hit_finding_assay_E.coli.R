
#==============================================================================================================================
#
# FUNCTIONS (DO NOT MOVE FROM HERE)
#
#==============================================================================================================================

number_to_name <- function(n) {
	# n is the number of the internal label
	names = c('KLG-4P', 'KLG-13A', 'KLG-9D', 'KLG-6M', 'KLG-12R', 
		'KLG-12T', 'KLG-11F', 'KLG-6I', 'KLG-9M', 'JK-8H', 'KLG-10B', 
		'KLG-14F', 'KLG-12J', 'KLG-6P', 'KLG-4A (CRUDE)', 'JM-3MTS', 
		'JM-4MTS', 'JM-2MTR', 'JM-23MTL', 'JM-11MIS', 'JM-11MTL', 'JM-12MTL', 
		'JM-4MFL', 'JM-6MFL', 'JM-1MFL', 'JM-9MIS', 'JM-16MIS', 'JM-5AMIS', 'JM-5MTL', 
		'JM-1BMIS', 'JM-8MIS', 'JM-15MIS', 'JM-MTR', 'JM-MIS', 'JM-MIR', 'JM-MFS', 
		'JM-MTL', 'JM-MIL', 'JM-MFR', 'JM-MFL', 'JM-MTS', 'MAN-SML-9596A', 'MAN-SML-105I', 
		'MAN-SML-105F', 'MAN-PES-B13', 'MAN-SML-105C', 'MAN-SMR-D180a', 'MAN-SML-110I', 
		'MAN-SMR-D75-78C23A', 'MAN-SMR-D108c', 'MAN-SML-117A', 'MAN-SMR-D8284CE12', 'MAN-SMR-D101H', 
		'MAN-SMR-D1416Y', 'MAN-SML-110M', 'MAN-SMR-D7C', 'MAN-SMR-D8B2', 'PERE', 'PELE', 'PESE', 
		'SMLE', 'SMLE-SMLE', 'MAN-SML-6465A (ACN LC-M5)', 'CLL', 'CKL', 'CLS', 'CKS', 'CKR', 'CLR', 'SMRE', 
		'KVBA-38B', 'KVBA-38A', 'KVBA-38C', 'KVPA-36C', 'KVPA-36L', 'KVBA-38D', 'KVBA-38F', 'NKR-2B','JM-8A MFL',	'JM-2MFR',	'JM-23MFL',	'JM-22 MFL',	'JM-16B MFL',	'JM-17 MFL',	'JM-4 MFR',	'JM-10B-1 MFL',	'JM-20 MFL',	'JM-12 MFL',	'JM-11 MFL',	
	'JM-14 MFL',	'JM-21 MFL',	'JM-19B MFL',	'JM-06MIS (MIXTURE OF TWO)',	'ET-CKL-1',	'CK-CF12-1CC',	'ET-CKL-7',	'CK-CF11-1D',	'ET-CKL-5',	'CK-CF10-1B',	
	'ET-CKL-6',	'ET-CKL-13',	'ET-CKL-12',	'ET-CKL-Crude',	'ET-CLS-Crude',	'ET-CLL-Crude',	'ET-CKS-crude',	'ET-CKR-Crude',	'ET-CLR-Crude', 'POS', 'NEG')

	return(names[n])
}


process_table <- function(filename, compound_numbers, verbose) {  
  # This function reads a table and returns a dataframe of compound names, means and standard deviations

  # Check that the argument is correct
  if (length(compound_numbers) < 32){print("ERROR! Missing compound numbers")}
  if (length(compound_numbers) > 32){print("ERROR! Too many compound numbers")}

	# We read the raw data from the file, skipping the first 14 rows which contain irrelevant information
	raw_data = read.table(filename,	skip=14, header=FALSE, sep="\t")

	# Now we construct a new data frame with the good data
	# For this, we remove the odd columns and rows which contain dummy data

	# Define the dimension of the cured dataset
	n_row = nrow(raw_data) / 2 
	n_col = (ncol(raw_data) - 1) / 2
	data = data.frame()
	for(i in 1:n_row) {
		for(j in 1:n_col) { 
			data[i,j] = raw_data[2*i-1, 2*j]
		}
	}

	if(verbose){print(data)}

	# We have n_m measurements for each of the n_c compounds
	n_m = 3 
	n_c = n_row * (n_col/n_m)

	# Create a new dataframe separating the measurements for each compound
	# The first column has the compound name, and the second column is the measurement
	compounds = data.frame()
	sum = 0
	for(i in 1:n_row) {
		for(j in 1:n_col) {
			# Obtain the index "i_c" for each compound
			sum = sum + 1
			i_c = as.integer( (sum-1)/3 ) + 1 # this truncates the number (no decimals)
			j_c = (j-1)%%3 + 1 # assign the same value three times
			compounds[n_m*(i_c-1)+j_c,1] = number_to_name(compound_numbers[i_c]) # name of the first column  
			compounds[n_m*(i_c-1)+j_c,2] = data[i,j]
		}
	}

	if(verbose){print(compounds)}

	# Add a header with the column names
	colnames(compounds) = c('Compound', 'Fluorescence')

	# Filter small values because maybe the plate is not full
	min_fluo = 1
	Compounds_filter <- subset(compounds, Fluorescence >= min_fluo, select=c(Compound, Fluorescence))

	# Calculate the mean and standard deviation
	mean_compounds <- tapply(Compounds_filter$Fluorescence,factor(Compounds_filter$Compound),mean)
	sd_compounds <- tapply(Compounds_filter$Fluorescence,factor(Compounds_filter$Compound),sd)

	# Normalize to the positive control and remove both controls from the list
	# Also store the name of the compounds that will be plotted
	mean_control_neg = mean_compounds["NEG"]
	mean_control_pos = mean_compounds["POS"] 
	mean_compounds_ratio = c()
	sd_compounds_ratio = c()
	name_compounds_plot = c()
	i = 0
	index_name = 0
	for(mean in mean_compounds){
		i = i + 1
		index_name = index_name + 1
		if(mean != mean_control_neg && mean != mean_control_pos){
			mean_compounds_ratio[i] = mean/mean_control_pos
			sd_compounds_ratio[i] = sd_compounds[i]/mean_control_pos
			name_compounds_plot[i] = names(mean_compounds)[index_name]} else i = i-1
	}

	# Create a dataframe and return it
	df = data.frame(
		"Comp"=name_compounds_plot,
		"Mean"=mean_compounds_ratio, 
		"SD"=sd_compounds_ratio
	)
	if(verbose){print(df)}

	return(df)
}


#==============================================================================================================================
#
# COMMANDS START HERE
#
#==============================================================================================================================

library("ggplot2")

verbose = TRUE # Set to TRUE to print stuff on the terminal

# I am using 19 and 110 for positive and negative controls, respectively (MUST BE CONSISTENT WITH FUNCTION number_to_name)

filename1 = "Hit_assay_samples_15_to_69_calibrated_equipment.csv"
#compound_numbers1 = c(15, 40, 79, 33, 41, 80, 34, 64, 80, 35, 65, 36, 66, 37, 67, 38, 68, 39, 69)
# 41 is bad, assume it is a negative control (80) to filter it
compound_numbers1 = c(15, 40, 109, 110, 33, 110, 110, 110, 34, 64, 110, 110, 35, 65, 110, 110, 36, 66, 110, 110, 37, 67, 110, 110, 38, 68, 110, 110, 39, 69, 110, 110)
df1 = process_table(filename1, compound_numbers1, verbose)

filename2 = "Hit_assay_samples_1_to_41_calibrated_equipment.csv"
#compound_numbers2 = c(1, 9, 18, 26, 2, 10, 19, 27, 3, 11, 20, 28, 4, 12, 21, 29, 5, 13, 22, 41, 6, 14, 23, 79, 7, 16, 24, 79, 8, 17, 25, 80)
# 22 and 23 are bad, assume they are negative controls (80) to filter them
compound_numbers2 = c(1, 9, 18, 26, 2, 10, 19, 27, 3, 11, 20, 28, 4, 12, 21, 29, 5, 13, 110, 41, 6, 14, 110, 109, 7, 16, 24, 109, 8, 17, 25, 110)
df2 = process_table(filename2, compound_numbers2, verbose)

filename3 = "Hit_assay_samples_22_to_71_calibrated_equipment_bottom_light.csv"
compound_numbers3 = c(22, 45, 53, 61, 23, 46, 54, 62, 30, 47, 55, 63, 31, 48, 56, 70, 32, 49, 57, 71, 42, 50, 58, 109, 43, 51, 59, 109, 44, 52, 60, 110)
df3 = process_table(filename3, compound_numbers3, verbose)

filename4 = "Hit_assay_samples_72_to_21_calibrated_equipment_bottom_light_Ecoli.csv"
compound_numbers4 = c(72, 110, 110, 110, 73, 110, 110, 110, 74, 110, 110, 110, 75, 110, 110, 110, 76, 110, 110, 110, 77, 110, 110, 110, 78, 110, 110, 110, 109, 110, 110, 110)
df4 = process_table(filename4, compound_numbers4, verbose)

filename5 = "Hit_assay_samples_79_to_108_Ecoli_nuevos_compuestos.csv"
compound_numbers5 = c(79, 88, 97, 105, 80, 90, 98, 106, 81, 91, 99, 107, 82, 92, 100, 108, 83, 93, 101, 109, 84, 94, 102, 109, 85, 95, 103, 110, 87, 96, 104, 110)
df5 = process_table(filename5, compound_numbers5, verbose)


# Add all the dataframes into a single one for plotting
df <- rbind(df1, df2, df3, df4, df5)

df_ordered = df[order(df$Mean),]


# Make the barplot
p <- ggplot(data=df, aes(x=reorder(Comp, Mean), y=Mean)) +
  geom_bar(stat="identity",fill = "darkolivegreen")+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2,
                 position=position_dodge(0.05)) + 
  							 theme_bw() +
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=0.5)) +
                 theme(legend.position="top") + 
                 theme(aspect.ratio=0.25) + 
                 theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
                 labs(title=element_blank(), x="Compound", y = "Relative Fluorescence")+
                 ggtitle("Escherichia coli")


#p + xlab("New X axis labe")

print(p)


ggsave("plot_compuestos_Ecoli_Nuevos_compuestos_poster.jpg")


ggsave("plot_compuestos_Ecoli_Nuevos_compuestos_poster.pdf")



#===================================================================================================================
# ANOTHER VERSION
