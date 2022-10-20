#https://cran.r-project.org/web/packages/growthcurver/vignettes/Growthcurver-vignette.html
library(growthcurver)
################################################################################################
#Read the input table, rename the data and eliminate the tw first columns
filename = "CURVA_CRECIMIENTO.csv"
raw_data = read.table(filename, skip=2, header=TRUE, sep=",")

colnames(raw_data) <- paste("Sample", colnames(raw_data), sep = "_")

raw_data = raw_data[,3:56]

#To replace the time in hours for minutes
time = seq(0, 2880, by=30) 
conc = "C"
replicate = 1
strain = "S_meliloti"

raw_data <- cbind(time, raw_data, conc, replicate, strain)

##################################################################################################
#Prepare the vector to filter the negative and positive controls
Sample = c()
for(i in 1:54) {
	if(i %in% seq(1, 54, by=1) ) {  
      i_c = seq(1, 54, by=1)
 Sample[i_c] = c(paste0("Sample_X", i_c))
}
else {
    next
}
}

Sample = Sample[!is.na(Sample)]

Sample = Sample[!Sample %in% c("Sample_X9", "Sample_X10","Sample_X19","Sample_X20","Sample_X29","Sample_X30","Sample_X39", "Sample_X40","Sample_X49", "Sample_X50", "Sample_X51","Sample_X52","Sample_X53","Sample_X54")] 

#add the positive control as the last sample
Sample = append(Sample,"Sample_X52")

#####################################################################################################
#Select the columns with time and Sample
raw_data_analisys2 <- subset(raw_data, select = c("time", Sample))

gc_out<- SummarizeGrowthByPlate(raw_data_analisys2)

#if you have a blank you can subtract that value 
#gc_out <- SummarizeGrowthByPlate(raw_data_analisys2, bg_correct = "blank")

gc_out <- SummarizeGrowthByPlate(raw_data_analisys2, plot_fit = TRUE, 
                                 plot_file = "gc_plots.pdf")

head(gc_out)

#Define the vector concentration for linear range 
concentration1 = c(130)
for(i in 2:8) {concentration1[i] = concentration1[i-1]-10} #ArgCl10

concentration2 = c(40)
for(i in 2:8) {concentration2[i] = concentration2[i-1]-5} #ArgCl12

concentration3 = c(65)
for(i in 2:8) {concentration3[i] = concentration3[i-1]-5} #C1

concentration4 = c(32)
for(i in 2:8) {concentration4[i] = concentration4[i-1]/2} #Ampicilina

concentration5 = c(32)
for(i in 2:8) {concentration5[i] = concentration5[i-1]/2} #Kanamicina

concentration = c(concentration1,concentration2,concentration3,concentration4,concentration5)

concentration = append(concentration, 0)
concentration = as.numeric(concentration)
     
gc_out= cbind(concentration,gc_out)


output_file_name <- "/home/paula/Back_up/pau/Desktop/GU 2021/Curva_de_crecimiento/myfilename.csv"
write.table(gc_out, file = output_file_name, 
            quote = FALSE, sep = "\t", row.names = FALSE)


#gc_out gives you the following parameters
#t_mid is the time at which the population density reaches 1/2 K (which occurs at the inflection point), 
#t_gen is the fastest possible generation time (also called the doubling time)
#r is the growth rate 
# k, n0, and r are the values of the parameters for the logistic equation that best fit the data
#I think k is the yield and n0 is equivalent to y0 of the other model that gives you the original number of cells
#######################################################################################################
raw_data_analisys <- subset(raw_data, select = c("time", Sample, "conc", "replicate", "strain"))

######################################################################################################
#I will prepare the data set to plot the models

raw_data_analisys_0 = list()
for(i in 1:length(Sample)) {
raw_data_analisys_0[[i]] = subset(raw_data_analisys, select = c("time", Sample[i], "conc", "replicate", "strain"))

}

     
   for (i in c(1:length(Sample))){ 
          raw_data_analisys_0[[i]]$conc[raw_data_analisys_0[[i]]$conc == "C"] <-concentration[i]
          #print(index)
  }

#############################################################################################
#this plot is with the value adjusted with the minimun value
gc_fit = list()
value_adjustment = c()
dat= list()
for (i in c(1:length(Sample))){ 
pdf(paste0("Growth_curve_DT_",i, ".pdf"))   
dat[[i]]<- raw_data_analisys_0[[i]]
dat[[i]]$Sample = as.numeric(dat[[i]]$Sample)
value_adjustment[i]= min(as.numeric(raw_data_analisys_0[[i]]$Sample))
dat[[i]]$Sample = as.numeric(dat[[i]]$Sample) - as.numeric(value_adjustment[i])+0.001
gc_fit[[i]] <- SummarizeGrowth(dat[[i]]$time, dat[[i]]$Sample)
plot(gc_fit[[i]])
dev.off()
}

#str_value = data.frame()
#for (i in c(1:length(Sample))){ 
#str_value[[i]] = str(gc_fit[[i]]$vals$r)
#str_value[[i]] = str(gc_fit[[i]]$vals$t_gen)
#str_value[[i]] = str(gc_fit[[i]]$vals$t_mid)
#str_value[[i]] = str(gc_fit[[i]]$vals$k)
#}





