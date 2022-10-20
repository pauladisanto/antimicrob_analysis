library(growthrates)
library(stringr)

filename = "growth_curve_1_S_meliloti_compuestos_ArgCl10_ArgCl12_C1.csv"
raw_data = read.table(filename, skip=13, header=FALSE, sep=",")

rep_str = c('Sample X'='Sample_X')
raw_data$V3 <- str_replace_all(raw_data$V3, rep_str)
print(raw_data)

raw_data = raw_data[,3:99]

#to define the fisrt row as header
raw_data <- setNames(data.frame(t(raw_data[ , - 1])), raw_data[ , 1])  # Transpose data
#eliminate the columne time
raw_data <- subset(raw_data, select = -c(Time))

#To replace the time in hours for minutes
time = seq(0, 2850, by=30) 
conc = "C"
replicate = 1
strain = "S_meliloti"

raw_data <- cbind(time, raw_data, conc, replicate, strain)

##################################################################################################
#This is because I have oredered the samples in the plate in a bad way
Sample = c()
for(i in 1:96) {
	if(i %in% seq(1, 91, by=3) ) {  
      i_c = seq(1, 91, by=3)
 Sample[i_c] = c(paste0("Sample_X", i_c))
}
else {
    next
}
}

Sample = Sample[!is.na(Sample)]

Sample = Sample[!Sample %in% c("Sample_X10", "Sample_X22","Sample_X34","Sample_X46","Sample_X58","Sample_X70","Sample_X82")] 

Sample = append(Sample,"Sample_X34")

raw_data_analisys <- subset(raw_data, select = c("time", Sample, "conc", "replicate", "strain"))
######################################################################################################
#I will prepare the data set to plot the models

raw_data_analisys_0 = list()
for(i in 1:length(Sample)) {
raw_data_analisys_0[[i]] = subset(raw_data_analisys, select = c("time", Sample[i], "conc", "replicate", "strain"))

}

#Define the vector concentration
concentration1 = c(250)
for(i in 2:8) {concentration1[i] = concentration1[i-1]/2} #compound 34

concentration2 = c(250)
for(i in 2:8) {concentration2[i] = concentration2[i-1]/2} #compound 35

concentration3 = c(125)
for(i in 2:8) {concentration3[i] = concentration3[i-1]/2} #DMSO



concentration = c()
for (i in c(1:8)){
   concentration = cbind(concentration, concentration1[i])
   concentration = cbind(concentration, concentration2[i])
   concentration = cbind(concentration, concentration3[i])
}

concentration = append(concentration, 2)
concentration = as.numeric(concentration)
     
   for (i in c(1:length(Sample))){ 
          raw_data_analisys_0[[i]]$conc[raw_data_analisys_0[[i]]$conc == "C"] <-concentration[i]
          #print(index)
  }
##############################################################
# Another way to define the vectors for concentration
#concentration= c(concentration1,concentration2,concentration3)
     
 #   for (i in c(1:length(Sample))){ 
   #         index = 1 + ((i+2)%%3)*8 + as.integer((i-1)/3) 
  #          raw_data_analisys_0[[i]]$conc[raw_data_analisys_0[[i]]$conc == "C"] <-concentration[index]
          #print(index)
   # }
##########################################################
p = c()
fit = list()
value_adjustment = c()
dat= list()
for (i in c(1:length(Sample))){ 
dat[[i]]<- raw_data_analisys_0[[i]]
dat[[i]]$Sample = as.numeric(dat[[i]]$Sample)
value_adjustment[i]= min(as.numeric(raw_data_analisys_0[[i]]$Sample))
dat[[i]]$Sample = as.numeric(dat[[i]]$Sample) - as.numeric(value_adjustment[i])+0.001
fit[[i]] <- fit_spline(y = dat[[i]]$Sample, time = dat[[i]]$time, spar = 0.5)
p[[i]] <- c(coef(fit[[i]]), K = max(as.numeric(dat[[i]]$Sample))) # create vector using fitted model
}


#########################################################################################
lower = c(y0 = 0, mumax = 0, K = 0)
upper = c(y0 = 2, mumax = 0.05, K = 15)

final_model = list()
for (i in c(1:length(Sample))){ 
pdf(paste0("Growth_curve_",i, ".pdf"))   
final_model[[i]] <- fit_growthmodel(FUN = grow_logistic, p = p[[i]], dat[[i]]$time, dat[[i]]$Sample,
lower = lower, upper = upper)
plot(final_model[[i]])
dev.off()
}

###################################################################################################
coeff = list()
for (i in c(1:length(Sample))){ 
coeff[[i]]= coef(final_model[[i]])
}
####################################################################################################

l = c()
for (i in c(1:length(Sample))){ 
l[i] = coeff[[c(i, 3)]]  #selects the elemnte three of the first element (in this case is K)
}


Yield = (l/l[25])*100

Inhibition = (1-(l/l[25]))*100


###################################################



      
