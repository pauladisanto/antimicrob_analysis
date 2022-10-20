library(growthrates)
library(stringr)

filename = "CURVA_CRECIMIENTO.csv"
raw_data = read.table(filename, skip=2, header=TRUE, sep=",")


colnames(raw_data) <- paste("Sample", colnames(raw_data), sep = "_")

#I remove the first two coloums tan have the time and a the blank [,3:] 
raw_data = raw_data[,3:56]


#To replace the time in hours for minutes and I add the coloumns that the library needs (conc, replicate, S_meliloti)
#IMP I think if you have two or three replicates in all the wells you can change replicate= 2

time = seq(0, 2880, by=30) 
conc = "C"
replicate = 1
strain = "S_meliloti"

raw_data <- cbind(time, raw_data, conc, replicate, strain)

##################################################################################################
#Prepare the vector to filter the negative and positive controls
#In that case depends wich wells you want to discard

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


#I and the positive control/s at the last part of the vector 
Sample = append(Sample,"Sample_X52")


raw_data_analisys <- subset(raw_data, select = c("time", Sample, "conc", "replicate", "strain"))
######################################################################################################
#I will prepare the data set to plot the models. You will have a data set for each well
raw_data_analisys_0 = list()
for(i in 1:length(Sample)) {
raw_data_analisys_0[[i]] = subset(raw_data_analisys, select = c("time", Sample[i], "conc", "replicate", "strain"))

}
#####################################################################################################
#Here you will set the concentrations you have in the data set 
#Define the vector concentration for linear range 
concentration1 = c(130)
for(i in 2:8) {concentration1[i] = concentration1[i-1]-10} #ArgCl10

concentration2 = c(40)
for(i in 2:8) {concentration2[i] = concentration2[i-1]-5} #ArgCl12

concentration3 = c(65)
for(i in 2:8) {concentration3[i] = concentration3[i-1]-5} #C1

#Define the vector concentration for two-fold range 

concentration4 = c(32)
for(i in 2:8) {concentration4[i] = concentration4[i-1]/2} #Ampicilina

concentration5 = c(32)
for(i in 2:8) {concentration5[i] = concentration5[i-1]/2} #Kanamicina

concentration = c(concentration1,concentration2,concentration3,concentration4,concentration5)


################################################################################################
#If you have to iterate the vectors of concentration
#concentration = c()
#for (i in c(1:8)){
 #  concentration = cbind(concentration, concentration1[i])
  # concentration = cbind(concentration, concentration2[i])
   #concentration = cbind(concentration, concentration3[i])
   #concentration = cbind(concentration, concentration4[i])
   #concentration = cbind(concentration, concentration5[i])

#}
###############################################################################################################
# Another way to define the vectors for concentration
#concentration= c(concentration1,concentration2,concentration3)
     
 #   for (i in c(1:length(Sample))){ 
   #         index = 1 + ((i+2)%%3)*8 + as.integer((i-1)/3) 
  #          raw_data_analisys_0[[i]]$conc[raw_data_analisys_0[[i]]$conc == "C"] <-concentration[index]
          #print(index)
   # }
##############################################################################################################

#I put a random number for the positive control in that case this random number is 2
concentration = append(concentration, 2)
concentration = as.numeric(concentration)

##############################################################################################################
#I paste the concentration vector in the data set of each well. Each data set corresponds to 1 well

   for (i in c(1:length(Sample))){ 
          raw_data_analisys_0[[i]]$conc[raw_data_analisys_0[[i]]$conc == "C"] <-concentration[i]
          #print(index)
  }
###############################################################################################################
#this is the first model that generates the "p-vector" p = c(). The p vector has elements with yO, mumax and K

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
#this is the second model that fits the parameters of the first model to adjust the growth curves

#It is better to define the lower and the upper limits, the lower is easy but the upper you need to check your data 
#and define those parameters according to your data, for example your positive control 
lower = c(y0 = 0, mumax = 0, K = 0)
upper = c(y0 = 0.1, mumax = 0.05, K = 3)

final_model = list()
for (i in c(1:length(Sample))){ 
pdf(paste0("Growth_curve_",i, ".pdf"))   
final_model[[i]] <- fit_growthmodel(FUN = grow_logistic, p = p[[i]], dat[[i]]$time, dat[[i]]$Sample,
lower = lower, upper = upper)
plot(final_model[[i]])
dev.off()
}

###################################################################################################
#After getting all the parameters of the second model I select the first element of the list and I store them in coeff 
coeff = list()
for (i in c(1:length(Sample))){ 
coeff[[i]]= coef(final_model[[i]])
}
####################################################################################################
#this is a simple calculation to estimate the yield of each well normalized by the positive control and the % of inhibition
l = c()
for (i in c(1:length(Sample))){ 
l[i] = coeff[[c(i, 3)]]  #selects the elemnte three of the first element (in this case is K)
}

Yield = (l/l[41])*100

Inhibition = (1-(l/l[41]))*100













###################################################



      
div <- aa$K/dmso$K






