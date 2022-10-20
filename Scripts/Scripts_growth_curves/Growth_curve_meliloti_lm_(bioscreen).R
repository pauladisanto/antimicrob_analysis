
filename = "CURVA_CRECIMIENTO.csv"
raw_data = read.table(filename, skip=2, header=TRUE, sep=",")


colnames(raw_data) <- paste("Sample", colnames(raw_data), sep = "_")

#this has to be setted by hand according to your data [,3:]
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

Sample = append(Sample,"Sample_X52")



gdata <- subset(raw_data, select = c("time", Sample))

##################################################################################
 

for(i in 2:length(Sample)) {

pdf(paste0("Growth_curve",i, ".pdf"))	
Plot_growth <- plot(gdata[,1], gdata[,i], main= paste0("Growth_curve_S meliloti_Sample_X",i) , xlab="minutes",
ylab="OD600")
dev.off()
}


#===============

for(i in 2:(length(Sample)+1)) {

pdf(paste0("Growth_curve",i, ".pdf"))	
Plot_growth <- plot(gdata[,1], log(gdata[,i]), main= paste0("Growth_curve_S meliloti_Sample_X",i) , xlab="minutes",
ylab="log(OD600)")
dev.off()
}


#============================================================================= 
fit=list()
# en "fit" quedan guardados TODOS los parámetros del ajuste para cada condición
for(i in 2:(length(Sample)+1)) {
	fit[i] = lm(log(gdata[3:18, i]) ~ gdata[3:18,1])
	print(i)
	#summary(fit[[i]])$coefficients[1]
	#summary(fit[[i]])$coefficients[2]
}

fit=c()
for(i in 2:(length(Sample)+1)) {
	a = lm(log(gdata[3:18, i]) ~ gdata[3:18,1])
	# the slope of the line can be found by looking at the second coefficient
	fit[i] = summary(a)$coefficients[2]  
	print(fit[i])
}

# use the slope to calculate the doubling time of the
# cells during their exponential growth phase
# the formula is log(2)/k where k is the growth rate (slope from the curve)

slopes = c()
slopes = log(2)/fit #implicit loop (note that i=1 element is empty!)



#T62nFTeGzZNyqnf
#T62nFTeGzZNyq_nf