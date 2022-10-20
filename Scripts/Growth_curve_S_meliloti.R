gdata <- read.table(file="growth_curve_1_S_meliloti.csv", sep="\t", header=T)

colnames(gdata)

#plot(gdata[,1], gdata[,2], main="Growth curve S meliloti", xlab="minutes", ylab="OD600")

#plot(gdata[,2], gdata[,5], main="OD versus Time", xlab="minutes", ylab="OD600")

# since this is a growth curve, we want to plot the log
# of the OD
#plot(gdata[,1], log(gdata[,2]), main="Growth curve S meliloti", xlab="minutes", ylab="log(OD600)")
 

for(i in 2:59) {

pdf(paste0("Growth_curve",i, ".pdf"))	
Plot_growth <- plot(gdata[,1], gdata[,i], main= paste0("Growth_curve_S meliloti_Sample_X",i) , xlab="minutes",
ylab="OD600")
dev.off()
}


#===============

for(i in 2:59) {

pdf(paste0("Growth_curve",i, ".pdf"))	
Plot_growth <- plot(gdata[,1], log(gdata[,i]), main= paste0("Growth_curve_S meliloti_Sample_X",i) , xlab="minutes",
ylab="log(OD600)")
dev.off()
}


#============================================================================= 
fit=list()
# en "fit" quedan guardados TODOS los parámetros del ajuste para cada condición
for(i in 2:58) {
	fit[i] = lm(log(gdata[3:18, i]) ~ gdata[3:18,1])
	print(i)
	#summary(fit[[i]])$coefficients[1]
	#summary(fit[[i]])$coefficients[2]
}

fit=c()
for(i in 2:58) {
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