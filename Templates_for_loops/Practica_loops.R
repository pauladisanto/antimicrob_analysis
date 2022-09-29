#For loop

for (var in vector) {
    statement(s)    
}

#======================================

# R Program to demonstrate
# the use of for loop
for (i in 1: 4)
{
    print(i ^ 2)
    print(i / 2)
    print(i * 2)
}

#======================================
# R Program to demonstrate the use of
# for loop along with concatenate
for (i in c(-8, 9, 11, 45))
{
    print(i)
}

#======================================

for (i in c(-8, 9, 11, 45))
{
    print(i+1)
}


#======================================
#este=es=igual=pero=defino=X=afuera===

x <- c(-8, 9, 11, 45)
for (i in x)
{
    print(i)
}

#======================================
# R Program to demonstrate the use of
# nested for loop
for (i in 1:3)
{
    for (j in 1:i)
    {
        print(i * j)
    }
}

#======================================
k = 6

for (i in 1:3)
{
    for (j in 1:k)
    {
        print(i + j+1)
    }
}

#======================================

# R Program to demonstrate the use of
# break in for loop
for (i in c(3, 6, 23, 19, 0, 21))
{
   if (i == 0)
   {
      break
   }
print(i)
}
print("Outside Loop")

#======================================

# R Program to demonstrate the use of
# break in for loop
for (i in c(3, 6, 23, 0, 10, 21))
{
   if (i == 10)
   {
      break
   }
print(i)
}
print("oh no")

#======================================
# R Program to demonstrate the use of
# next in for loop
for (i in c(3, 6, 23, 19, 0, 21))
{
   if (i == 0)
   {
      next
   }
   print(i)
}
print('Outside Loop')

#======================================
# R program to illustrate while loop

result <- c("Hello World")
i <- 1

# test expression
while (i < 6) {

print(result)
   
# update expression
i = i + 2
}

#======================================
# R program to illustrate while loop

result <- 1
i <- 1

# test expression
while (i < 6) {

print(result)
   
# update expression
i = i + 1
result = result + 1
}
#======================================
# R program to illustrate while loop

result <- c("Hello World")
i <- 1

# test expression
while (i < 6) {

   print(result)
   
   if( i == 3){
      break}
   # update expression
   i = i + 1
}

#======================================
# R program to illustrate repeat loop
  
result <- c("Hello World")
i <- 1
  
# test expression 
repeat {
  
   print(result)
     
   # update expression 
   i <- i + 1
     
   # Breaking condition
   if(i >6) {
      break
   }
}

#====================================== 
# R program to illustrate repeat loop
  
result <- 1
i <- 1
  
# test expression 
repeat {
  
   print(result)
     
   # update expression 
   i <- i + 1
   result = result + 1
  
   # Breaking condition
   if(i > 5) {
      break
   }
}

#======================================
#Program to check for even and odd numbers
a <- 6.5
if ((a %% 2) == 0)
{ 
    print("even") 
} else {
    print("odd")
}
#======================================
#Program to check for prime numbers # no funciona el else
a <- 8
b <- a/2
flag <- 0
i <- 2
repeat
{
   if ((a %% i)== 0)
   {
      flag <- 1
      break
   }
}

if (flag == 1)
{
   print("composite")
} else {
   print("prime")
}

#======================================
# R program for break statement in For-loop

no <- 1:10

for (val in no)
{
   if (val == 5)
   {
      print(paste("Coming out from for loop Where i = ", val))
      break
   }
   print(paste("Values are: ", val))
}
#======================================

# R Break Statement Example
a<-1  
while (a < 10)
{  
   print(a) 
   if(a==7) 
      break 
   a = a + 1   
}  

#el a=a+1 es para que no quede corriendo por siempre

#======================================

#The next statement is used to skip the current iteration in 
#the loop and move to the next iteration without exiting from the loop itself.
# R Next Statement Example

no <- 1:10

for (val in no)
{
   if (val == 6)
   {
      print(paste("Skipping for loop Where i = ", val))
      next
   }
   print(paste("Values are: ", val))
}

#if esta anidado en el for y cuando val = 6  entonces printea skipping..

#======================================
# R Next Statement Example
x <- 0
while(x < 15)
{
   x <- x + 1;
   if (x == 3)
      next;
   print(x);
}

#se imprimen todos los numeros expcepto el 3 (ya que x vale 0)

