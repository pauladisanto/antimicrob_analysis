#Functions in R Programming

# Find sum of numbers 4 to 6.
print(sum(4:6))

# Find max of numbers 4 and 6.
print(max(4:6))

# Find min of numbers 4 and 6.
print(min(4:6))


#====================================

# A simple R function to check
# whether x is even or odd

evenOdd = function(x){
if(x %% 2 == 0)
	return("even")
else
	return("odd")
}

print(evenOdd(4))
print(evenOdd(1))

#====================================

# A simple R function to calculate
# area of a circle

areaOfCircle = function(radius){
area = pi*radius^2
return(area)
}

print(areaOfCircle(5))
#el print le dice que x usar para calcular el área

#====================================

# A simple R function to calculate
# area and perimeter of a rectangle

Rectangle = function(length, width){
area = length * width
perimeter = 2 * (length + width)

# create an object called result which is
# a list of area and perimeter
result = list("Area" = area, "Perimeter" = perimeter)
return(result)
}

resultList = Rectangle(2, 3)
print(resultList["Area"])
print(resultList["Perimeter"])


#====================================
# A simple R program to
# demonstrate the inline function


f = function(x) x^2*4+x/3

print(f(4))
print(f(-2))
print(0)


#====================================


# A simple R program to demonstrate
# Lazy evaluations of functions

Cylinder = function(diameter, length, radius ){
volume = pi*diameter^2*length/4
return(volume)
}

# This'll execute because this
# radius is not used in the
# calculations inside the function.
print(Cylinder(5, 10))

#no es necesario escribir el radio porque no se utilza en los calculos

#============================================

# Function definition
# To check n is divisible by 5 or not
divisbleBy5 <- function(x){
if(x %% 5 == 0)
{
	return("number is divisible by 5")
}
else
{
	return("number is not divisible by 5")
}
}

# Function call
divisbleBy5(100)
divisbleBy5(4)
divisbleBy5(20.0)

#============================================


# Function definition
# To check a is divisible by b or not
divisible <- function(a, b){
if(a %% b == 0)
{
	return(paste(a, "is divisible by", b))
}
else
{
	return(paste(a, "is not divisible by", b))
}
}

# Function call
divisible(7, 3)
divisible(36, 6)
divisible(9, 2)

#============================================

# Function definition of dots operator
fun <- function(n, ...){
l <- list(n, ...)
paste(l, collapse = " ")
}

# Function call
fun(5, 1L, 6i, TRUE, "GeeksForGeeks", "Dots operator")


