# Modeling migratory timing for Upper Cook Inlet sockeye run using R 
# model: y = 1/(1+exp(-a - b*d))
# Mundy (1979) and ADFG report (e.g. Aaron Dupuis and Mark Willette 2116 Report)
#It is the same as old FORTRAN program:
#O F F S H O R E   T E S T   F I S H I N G
#N 0 N - L I N E A R   C U R V E   F I T   P R O G R A M 

obs.cpue<-read.table("obsCPUE2022.csv", sep = ",", header=T) #input observed daily cpue from 2022 OTF test fishery table
x<-obs.cpue$dailyCPUE #to calculate cumulated cpue below
obs.cpue$ccumCPUE[!is.na(x)] <- cumsum(na.omit(x)) #calculate accumulated cpue from daily cpue
#n<-length(obs.cpue$d) 
n<-length(obs.cpue$d) - 1 #no data on the last day (7/31); remove it.
ccpuef.D<-with(obs.cpue, ccumCPUE[n]) #ccpuef on day D 
obs.cpue$yd<-with(obs.cpue, ccumCPUE/ccpuef.D) #obs.proportion of ccpue by day d. It is equation (7).

#Method 1
a_start<- 0 
b_start<- 0
#the formula for the model:yd ~ 1/(1 + exp(-a - b*d
y_formula <- formula(yd ~ 1/(1 + exp(-a - b*d)))
#fit the model
m<-nls(y_formula,start=list(a=a_start,b=b_start), data=obs.cpue)
#estimated parameters
summary(m)
coeffs <- coefficients(m)
coeffs

#The midpoint of the run, defined as the day that approximately 50% of 
#the total run has passed the southern OTF transect from day 1 (June 24)
midrun<- -coeffs[1]/coeffs[2] # -a/b
names(midrun) <- NULL
midrun
#the date of 50% of the total run have passed
as.Date("06/24", '%m/%d') + round(midrun) -1

#Method 2:
#run_timing model function
#model: y = 1/(1+exp(-a + b*d))
mydata<-obs.cpue
run_timing<- function (mydata) 
{
  d <- mydata$d
  y <- mydata$yd
  a <- -1  #initial vaule for running 2021
  b <- 0.1 #initial vaule
  input <- c(a, b) #model parameters to be estimated by minimizing SSE below in ssefn function
  ssefn <- function(input) {
    a <- input[1]
    b <- input[2]
    y.hat <- NULL 
    days <- 1:length(d) 
    for (i in days) {
      y.hat[i] <- 1/(1 + exp(-a - b*d[i]))
    }
    SSE <- sum((y - y.hat)^2, na.rm = TRUE)
    return(SSE) 
  }
  estA <- nlm(ssefn, input, typsize = input, iterlim = 1000)
  estA <- nlm(ssefn, estA$est, typsize = input, iterlim = 1000)
  a <- estA$est[1]
  b <- estA$est[2]

  out <- list(a = a, b = b, midrun = -a/b)
}

s<-run_timing(mydata) #call run_timing model function
s
#the date of 50% of the total run have passed
as.Date("06/24", '%m/%d') + round(s$midrun) -1
