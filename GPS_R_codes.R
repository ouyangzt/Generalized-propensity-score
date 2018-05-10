###################################################################################
# Estimating generalized propensity score#################################################
# Set work directory and read data
setwd("C:/Data")
xstream<-read.table("stream.csv",header=T,sep=",")

#########################################################################################################################Functions#########################################################
##Function theta estimate the dose-response function with covariates
theta<-function (x,xstream) #x is a vector, a number sequence from 1 to the length of xstream
  #xstream is the dataframe with all covariates, the treatment (logTN), and the response (IR) included
{
  stream=xstream[x,]
  #fit a multiple linear model to predict logTN, This is the implementation of equation (2).
  pre_mlr<-lm(logTN~ELEV+longitude+logprecip+logAREA+logCL+logHCO3+logSO4+SED+STRMTEMP+Percent.AGT+Percent.URB+Percent.Canopy+Riparian.Disturb.,data=stream)
  pre_mlr_sum=summary(pre_mlr)
  #get the standard deviation of residuals
  sdd=pre_mlr_sum$sigma
  #get the residuals of the fitted model
  resi=residuals(pre_mlr)
  #estimate the propensity score at each treatment level with its associated covariates. This is the implementation of equation (3)
  ps=dnorm(resi,0,sdd)
  # add ps to the original data matrix and save the matrix
  stream_ps=stream
  stream_ps$ps=ps  
  #write.csv(stream_ps,"stream_ps.csv")
  #get the predicted logTN, which will be used in the estimation of the final dose-function
  pre=predict(pre_mlr)
  
  #define and assign variables for the polynomial model, e.g. equation 4
  ps_sq=ps^2
  logTN=stream$logTN
  logTN_sq=logTN^2
  Taxonrich=stream$Taxonrich
  lmps<-lm(Taxonrich~logTN+logTN_sq+ps+ps_sq+logTN*ps) # this is the implementation of equation (4)
  EY=rep(0,length(x)) # define the vector to put in the expected outcome (IR) at each treatment level
  # This is the implementation of equation (5), the inside loop means, for each site, if observed covariates from another site would have been observed with it, what would be the probability of the treatment associated with it
  for (i in 1:length(x)){for(j in 1:length(x)){EY[i]=EY[i]+lmps$coefficients[1]+lmps$coefficients[2]*logTN[i]+lmps$coefficients[3]*logTN_sq[i]+lmps$coefficients[4]*dnorm(logTN[i]-pre[j],0,sdd)+lmps$coefficients[5]*dnorm(logTN[i]-pre[j],0,sdd)*dnorm(logTN[i]-pre[j],0,sdd)+lmps$coefficients[6]*logTN[i]*dnorm(logTN[i]-pre[j],0,sdd)}}
  EY=EY/length(x)
  return(EY) # return the expected outcome at each treatment level averaged over the score distribution
  
}

#############################################################################


## Estimate the dose-response function and draw a figure to show the result
EY=theta (1:670,xstream)

par(mfrow=c(1,1),mar=c(4.2,4.2,1,0.5),oma=c(0,0,0,0))
plot(xstream$logTN,xstream$Taxonrich,xlim=c(1.0,4.1),
     ylab="Invertabrate Richness",xlab=expression(paste("Log TN (",mu,"g/L)")),cex.lab=1.2)
lines(spline(xstream$logTN,EY),col="blue",lty=1,lwd=2)


