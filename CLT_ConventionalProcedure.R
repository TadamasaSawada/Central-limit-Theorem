#if (!require("tidyverse")) install.packages("tidyverse")
#library(tidyverse)
#if (!require("scales")) install.packages("scales")
#library(scales)

################################################
## Conventional method (Bimodal distribution) ##
################################################

nSamples1 = 10^6
nSamples2 = c(4,32)
lns2 = length(nSamples2)

mHist1 = matrix(rep(0,lns2*nSamples1),lns2,nSamples1);

for(i1 in 1:nSamples1)
{
  for(i2 in 1:lns2)
  {
    myu = rnorm(1);
    if(rnorm(1)<0)
      myu = myu-2
    else
      myu = myu+2
    sig = rlnorm(1,meanlog=1,sdlog=0.5)
    
    averaged_var = mean(rnorm(nSamples2[i2], mean=myu, sd=sig))
    
    mHist1[i2,i1] = averaged_var
  }
}

x=seq(-6,+6,0.01)
y=dnorm(x,mean=-2,sd=1)/2+dnorm(x,mean=+2,sd=1)/2

i2=1
hist(mHist1[i2,], xlim=c(-6,+6), ylim=c(0,0.25), breaks=c(-Inf,seq(-6,+6, by=0.1),+Inf))
lines(x, y, lwd=2, col="blue")

i2=2
hist(mHist1[i2,], xlim=c(-6,+6), ylim=c(0,0.25), breaks=c(-Inf,seq(-6,+6, by=0.1),+Inf))
lines(x, y, lwd=2, col="blue")
