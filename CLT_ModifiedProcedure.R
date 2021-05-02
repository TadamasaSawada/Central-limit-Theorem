#if (!require("sn")) install.packages("sn")
#library(sn)
#if (!require("tidyverse")) install.packages("tidyverse")
#library(tidyverse)
#if (!require("scales")) install.packages("scales")
#library(scales)
#if (!require("moments")) install.packages("moments")
#library(moments)
#install.packages("dplyr")
#library(dplyr)

###########################
## Skewness and Kurtosis ##
###########################

nSamples0 = 10^4 # num MC-sessions
nSamples1 = 10^3 # num trials for each participant
nSamples2 = c(4,32,256) # num Participants
lns2 = length(nSamples2)

count_psw = rep(0,lns2)
count_pks = rep(0,lns2)
estimated_skw = matrix(rep(0,lns2*nSamples0),lns2,nSamples0)
estimated_krt = matrix(rep(0,lns2*nSamples0),lns2,nSamples0)

for(i2 in 1:lns2)
{
  for(i0 in 1:nSamples0)
  {
    myu = rnorm(nSamples2[i2])
    myu[myu<=0] = -2
    myu[myu> 0] = +2
    myu = myu + rnorm(nSamples2[i2])
    
    sig = rlnorm(nSamples2[i2],meanlog=1,sdlog=0.5)
    skw = runif(nSamples2[i2],min=2,max=4)
    
    averaged_var = rep(0,nSamples1)
  
    for(i3 in 1:nSamples2[i2])
    {
      sampled_var = rsn(nSamples1, xi=myu[i3], omega=sig[i3], alpha=skw[i3], tau=0)
      averaged_var = averaged_var + sampled_var
    }
    averaged_var = averaged_var/nSamples2[i2]
    
    resultKS = ks.test(x=averaged_var,y="pnorm", mean=mean(averaged_var), sd=sd(averaged_var))
    resultSW = shapiro.test(averaged_var)
    
    estimated_skw[i2,i0] = skewness(averaged_var)
    estimated_krt[i2,i0] = kurtosis(averaged_var)
    
    if(resultSW$p.value<0.05)
      count_psw[i2] = count_psw[i2] + 1
    if(resultKS$p.value<0.05)
      count_pks[i2] = count_pks[i2] + 1
  }
}

count_psw = count_psw/nSamples0
count_pks = count_pks/nSamples0

dModified = data.frame(Skew = as.vector(t(estimated_skw)),
                       Kurt = as.vector(t(estimated_krt)),
                       SampleSize=c(rep(nSamples2[1],nSamples0),
                                    rep(nSamples2[2],nSamples0),
                                    rep(nSamples2[3],nSamples0))
                      )
dModified$SampleSize = as.factor(dModified$SampleSize)


skwQ123 = c()
krtQ123 = c()
for(i2 in 1:lns2)
{
  skwQA = quantile(estimated_skw[i2,])
  skwQ123 = c(skwQ123,skwQA[2:4])
  krtQA = quantile(estimated_krt[i2,])
  krtQ123 = c(krtQ123,krtQA[2:4])
}
dModifiedQ = data.frame(SkewQ = skwQ123,
                        KurtQ = krtQ123,
                        Q123 = rep(c("Q1","Q2","Q3"),lns2),
                        SampleSize=c(rep(nSamples2[1],3),
                                     rep(nSamples2[2],3),
                                     rep(nSamples2[3],3))
                       )
dModifiedQ$SampleSize = as.factor(dModifiedQ$SampleSize)


## Plotting Graphs ##
g = ggplot(data = dModified, mapping = aes(x=SampleSize, y=Skew))
g = g + geom_violin(width=1)
g = g + geom_point(data=dModifiedQ, mapping = aes(x=SampleSize, y=SkewQ), size=0.1)
g = g + ylim(-0.5,1)
g + theme_bw()

g = ggplot(data = dModified, mapping = aes(x=SampleSize, y=Kurt))
g = g + geom_violin(width=1)
g = g + geom_point(data=dModifiedQ, mapping = aes(x=SampleSize, y=KurtQ), size=0.1)
g = g + ylim(2.5,4)
g + theme_bw()

