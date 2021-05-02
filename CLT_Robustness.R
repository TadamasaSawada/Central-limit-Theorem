#if (!require("gnorm")) install.packages("gnorm")
library(gnorm)
#if (!require("sn")) install.packages("sn")
library(sn)
#if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
#if (!require("scales")) install.packages("scales")
library(scales)

##########################################
## Simulation 1: Descriptive Statistics ##
##########################################

nSessions = 10^6
nSamples = c(4,8,16,32)
lns = length(nSamples)

delta = 4/sqrt(1+4^2)

true_md_skew = qsn(0.5,alpha=4)
true_mn_skew = delta*sqrt(2/pi)

true_sd_skew = sqrt( 1-2*delta^2/pi )
true_sd_pow1 = sqrt( 2*gamma(3/1)/gamma(1/1) )
true_sd_pow4 = sqrt( 2*gamma(3/4)/gamma(1/4) )

mn_norm = rep(NaN,lns*nSessions); mn_skew = rep(NaN,lns*nSessions); mn_pow1 = rep(NaN,lns*nSessions); mn_pow4 = rep(NaN,lns*nSessions)
sd_norm = rep(NaN,lns*nSessions); sd_skew = rep(NaN,lns*nSessions); sd_pow1 = rep(NaN,lns*nSessions); sd_pow4 = rep(NaN,lns*nSessions)

for(i in 1:lns) for(s in 1:nSessions)
{
  j=i-1
  
  sampledValues = rnorm(nSamples[i])
  mn_norm[j*nSessions+s] =   mean(sampledValues)
  sd_norm[j*nSessions+s] =     sd(sampledValues)
  
  sampledValues = rsn(nSamples[i], alpha=4)
  mn_skew[j*nSessions+s] =   mean(sampledValues)
  sd_skew[j*nSessions+s] =     sd(sampledValues)
  
  sampledValues = rgnorm(nSamples[i], alpha=sqrt(2), beta=1)
  mn_pow1[j*nSessions+s] =   mean(sampledValues)
  sd_pow1[j*nSessions+s] =     sd(sampledValues)
  
  sampledValues = rgnorm(nSamples[i], alpha=sqrt(2), beta=4)
  mn_pow4[j*nSessions+s] =   mean(sampledValues)
  sd_pow4[j*nSessions+s] =     sd(sampledValues)
}

colSampleSize = c()
for(i in 1:lns)
  colSampleSize = append(colSampleSize, rep(nSamples[i],nSessions))


## Normalization ##
es_norm =  mn_norm / sd_norm
es_skew = (mn_skew-true_mn_skew) / sd_skew
es_pow1 =  mn_pow1 / sd_pow1
es_pow4 =  mn_pow4 / sd_pow4

mn_skew = (mn_skew-true_mn_skew) / true_sd_skew
mn_pow1 =  mn_pow1 / true_sd_pow1
mn_pow4 =  mn_pow4 / true_sd_pow4

sd_skew = sd_skew / true_sd_skew
sd_pow1 = sd_pow1 / true_sd_pow1
sd_pow4 = sd_pow4 / true_sd_pow4

dDist = data.frame(Mn = c(mn_norm, mn_skew, mn_pow1, mn_pow4),
                   Sd = c(sd_norm, sd_skew, sd_pow1, sd_pow4),
                   Es = c(es_norm, es_skew, es_pow1, es_pow4),
                   SampleSize=rep(colSampleSize,4),
                   Distribution=c( rep("Norm",lns*nSessions),
                                   rep("Skew",lns*nSessions),
                                   rep("Pow1",lns*nSessions),
                                   rep("Pow4",lns*nSessions))
                  )
dDist$SampleSize = as.factor(dDist$SampleSize)

## Quartiles ##
q1_mn_norm = rep(NaN,lns); q1_mn_skew = rep(NaN,lns); q1_mn_pow1 = rep(NaN,lns); q1_mn_pow4 = rep(NaN,lns);
q2_mn_norm = rep(NaN,lns); q2_mn_skew = rep(NaN,lns); q2_mn_pow1 = rep(NaN,lns); q2_mn_pow4 = rep(NaN,lns);
q3_mn_norm = rep(NaN,lns); q3_mn_skew = rep(NaN,lns); q3_mn_pow1 = rep(NaN,lns); q3_mn_pow4 = rep(NaN,lns);

q1_sd_norm = rep(NaN,lns); q1_sd_skew = rep(NaN,lns); q1_sd_pow1 = rep(NaN,lns); q1_sd_pow4 = rep(NaN,lns);
q2_sd_norm = rep(NaN,lns); q2_sd_skew = rep(NaN,lns); q2_sd_pow1 = rep(NaN,lns); q2_sd_pow4 = rep(NaN,lns);
q3_sd_norm = rep(NaN,lns); q3_sd_skew = rep(NaN,lns); q3_sd_pow1 = rep(NaN,lns); q3_sd_pow4 = rep(NaN,lns);

q1_es_norm = rep(NaN,lns); q1_es_skew = rep(NaN,lns); q1_es_pow1 = rep(NaN,lns); q1_es_pow4 = rep(NaN,lns);
q2_es_norm = rep(NaN,lns); q2_es_skew = rep(NaN,lns); q2_es_pow1 = rep(NaN,lns); q2_es_pow4 = rep(NaN,lns);
q3_es_norm = rep(NaN,lns); q3_es_skew = rep(NaN,lns); q3_es_pow1 = rep(NaN,lns); q3_es_pow4 = rep(NaN,lns);

for(i in 1:lns)
{
 # Mean
  qAll = quantile(mn_norm[colSampleSize==nSamples[i]])
  q1_mn_norm[i] = qAll[2]
  q2_mn_norm[i] = qAll[3]
  q3_mn_norm[i] = qAll[4]
  
  qAll = quantile(mn_skew[colSampleSize==nSamples[i]])
  q1_mn_skew[i] = qAll[2]
  q2_mn_skew[i] = qAll[3]
  q3_mn_skew[i] = qAll[4]
  
  qAll = quantile(mn_pow1[colSampleSize==nSamples[i]])
  q1_mn_pow1[i] = qAll[2]
  q2_mn_pow1[i] = qAll[3]
  q3_mn_pow1[i] = qAll[4]
  
  qAll = quantile(mn_pow4[colSampleSize==nSamples[i]])
  q1_mn_pow4[i] = qAll[2]
  q2_mn_pow4[i] = qAll[3]
  q3_mn_pow4[i] = qAll[4]
  
 # Standard Deviation
  qAll = quantile(sd_norm[colSampleSize==nSamples[i]])
  q1_sd_norm[i] = qAll[2]
  q2_sd_norm[i] = qAll[3]
  q3_sd_norm[i] = qAll[4]
  
  qAll = quantile(sd_skew[colSampleSize==nSamples[i]])
  q1_sd_skew[i] = qAll[2]
  q2_sd_skew[i] = qAll[3]
  q3_sd_skew[i] = qAll[4]
  
  qAll = quantile(sd_pow1[colSampleSize==nSamples[i]])
  q1_sd_pow1[i] = qAll[2]
  q2_sd_pow1[i] = qAll[3]
  q3_sd_pow1[i] = qAll[4]
  
  qAll = quantile(sd_pow4[colSampleSize==nSamples[i]])
  q1_sd_pow4[i] = qAll[2]
  q2_sd_pow4[i] = qAll[3]
  q3_sd_pow4[i] = qAll[4]
  
 # Effect Size
  qAll = quantile(es_norm[colSampleSize==nSamples[i]])
  q1_es_norm[i] = qAll[2]
  q2_es_norm[i] = qAll[3]
  q3_es_norm[i] = qAll[4]
  
  qAll = quantile(es_skew[colSampleSize==nSamples[i]])
  q1_es_skew[i] = qAll[2]
  q2_es_skew[i] = qAll[3]
  q3_es_skew[i] = qAll[4]
  
  qAll = quantile(es_pow1[colSampleSize==nSamples[i]])
  q1_es_pow1[i] = qAll[2]
  q2_es_pow1[i] = qAll[3]
  q3_es_pow1[i] = qAll[4]
  
  qAll = quantile(es_pow4[colSampleSize==nSamples[i]])
  q1_es_pow4[i] = qAll[2]
  q2_es_pow4[i] = qAll[3]
  q3_es_pow4[i] = qAll[4]
}

dQuart = data.frame(MnQ1 = c(q1_mn_norm, q1_mn_skew, q1_mn_pow1, q1_mn_pow4),
                    MnQ2 = c(q2_mn_norm, q2_mn_skew, q2_mn_pow1, q2_mn_pow4),
                    MnQ3 = c(q3_mn_norm, q3_mn_skew, q3_mn_pow1, q3_mn_pow4),
                     SdQ1 = c(q1_sd_norm, q1_sd_skew, q1_sd_pow1, q1_sd_pow4),
                     SdQ2 = c(q2_sd_norm, q2_sd_skew, q2_sd_pow1, q2_sd_pow4),
                     SdQ3 = c(q3_sd_norm, q3_sd_skew, q3_sd_pow1, q3_sd_pow4),
                      EsQ1 = c(q1_es_norm, q1_es_skew, q1_es_pow1, q1_es_pow4),
                      EsQ2 = c(q2_es_norm, q2_es_skew, q2_es_pow1, q2_es_pow4),
                      EsQ3 = c(q3_es_norm, q3_es_skew, q3_es_pow1, q3_es_pow4),
                       SampleSize=rep(nSamples,4),
                    Distribution=c( rep("Norm",lns),
                                    rep("Skew",lns),
                                    rep("Pow1",lns),
                                    rep("Pow4",lns))
                   )
dQuart$SampleSize = as.factor(dQuart$SampleSize)

## Plotting Graphs ##
g = ggplot(data = dDist, mapping = aes(x=SampleSize, y=Mn, color=Distribution))
g = g + geom_violin(width=1, position = position_dodge(width = 0.6))
g = g + scale_color_manual(values = c("black", "blue", "pink", "red"))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=MnQ1), size=0.1, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=MnQ2), size=1.0, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=MnQ3), size=0.1, position = position_dodge(width = 0.6))
#g = g + ylim(-2,2)
g + theme_bw()

g = ggplot(data = dDist, mapping = aes(x=SampleSize, y=Sd, color=Distribution))
g = g + geom_violin(width=1, position = position_dodge(width = 0.6))
g = g + scale_color_manual(values = c("black", "blue", "pink", "red"))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=SdQ1), size=0.1, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=SdQ2), size=1.0, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=SdQ3), size=0.1, position = position_dodge(width = 0.6))
#g = g + ylim(0.1,2)
g + theme_bw() + scale_y_continuous(trans = log2_trans(),limits = c(0.125, 4))

g = ggplot(data = dDist, mapping = aes(x=SampleSize, y=Es, color=Distribution))
g = g + geom_violin(width=1, position = position_dodge(width = 0.6))
g = g + scale_color_manual(values = c("black", "blue", "pink", "red"))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=EsQ1), size=0.5, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=EsQ2), size=1.0, position = position_dodge(width = 0.6))
g = g + geom_point(data=dQuart, mapping = aes(x=SampleSize, y=EsQ3), size=0.5, position = position_dodge(width = 0.6))
g = g + ylim(-2,2)
g + theme_bw()
