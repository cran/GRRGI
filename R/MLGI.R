`MLGI` <-
function(data,prob,N) { #***********************MLGI*************************************************************************
#
#
#          use data and prob to do all lmer (max likelihood) estimates and Generalized Inference Confidence limits
#            inputs are 'data', 'prob' level and 'N' the number of simulations
#            'data' needs columns named 'part', 'operator', and 'resp' (for response)
#            Analysis CROSSED, did not do nested
#            output will be the matrix out.gi.i
#            example: result<-MLGI(data,0.95,10000)
#
#
#
suppressPackageStartupMessages(require(nlme))   # for lmer groupeddata, need package installed
suppressPackageStartupMessages(require(lme4))   # for VarCorr, need package installed
suppressPackageStartupMessages(require(arm))    # for sigma.hat, need to have car package installed for arm
#
#
#
#
#--------------- Create interaction variable ----------------------
#
partOp<-data$part*100+data$operator          # part in 100s place, operator in 1s place
data<-cbind(data,partOp)                     # add to data
#data                                        # look at data (commented out)
#
#--------------- make factor variables ----------------------------
#
resp<-data$resp
part<-as.factor(data$part)
operator<-as.factor(data$operator)
partop<-as.factor(data$partOp)
#
#--------------- get levels ---------------------------------------
#
p<-nlevels(part)                            # levels for part
o<-nlevels(operator)                        # levels for operator
n<-nrow(data)                                # total number of data points
m<-n/nlevels(partop)                        # repeats, m = total / unique points
#
#
#--------------- check if balanced design -------------------------
#  if the design is unbalanced m will not be an integer
#  p*o*m = n   but p*o*round(n) = n only if n is an integer (and the design is balanced)
#
if (p*o*round(m) != n) { print('NOT a Balanced Design')}
#
#
#-------------- conf int level ------------------------------------
#
alpha<- 1-prob
#
#
#------------------------------Test for significance of interaction -----------------------
#
#   do two lmer models, one with interaction and one without interaction (nested models)
#
lmer.i<-lmer(resp~1+(1|part) + (1|operator) + (1|partOp),data)     # with interaction
lmer.ni<-lmer(resp~1+(1|part) + (1|operator),data)                 # without interaction
#
#   do anova on the nested models and extract Pr[2], which is the
#   likelihood ratio for significance of the interaction model over the reduced model
#
Pint.lmer<-anova(lmer.i,lmer.ni)$Pr[2]
#
#
#
#---------------------------------------------------------------------------------------------------
#                             if significant go on here and do with interaction                     |
#                             if not significant skip ahead to do no interaction                    |
#                             a level of 0.25 is used as being significant                          |
#---------------------------------------------------------------------------------------------------
#                                                                                                  
if(Pint.lmer<=0.25) {           # if significant go on, if not skip this
#
#
#
#
# --------------define matrix for storing results----------------------------------------
#
#
out.gi.i<-matrix(NA,17,4)
#
rownames(out.gi.i)<-c("SD:Part","SD:Operator","SD:PartOp","SD:Repeat","SD:Reproduce","SD:Gauge","SD:Total",
   "Var:Part","Var:Operator","Var:PartOp","Var:Repeat","Var:Reproduce","Var:Gauge","Var:Total",
   "Gauge/Tot","Gauge/Parts","Repeat/Gauge")
colnames(out.gi.i)<-c("ML","Lower GI","Mid GI","Upper GI")
#
#
#
#------------------lmer with interaction------------------------------------------------
#
lmer.i<-lmer(resp~1+(1|part) + (1|operator) + (1|partOp),data)
#
#------------------take estimates using VarCorr and sigma.hat --------------------------
#
out.gi.i[8,1]<-VarCorr(lmer.i)$part[1]            # lmer part var
out.gi.i[9,1]<-VarCorr(lmer.i)$operator[1]        # lmer operator var
out.gi.i[10,1]<-VarCorr(lmer.i)$partOp[1]         # lmer PartOp var
#                                                   lmer error will come from sigma.hat below
#
out.gi.i[4,1]<-sigma.hat(lmer.i)$sigma$data       # lmer error stdev - needs 'arm', which needs car installed
#
out.gi.i[11,1]<-out.gi.i[4,1]^2                   # lmer error var from the stdev
#
#-----------------calc other estimates --------------------------------
#
out.gi.i[12,1]<-out.gi.i[9,1]+out.gi.i[10,1]       # lmer repro = Op + PartOp
out.gi.i[13,1]<-out.gi.i[12,1]+out.gi.i[11,1]      # lmer Gauge = repro + error
out.gi.i[14,1]<-out.gi.i[13,1]+out.gi.i[8,1]       # lmer total = Gauge + part
#
out.gi.i[15,1]<-out.gi.i[13,1]/out.gi.i[14,1]      # lmer gauge / total
out.gi.i[16,1]<-out.gi.i[13,1]/out.gi.i[8,1]       # lmer gauge / part
out.gi.i[17,1]<-out.gi.i[11,1]/out.gi.i[13,1]      # lmer repeat / gauge
#
#
#-------------------------------------------------------------------------------------------------------------
#                     Generalized Inference w/interaction
#-------------------------------------------------------------------------------------------------------------
#
#
#---------------need SS and df from ANOVA ----------------------------------
#
lm.anova.i<-anova(lm(resp~part*operator))       # anova of linear regression w/interaction
#
dfi.p<-lm.anova.i$Df[1]                         # get degrees of freedom
dfi.o<-lm.anova.i$Df[2]
dfi.po<-lm.anova.i$Df[3]
dfi.e<-lm.anova.i$Df[4]
#
ssi.p<-lm.anova.i$Sum[1]                        # get sum of squares
ssi.o<-lm.anova.i$Sum[2]
ssi.po<-lm.anova.i$Sum[3]
ssi.e<-lm.anova.i$Sum[4]
#
#
#--------------------define matrices for storing simulated values------------
#
Ti<-matrix(0,N,10)
colnames(Ti)<-c("Part","Op","PartOp","Rep","Repro","Gauge","Total",
                "gauge/total","gauge/part","repeat/gauge")
#
#
#---------------------------------------------------------------------------
#           calculate simulated distribution Ti (with interaction)
#---------------------------------------------------------------------------
#
# set Ci coefficients so the 1st seven columnx of Ti will be:
#
# 1st column will be "part"
# 2nd column will be "operator"
# 3rd column will be "partOp"  interaction
# 4th column will be "repeat"
# 5th column will be "reproducibility"
# 6th column will be "total Gauge or test variability"
# 7th column will be "total variability"
# see ratio columns 8,9,10 below
# 
Ci0<-c(1,0,0,0,0,0,1)
Ci1<-c(0,1,0,0,1,1,1)
Ci2<-c(0,0,1,0,1,1,1)
Ci3<-c(0,0,0,1,0,1,1)
#
#
for(j in 1:7) {                            # one loop for each column of C
   for(i in 1:N) {                            # one loop for each simulation
      #
      Y0<-rchisq(1,dfi.p)                       # generate random variables
      Y1<-rchisq(1,dfi.o)
      Y2<-rchisq(1,dfi.po)
      Y3<-rchisq(1,dfi.e)
      #
      #       calculate the random variable for the estimate (Eq 6 in paper)
      Ti[i,j]<-
        (Ci0[j]*ssi.p)/(m*o*Y0) + 
        (Ci1[j]*ssi.o)/(m*p*Y1) + 
        (1/m)*(Ci2[j]-Ci0[j]/o-Ci1[j]/p)*(ssi.po/Y2) + 
        (Ci3[j]-Ci2[j]/m)*(ssi.e/Y3)
      #
   }
}
#
#------------------ ratios ------------------------------------------------
#
# TOP of ratio
# 1st top column will be "gauge"
# 2nd top column will be "gauge"
# 3rd top column will be "repeat"
# (other ratios can be chosen by adjusting a's and b's)
# (a specific element of a0 to a3 forms an "a" vector as in the paper)
#
ai0<-c(0,0,0)
ai1<-c(1,1,0)
ai2<-c(1,1,0)
ai3<-c(1,1,1)
#
# BOTTOM of ratio
# 1st bottom column will be "total"
# 2nd bottom column will be "Part or process"
# 3rd bottom column will be "gauge"
#
bi0<-c(1,1,0)
bi1<-c(1,0,1)
bi2<-c(1,0,1)
bi3<-c(1,0,1)
#
for(j in 1:3) {                            # loop for the 3 different ratios
   for(i in 1:N) {                            # loop for the N simulated values
      #
      Y0t<-rchisq(1,dfi.p)                       # random vars for top of ratio
      Y1t<-rchisq(1,dfi.o)
      Y2t<-rchisq(1,dfi.po)
      Y3t<-rchisq(1,dfi.e)
      #
      Y0b<-rchisq(1,dfi.p)                       # random vars for bottom of ratio
      Y1b<-rchisq(1,dfi.o)
      Y2b<-rchisq(1,dfi.po)
      Y3b<-rchisq(1,dfi.e)
      #
      #       calculate the random variable for the estimate (Eq 7 in paper) 
      numerator<-
        (ai0[j]*ssi.p)/(m*o*Y0t) + 
        (ai1[j]*ssi.o)/(m*p*Y1t) + 
        (1/m)*(ai2[j]-ai0[j]/o-ai1[j]/p)*(ssi.po/Y2t) + 
        (ai3[j]-ai2[j]/m)*(ssi.e/Y3t)
      #
      denominator<-
        (bi0[j]*ssi.p)/(m*o*Y0b) + 
        (bi1[j]*ssi.o)/(m*p*Y1b) + 
        (1/m)*(bi2[j]-bi0[j]/o-bi1[j]/p)*(ssi.po/Y2b) + 
        (bi3[j]-bi2[j]/m)*(ssi.e/Y3b)
      #
      Ti[i,j+7]<-numerator/denominator       # calculate the ratio simulation, store in columns of Ti after 1st 7
   }
}
#
#-------------------sort the values --------------------------------------------------------------------
#
for(j in 1:10) { Ti[,j]<-sort(Ti[,j]) }                # sort column by column
#
#
#---------------find the first positive value and save that index - ip[j]-------------------------------
#
ip<-c(1,1,1,1,1,1,1,1,1,1)                              # ip will be the vector of indices of 1st positive value
for (j in 1:10) {
    while( Ti[ ip[j] ,j] < 0 ) { ip[j]<-ip[j]+1 }       # increase index as long as the value is negative
}
#
#
#-------------------pick out the confidence limits (and midpoint) ----------------------------------------
#    (alpha/2) quantile to (1-alpha/2) quantile of all positive random estimates
#
#  Ti columns are
#      1="Part", 2="Op", 3="PartOp", 4="Rep", 5="Repro", 6="Gauge", 7="Total", 8="gauge/total", 9="gauge/Part", 10="Repeat/Gauge"
#
#
#  out.gi.i rows are (these are variances, the stddev's for the non-ratios will be added to rows 1 to 7):  
#      8=Part, 9=Oper, 10=PartOp, 11=repeat, 12=Repro, 13=Gauge, 14=Total, 15=Gauge/Tot, 16=Gauge/part, 17=Repeat/Gauge
#
for(j in 1:10) {
   #
   #         index the (1st pos) + (alpha/2 quantile of all pos #'s)
   #
   out.gi.i[j+7,2]<-Ti[ (ip[j]-1) + round( (alpha/2)*(N-(ip[j]-1)) ),j]              # lower Conf int in 2nd col
   out.gi.i[j+7,3]<-Ti[ (ip[j]-1) + round( 0.5*(N-(ip[j]-1)) ),j]                    # Mid point in 3rd col
   out.gi.i[j+7,4]<-Ti[ (ip[j]-1) + round( (1-alpha/2)*(N-(ip[j]-1)) ),j]            # Upper Conf int in 4th col
}
#
# -----------------calc the StDev's --------------------------------------------
#
for(i in 1:7) { for(j in 1:4) { out.gi.i[i,j]<- suppressWarnings( sqrt(out.gi.i[i+7,j]) )  }}
#
#
#------------------------------------------------------------
#        histogram of results - commented out
#
#  index  1=Part, 2=Oper, 3=PartOp, 4=repeat, 5=Repro, 6=gauge, 7=Total, 8=gauge/Tot, 9=Rep/gauge, 10=Repro/gauge
#
#i<-1
#hist(Ti[,i],breaks=100,xlim=c(Ti[round((alpha/2)*N),i],Ti[round((1-alpha/2)*N),i]))
#
#-------------------------------------------------------------
#
#print(out.gi.i, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
#
#
#-----------------------------------------------------------------------------------------------------------------
}            # end "if" of with Interaction
#-----------------------------------------------------------------------------------------------------------------
#
#
#
#
#
#                     now for when the interaction term is NOT significant
#-----------------------------------------------------------------------------------------------------------------
#
#
if(Pint.lmer>0.25) {                    # interaction is not significant 
#
#
#
# --------------define matrix for storing output results------------------------
#
out.gi.ni<-matrix(NA,17,4)
rownames(out.gi.ni)<-c("SD:Part","SD:Operator","SD:PartOp","SD:Repeat","SD:Reproduce","SD:Gauge","SD:Total",
   "Var:Part","Var:Operator","Var:PartOp","Var:Repeat","Var:Reproduce","Var:Gauge","Var:Total",
   "Gauge/Total","Gauge/Part","Repeat/Gauge")
colnames(out.gi.ni)<-c("ML","Lower GI","Mid GI","Upper GI")
#
#------------------lmer without interaction---------------------------------
#
lmer.ni<-lmer(resp~1+(1|part) + (1|operator),data)
#
#------------------take estimates--------------------------------------
#
out.gi.ni[8,1]<-VarCorr(lmer.ni)$part[1]            # lmer part var
out.gi.ni[9,1]<-VarCorr(lmer.ni)$operator[1]        # lmer operator var
#                                                     lmer error will come from sigma.hat below
#
out.gi.ni[4,1]<-sigma.hat(lmer.ni)$sigma$data       # lmer error sd - needs require(arm), which needs car installed
#
out.gi.ni[11,1]<-out.gi.ni[4,1]^2                   # lmer error var
#
#-----------------calc other estimates --------------------------------
#
out.gi.ni[12,1]<-out.gi.ni[9,1]                       # lmer repro = Oper
out.gi.ni[13,1]<-out.gi.ni[12,1]+out.gi.ni[11,1]      # lmer gauge = repro + error
out.gi.ni[14,1]<-out.gi.ni[13,1]+out.gi.ni[8,1]       # lmer total = gauge + part
#
out.gi.ni[15,1]<-out.gi.ni[13,1]/out.gi.ni[14,1]      # lmer gauge / total
out.gi.ni[16,1]<-out.gi.ni[13,1]/out.gi.ni[8,1]       # lmer gauge / part
out.gi.ni[17,1]<-out.gi.ni[11,1]/out.gi.ni[13,1]      # lmer repeat / gauge
#
#
#
#----------------------------------------------------------------------------
#                     Generalized Inference w/NO interaction
#----------------------------------------------------------------------------
#
#
#---------------need SS and df from anova for GI---------------------------
#
lm.anova.ni<-anova(lm(resp~part+operator))     # anova of linear regression w/no int
#
df.p<-lm.anova.ni$Df[1]                         # get degrees of freedom
df.o<-lm.anova.ni$Df[2]
df.e<-lm.anova.ni$Df[3]
#
ss.p<-lm.anova.ni$Sum[1]                        # get sum of squares
ss.o<-lm.anova.ni$Sum[2]
ss.e<-lm.anova.ni$Sum[3]
#
#
#N<-1000            # number of simulations for GI or get from function call
#
#-------------------------------------------------------------------------------
#       calculate Tni (with NO interaction) simulated distribution 
#-------------------------------------------------------------------------------
#
#
#         define matrices for storing simulated values
#
#                                     no interaction (Repro will be Op)
Tni<-matrix(0,N,8)
colnames(Tni)<-c("Part","Op","Rep","gauge","Total",
                 "gauge/total","gauge/part","repeat/gauge")
#
#
#----------set C coefficients so Ti will be matrix with columns:-----------------
#
# 1st column will be "part"
# 2nd column will be "operator" same as reproducibility (no interaction)
# 3rd column will be "repeat"
# 4th column will be "total gauge or test variability"
# 5th column will be "total variability"
#
C0<-c(1,0,0,0,1)
C1<-c(0,1,0,1,1)
C3<-c(0,0,1,1,1)
#
#
for(j in 1:5) {                            # one loop for each estimate
   for(i in 1:N) {                            # one loop for each simulation
   #
   Y0<-rchisq(1,df.p)
   Y1<-rchisq(1,df.o)
   Y3<-rchisq(1,df.e)
   #
   #   This equation similar to Eq 6 of paper but derived for No interaction
   #
   Tni[i,j]<-(C0[j]*ss.p)/(m*o*Y0) + (C1[j]*ss.o)/(m*p*Y1) + (C3[j]-C0[j]/(o*m)-C1[j]/(p*m))*(ss.e/Y3)
   #
   }
}
#
#---------------------------------- for ratios ------------------------------------------------
#
# 1st top column will be "gauge"
# 2nd top column will be "gauge"
# 3rd top column will be "repeat"
#
a0<-c(0,0,0)
a1<-c(1,1,0)
a3<-c(1,1,1)
#
# 1st bottom column will be "total"
# 2nd bottom column will be "part"
# 3rd bottom column will be "gauge"
#
b0<-c(1,1,0)
b1<-c(1,0,1)
b3<-c(1,0,1)
#
for(j in 1:3) {                           # for each of the three ratios
   for(i in 1:N) {                           # N simulations
   #
   Y0t<-rchisq(1,df.p)                       # random vars for top of ratio
   Y1t<-rchisq(1,df.o)
   Y3t<-rchisq(1,df.e)
   #
   Y0b<-rchisq(1,df.p)                       # random vars for bottom of ratio
   Y1b<-rchisq(1,df.o)
   Y3b<-rchisq(1,df.e)
   #
   numerator<-(a0[j]*ss.p)/(m*o*Y0t) + (a1[j]*ss.o)/(m*p*Y1t) + (a3[j]-a0[j]/(o*m)-a1[j]/(p*m))*(ss.e/Y3t)
   denominator<-(b0[j]*ss.p)/(m*o*Y0b) + (b1[j]*ss.o)/(m*p*Y1b) + (b3[j]-b0[j]/(o*m)-b1[j]/(p*m))*(ss.e/Y3b)
   #
   Tni[i,j+5]<-numerator/denominator
   }
}
#
#-------------------sort the values ----------------------------------------
#
for(j in 1:8) { Tni[,j]=sort(Tni[,j]) }                # sort col by col
#
#
#---------------find the first positive value and save that index - ip[j]---------
#
ip<-c(1,1,1,1,1,1,1,1)
#
for (j in 1:8) {
   while( Tni[ ip[j] ,j] < 0 ) { ip[j]<-ip[j]+1 }
}
#
#
#-------------------pick out the confidence limits (and midpoint) ------------
#
# j index for Tni 1=Part, 2=Oper,             3=repeat,            4=gauge,  5=Total,  6=gauge/Tot,  7=Gauge/Part,   8=Repeat/Gauge
# k index for out 8=Part, 9=Oper, 10=PartOp, 11=repeat, 12=Repro, 13=gauge, 14=Total, 15=Gauge/Tot, 16=Gauge/Part,  17=Repeat/Gauge
#
for(j in 1:8) {
   #
   if(j<3) k<-j+7               # use if statements to put results in correct rows
   if(j==3) k<-j+8
   if(j>3) k<-j+9
   #
   out.gi.ni[k,2]<-Tni[ (ip[j]-1) + round( (alpha/2)*(N-(ip[j]-1)) ),j]              # lower in 2nd col
   out.gi.ni[k,3]<-Tni[ (ip[j]-1) + round( 0.5*(N-(ip[j]-1)) ),j]                    # Middle in 3rd col
   out.gi.ni[k,4]<-Tni[ (ip[j]-1) + round( (1-alpha/2)*(N-(ip[j]-1)) ),j]            # Upper in 4th col
}
out.gi.ni[12,]<-out.gi.ni[9,]                        # repro is same as oper
#
#
# -----------------calc the StDev's -------------------------------------------------------------
#
for(i in 1:7) { for(j in 1:4) { out.gi.ni[i,j]<-suppressWarnings( sqrt(out.gi.ni[i+7,j]) ) }}
#
#
#-----------------------------Print results------------------------------------------------------
#
#print(out.gi.ni, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
#
out.gi.i<-out.gi.ni    # use one print function at end so it outputs
#
#
}  # end of 'if' with NO Interaction   -----------------------------------------------------------------------------
#
#
#   Print out the final results
#
print(out.gi.i, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
#
}  # end of function MLGI     ***********************************************************

