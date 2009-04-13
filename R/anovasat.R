`anovasat` <-
function(data,prob) { #***********************anovasat ***********************************
#
#          use data and prob to do Gauge R&R ANOVA estimates and satterthwait limits
#            'data' needs columns named 'part', 'operator', and 'resp' (for response)
#
#            output will be the matrix 'out.A.i'
#            example: result<-anovasat(data,0.95)
#
#
#
require(nlme)# need package installed
#
#
#
#--Create an interaction variable ---------------------------------------------------------------
#  by placing the part number in the 1000s place and the operator number in the ones place
#
partOp<-data$part*1000+data$operator          # part in 1000s place, operator in 1s place
data<-cbind(data,partOp)                     # add to data
#data                                        # look at data, commented out
#
#
#--Make the variables into factor variables-------------------------------------------------------
#
resp<-data$resp
part<-as.factor(data$part)
operator<-as.factor(data$operator)
partop<-as.factor(data$partOp)
#
#
#--Extract p,o,n and calculate m; the number of levels for each factor. --------------------------
#
p<-nlevels(part)                            # levels for part
o<-nlevels(operator)                        # levels for operator
n<-nrow(data)                               # total number of data points
m<-n/nlevels(partop)                        # repeats, m = total / unique points
#
#
#--check if balanced design ----------------------------------------------------------------------
#  by using the fact that p*o*m = n only if parts and operators are balanced and the 
#  fact that m is only an integer if the repeats are balanced. So p*o*round(m) = n only 
#  if m is an integer and the entire design is balanced. 
#  
#
if (p*o*round(m) != n) { return(print('NOT a Balanced Design'))}
#
#
#--Get the alpha level from the desired probability level ----------------------------------------
#  which was passed in with the function call
#
alpha<- 1-prob
#
#
#--------------------------------------------------------------------------------------------------
#  Do an ANOVA and test for significance of the interaction term. If it is significant go on
#  here and do the analysis with an interaction term, if it is not significant, skip ahead to 
#  do the 'no' interaction analysis. A level of P < 0.25 is used as being significant.
#--------------------------------------------------------------------------------------------------
#
lm.anova.i<-anova(lm(resp~part*operator))     # anova of linear regression w/interaction            
#                                                                                                   
if(lm.anova.i$Pr[3]<=0.25) {  # compare interaction term to 0.05, if signif go on, if not skip this 
#
#
#******************************ANOVA w/Interaction*******************************************************
#
#
#--Define a matrix for storing final ANOVA results--------------------------------------------------
#
out.A.i<-matrix(NA,17,3)
rownames(out.A.i)<-c("SD:Part","SD:Operator","SD:PartOp","SD:Repeat","SD:Reproduce","SD:Gauge","SD:Total",
   "Var:Part","Var:Operator","Var:PartOp","Var:Repeat","Var:Reproduce","Var:Gauge","Var:Total",
   "Gauge/Tot","Gauge/Parts","Repeat/Gauge")
colnames(out.A.i)<-c("ANOVA est","Lower STH","Upper STH")
#
#
#--Extract degrees of freedom, sum of squares and mean squares from ANOVA----------------------------
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
msi.p<-lm.anova.i$Mean[1]                       # get mean squares
msi.o<-lm.anova.i$Mean[2]
msi.po<-lm.anova.i$Mean[3]
msi.e<-lm.anova.i$Mean[4]
#
#
#--Calculate the ANOVA variance estimates-------------------------------------------------------------
#  'part, operator, part*operator, error, reproducibility, gauge, and total study 
#  and put results into the output matrix (replace negatives with zero)
#  and also put into output matrix
#
vari.p<-max((msi.p - msi.po)/(o*m),0)             # value or 0 if neg
out.A.i[8,1]<-vari.p
#
vari.o<-max((msi.o - msi.po)/(p*m),0)             # value or 0 if neg
out.A.i[9,1]<-vari.o
#
vari.po<-max((msi.po - msi.e)/m,0)                # value or 0 if neg
out.A.i[10,1]<-vari.po
#
vari.e<-msi.e
out.A.i[11,1]<-vari.e
#
#
vari.repro<-vari.o+vari.po
out.A.i[12,1]<-vari.repro
#
vari.gauge<-vari.repro+vari.e
out.A.i[13,1]<-vari.gauge
#
vari.total<-vari.gauge+vari.p
out.A.i[14,1]<-vari.total
#
#--calculate all std dev's --------------------------------------------------------------------------
#  Take the square root of all the variances and place them into the first 7 rows of out.A.i.
#
sdi.p<-sqrt(vari.p)
sdi.o<-sqrt(vari.o)
sdi.po<-sqrt(vari.po)
sdi.e<-sqrt(vari.e)
#
sdi.repro<-sqrt(vari.repro)
sdi.gauge<-sqrt(vari.gauge)
sdi.total<-sqrt(vari.total)
#
for(i in 1:7){ out.A.i[i,1]<-sqrt( out.A.i[i+7,1] ) } 
#
#--Calculate the ratios------------------------------------------------------------------------------
#  gauge/total, gauge/part, and error/gauge; place them into the output matrix
#
out.A.i[15,1]<-vari.gauge/vari.total
out.A.i[16,1]<-vari.gauge/vari.p
out.A.i[17,1]<-vari.e/vari.gauge
#
#
#****************************************************************************************************
#             confidence intervals - with interaction
#***************************************************************************************************
#
#
#--Repeatability------------------------------------------------------------------------------------
#
xupper<-qchisq(1-alpha/2,dfi.e)             # upper chisquare alpha/2
xlower<-qchisq(alpha/2,dfi.e)               # lower chisquare alpha/2
#
L95.vari.e<-dfi.e*vari.e/xupper             # lower CL for Var
out.A.i[11,2]<-L95.vari.e
#
L95.sdi.e<-sqrt(L95.vari.e)                     # lower CL for SD
out.A.i[4,2]<-L95.sdi.e
#
U95.vari.e<-dfi.e*vari.e/xlower             # upper CL for var
out.A.i[11,3]<-U95.vari.e
#
U95.sdi.e<-sqrt(U95.vari.e)                     # upper CL for SD
out.A.i[4,3]<-U95.sdi.e
#
#
#
#--Reproducibility-----------------------------------------------------------------------------------
#
v1<- (msi.o/(m*p))^2/(o-1)                      # bottom 1st term
v2<- ( ( (p-1)/(m*p) ) * msi.po )^2 / ( (o-1)*(p-1) )       # 2nd term
v3<- (msi.e/m)^2 / ( o * p * (m-1) )            # 3rd term
v<- (vari.repro)^2 / ( v1 + v2 + v3 )           # estimated df for repro
#
xupper<-qchisq(1-alpha/2,v)                     # upper chisquare alpha/2
xlower<-qchisq(alpha/2,v)                       # lower chisquare alpha/2
#
L95.vari.repro<-v*vari.repro/xupper             # lower CL for Var
out.A.i[12,2]<-L95.vari.repro
#
L95.sdi.repro<-sqrt(L95.vari.repro)             # lower CL for SD
out.A.i[5,2]<-L95.sdi.repro
#
U95.vari.repro<-v*vari.repro/xlower             # upper CL for var
out.A.i[12,3]<-U95.vari.repro
#
U95.sdi.repro<-sqrt(U95.vari.repro)             # upper CL for SD
out.A.i[5,3]<-U95.sdi.repro
#
#
#--Total Gauge--------------------------------------------------------------------------------------
#
u1<- (msi.o/(m*p))^2/(o-1)                      # bottom 1st term
u2<- ( ( (p-1)/(m*p) ) * msi.po )^2 / ( (o-1)*(p-1) )    # 2nd term
u3<- ( (m-1) * msi.e/m)^2 / ( o * p * (m-1) )            # 3rd term
u<- (vari.gauge)^2 / ( v1 + v2 + v3 )           # estimated df for repro
#
xupper<-qchisq(1-alpha/2,u)                     # upper chisquare alpha/2
xlower<-qchisq(alpha/2,u)                       # lower chisquare alpha/2
#
L95.vari.gauge<-u*vari.gauge/xupper             # lower CL for Var
out.A.i[13,2]<-L95.vari.gauge
#
L95.sdi.gauge<-sqrt(L95.vari.gauge)             # lower CL for SD
out.A.i[6,2]<-L95.sdi.gauge
#
U95.vari.gauge<-u*vari.gauge/xlower             # upper CL for var
out.A.i[13,3]<-U95.vari.gauge
#
U95.sdi.gauge<-sqrt(U95.vari.gauge)             # upper CL for SD
out.A.i[6,3]<-U95.sdi.gauge
#
#
#
#print(out.A.i, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
#
} # end of if  *******************************************************************************************
#                                    end of interaction
#                     now for when the interaction term is NOT significant
#*********************************************************************************************************
#
#
#
if(lm.anova.i$Pr[3]>0.25) {                            # if interaction term is not significant
#
#
#**********************************************************************************************************
#                                       NO interaction
#**********************************************************************************************************
#
#--Define a matrix for storing final ANOVA results------------------------------------------------------
#   
#
out.A.ni<-matrix(NA,17,3)
rownames(out.A.ni)<-c("SD:Part","SD:Operator","SD:PartOp","SD:Repeat","SD:Reproduce","SD:Gauge","SD:Total",
   "Var:Part","Var:Operator","Var:PartOp","Var:Repeat","Var:Reproduce","Var:Gauge","Var:Total",
   "Gauge/Total","Gauge/Part","Repeat/Gauge")
colnames(out.A.ni)<-c("ANOVA est","Lower STH","Upper STH")
#
#
#--ANOVA with no interaction-----------------------------------------------------------------------------
#  and Extract degrees of freedom, sum of squares and mean squares from ANOVA
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
ms.p<-lm.anova.ni$Mean[1]                       # get mean squares
ms.o<-lm.anova.ni$Mean[2]
ms.e<-lm.anova.ni$Mean[3]
#
#--Calculate the ANOVA variance estimates-------------------------------------------------------------
#
var.p<-max((ms.p - ms.e)/(o*m),0)                    # value or 0 if neg
out.A.ni[8,1]<-var.p
#
var.o<-max((ms.o - ms.e)/(p*m),0)                    # value or 0 if neg
out.A.ni[9,1]<-var.o
#
var.e<-ms.e
out.A.ni[11,1]<-var.e
#
#
var.repro<-var.o
out.A.ni[12,1]<-var.repro
#
var.gauge<-var.repro+var.e
out.A.ni[13,1]<-var.gauge
#
var.total<-var.gauge+var.p
out.A.ni[14,1]<-var.total
#
#
#--calculate all std dev's --------------------------------------------------------------------------
#
sd.p<-sqrt(var.p)
sd.o<-sqrt(var.o)
sd.e<-sqrt(var.e)
#
sd.repro<-sqrt(var.repro)
sd.gauge<-sqrt(var.gauge)
sd.total<-sqrt(var.total)
#
for(i in 1:7){ out.A.ni[i,1]<-sqrt( out.A.ni[i+7,1] ) }
#
#--Calculate the ratios------------------------------------------------------------------------------
#
out.A.ni[15,1]<-var.gauge/var.total
out.A.ni[16,1]<-var.gauge/var.p
out.A.ni[17,1]<-var.e/var.gauge
#
#
#
#
#
#***************************************************************************
#             confidence intervals - no interaction
#***************************************************************************
#
#
#--Repeatability------------------------------------------------------------------------------------
#
xupper<-qchisq(1-alpha/2,df.e)             # upper chisquare alpha/2
xlower<-qchisq(alpha/2,df.e)               # lower chisquare alpha/2
#
L95.var.e<-df.e*var.e/xupper                   # lower CL for Var
out.A.ni[11,2]<-L95.var.e
#
L95.sd.e<-sqrt(L95.var.e)                     # lower CL for SD
out.A.ni[4,2]<-L95.sd.e
#
U95.var.e<-df.e*var.e/xlower                   # upper CL for var
out.A.ni[11,3]<-U95.var.e
#
U95.sd.e<-sqrt(U95.var.e)                     # upper CL for SD
out.A.ni[4,3]<-U95.sd.e
#
#
#--Reproducibility-----------------------------------------------------------------------------------
#
xupper.o<-qchisq(1-alpha/2,df.o)               # also need chisquare with operator df
xupper<-qchisq(1-alpha/2,df.e)
xlower.o<-qchisq(alpha/2,df.o)
xlower<-qchisq(alpha/2,df.e)
#
c1<- df.o * ms.o / xupper.o                    # 1st 'c' top term
c2<- df.e * ms.e / xlower                      # 2nd 'c' top term, uses df from error
c<- max( (c1 - c2) / (p*m), 0) 
#
d1<- df.o * ms.o / xlower.o                    # 1st 'd' top term
d2<- df.e * ms.e / xupper                      # 2nd 'd' top term
d<- ( d1-d2 ) / (p*m)
#
#
L95.var.repro<-c                              # lower CL for Var
out.A.ni[12,2]<-L95.var.repro
#
L95.sd.repro<-sqrt(L95.var.repro)             # lower CL for SD
out.A.ni[5,2]<-L95.sd.repro
#
U95.var.repro<-d                              # upper CL for var
out.A.ni[12,3]<-U95.var.repro
#
U95.sd.repro<-sqrt(U95.var.repro)             # upper CL for SD
out.A.ni[5,3]<-U95.sd.repro
#
#
#--Total Gauge--------------------------------------------------------------------------------------
#
w1<- (ms.o/(m*p))^2/(o-1)                       # bottom 1st term
w2<- ( (m*p-1) * ms.e / (m*p) )^2 / ( df.e )    # bottom 2nd term
w<- (var.gauge)^2 / ( w1 + w2 )                 # estimated df for gauge
#
xupper<-qchisq(1-alpha/2,w)                     # upper chisquare alpha/2
xlower<-qchisq(alpha/2,w)                       # lower chisquare alpha/2
#
L95.var.gauge<-w*var.gauge/xupper               # lower CL for Var
out.A.ni[13,2]<-L95.var.gauge
#
L95.sd.gauge<-sqrt(L95.var.gauge)               # lower CL for SD
out.A.ni[6,2]<-L95.sd.gauge
#
U95.var.gauge<-w*var.gauge/xlower               # upper CL for var
out.A.ni[13,3]<-U95.var.gauge
#
U95.sd.gauge<-sqrt(U95.var.gauge)               # upper CL for SD
out.A.ni[6,3]<-U95.sd.gauge
#
#
#
#
#print(out.A.ni, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
out.A.i<-out.A.ni     # also put it into the with/interaction matrix so can print out one for everything at end
#
#
#
} #--------------- end of "if" for no interaction   -------------------------------------------------
#
#
#
#
#-----------------------------  print out -----------------------------------------------------------
#
print(out.A.i, digits = 3, quote = FALSE, na.print = "-", print.gap = 3, right = TRUE)
#
#
#
#
} # end of function anovasat *********************************************************************************************************************************

