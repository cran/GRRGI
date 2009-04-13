`intsignif` <-
function(data) {  #***********************intsignif********************************************************
#
#
#            calculate and report interaction significance for ANOVA and for LMER
#            only "data" needs to be passed
#
#
#
require(nlme,warn.conflicts = FALSE)           # for lmer groupeddata, need package installed
require(lme4,warn.conflicts = FALSE)           #
#
#--Create an interaction variable---------------------------------------------------------------
#  place the part number in the 1000s place and the operator number in the ones place
#
partOp<-data$part*1000+data$operator           # part in 1000s place, operator in 1s place
data<-cbind(data,partOp)                       # add to data
#
#--Make the variables into factor variables-------------------------------------------------------
#
resp<-data$resp
part<-as.factor(data$part)
operator<-as.factor(data$operator)
partop<-as.factor(data$partOp)
#
#
#--Run ANOVA and extract the P value for the interaction factor-----------------------------------
#
lm.anova.i<-anova(lm(resp~part*operator))     # anova of linear regression w/interaction
#
Pint.lm<-lm.anova.i$Pr[3]
#
#
#--Test for significance of the interaction term---------------------------------------------------
#  fit two lmer (linear mixed effects) models, one with and one without the 
#  interaction term. Do ANOVA on the nested models and extract Pr[2], which is the 
#  likelihood ratio for significance of the full model (with interaction) over the 
#  reduced model.
#
lmer.i<-lmer(resp~1+(1|part) + (1|operator) + (1|partOp),data)
#
lmer.ni<-lmer(resp~1+(1|part) + (1|operator),data)
#
Pint.lmer<-anova(lmer.i,lmer.ni)$Pr[2]
#
#
#-----------------results------------------------------------------------
#
result<-matrix(0,2,1)
rownames(result)<-c("ANOVA Interaction P value","Likelihood Ratio Test for Interaction")
colnames(result)<-(" ")
#
result[1,1]<-Pint.lm
result[2,1]<-Pint.lmer
result
} # end of  intsignif   ------------------------------------------------------------------------------------------

