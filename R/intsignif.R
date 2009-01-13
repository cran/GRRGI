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
#--------------- Create interaction variable ----------------------
#
partOp<-data$part*100+data$operator            # part in 100s place, operator in 1s place
data<-cbind(data,partOp)                       # add to data
#
#--------------- make factor variables ----------------------------
resp<-data$resp
part<-as.factor(data$part)
operator<-as.factor(data$operator)
partop<-as.factor(data$partOp)
#
#
#-------------------------------do the ANOVA---------------------
#
lm.anova.i<-anova(lm(resp~part*operator))     # anova of linear regression w/interaction
#
Pint.lm<-lm.anova.i$Pr[3]
#
#
#------------------------------do the lmer-----------------------
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
} # end of  intsignif   ----------------------------------------------------

