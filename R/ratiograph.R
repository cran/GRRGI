`ratiograph` <-
function(sat1,GI1) {  #*****************ratiograph**********************************************
#
#           Graph confidence intervals of ratios from anovasat and MLGI
#           need a matrix of results from anovasat and MLGI
#
#
#
plot(sat1[15,],c(1,1,1),type="b",xlim=c(0,1),
ylim=c(0,6),xlab="(graph truncated at 1.0)",ylab="         Gauge/Total            Repeat/Gauge           Repro/Gauge",
main="Ratio Confidence limits",yaxt="n")
#
lines(GI1[15,c(1,2,4)],c(2,2,2),col=1)
points(GI1[15,c(1,2,4)],c(2,2,2),col=1)
#
for(i in 1:2){
lines(sat1[i+15,],c(2*i+1,2*i+1,2*i+1),col=i+1)
points(sat1[i+15,],c(2*i+1,2*i+1,2*i+1),col=i+1)
#
lines(GI1[i+15,c(1,2,4)],c(2*i+2,2*i+2,2*i+2),col=i+1)
points(GI1[i+15,c(1,2,4)],c(2*i+2,2*i+2,2*i+2),col=i+1)
}
Axis(x =NULL, at = c(1,2,3,4,5,6),side=2, labels = c("ANOVA","GI","ANOVA","GI","ANOVA","GI"))
#
#
}   # end of ratiograph  ***********************************************************************

