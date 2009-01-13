`RRgraph` <-
function(sat1,GI1) {  #******************RRgraph************************************************
#
#           Graph confidence intervals of anovasat and MGGI outputs
#           need a matrix of results from anovasat and MLGI
#
#
#
#               plot repeat, repro, gauge   with interaction
#
plot(sat1[11,],c(1,1,1),type="b",xlim=c(0,max(GI1[11:13,])),
ylim=c(0,6),xlab="",ylab="        repeat                  reproduce                  Gauge",
main="Confidence limits",yaxt="n")
#
lines(GI1[11,c(1,2,4)],c(2,2,2),col=1)
points(GI1[11,c(1,2,4)],c(2,2,2),col=1)
#
for(i in 1:2){
   lines(sat1[i+11,],c(2*i+1,2*i+1,2*i+1),col=i+1)
   points(sat1[i+11,],c(2*i+1,2*i+1,2*i+1),col=i+1)
   #
   lines(GI1[i+11,c(1,2,4)],c(2*i+2,2*i+2,2*i+2),col=i+1)
   points(GI1[i+11,c(1,2,4)],c(2*i+2,2*i+2,2*i+2),col=i+1)
}
Axis(x =NULL, at = c(1,2,3,4,5,6),side=2, labels = c("Satt","GI","Satt","GI","Satt","GI"))
#
}   # end of RRgraph *****************************************************************

