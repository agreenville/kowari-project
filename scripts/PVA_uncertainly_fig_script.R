###############################
# coe for fig 2.
# Ellner S. P. & Holmes E. E. (2008)
# Commentary on Holmes et al. (2007): resolving the debate on when extinction risk is predictable.
# Ecol. Lett. 11, E1-E5.

Pext=function(U,V) {pnorm(U-V)+exp(2*U*V+pnorm(-(U+V),log.p=T))}
TMU=function(T,RelNe,side="one") {
  qt=switch(EXPR = side,one=qnorm(0.95),two=qnorm(0.975));
  a=log(1/RelNe);
  U =-mu*sqrt(T/sigma2.b);
  V = a/( sqrt(sigma2.b*T) );
  upci = Pext(U + qt*sqrt(T/N),V);
  lowci = Pext(U - qt*sqrt(T/N),V);
  return(list(upper=upci,lower=lowci))
}
TMUpper=function(T,RelNe,side="one") TMU(T,RelNe,side)$upper
TMLower=function(T,RelNe,side="one") TMU(T,RelNe,side)$lower
CIWidth=function(T,RelNe) {
  qt=qnorm(0.975);
  a=log(1/RelNe);
  U =-mu*sqrt(T/sigma2.b);
  V = a/( sqrt(sigma2.b*T) );
  upci = Pext(U + qt*sqrt(T/N),V);
  lowci = Pext(U - qt*sqrt(T/N),V);
  return(upci-lowci)
}
ngrid=100; Tvals=seq(1,110,length=ngrid);
win.graph(6,10);
par(mfrow=c(2,1),cex.axis=1.35,cex.lab=1.35,yaxs="i",xaxs="i");

N=kow.x.ZC.1.TN.r[1] #kow.x.ZC.1.TN[1] #24;
mu= mean.u #-0.1;
sigma2.b= mean.Q  #0.08^2;

xlabs=expression(paste("Projection interval ",italic(T))); ylabs="";
safe.limits=numeric(ngrid); dead.limits=numeric(ngrid)
minval=1e-16; maxval=1-minval;
for(j in 1:ngrid) {
  T=Tvals[j];
  safe.limits[j]=uniroot( function(x) {TMUpper(T,x)-0.05},
                          lower=minval,upper=maxval)$root;
  dead.limits[j]=uniroot( function(x) {TMLower(T,x)-0.95},
                          lower=minval,upper=maxval)$root;
}
matplot(Tvals,log10(cbind(safe.limits,dead.limits)),ylim=c(-
                                                             3,0),type="l",lty=1,col="white",
        xlim=c(1,100),xlab="",ylab=ylabs)
polygon(c(Tvals,rev(Tvals)),log10(c(safe.limits,rev(dead.limits))),
        col="grey85");
polygon(c(Tvals,rev(Tvals)),log10(c(dead.limits,rep(1,ngrid))),
        col="white");
logRelNe=seq(-3.125,-0.001,length=ngrid); RelNe=10^(logRelNe);
CIW=outer(Tvals,RelNe,CIWidth);
z=contourLines(Tvals,log10(RelNe),CIW,levels=0.8)
xvals=z[[1]]$x; yvals=z[[1]]$y;
polygon(c(xvals,100),c(yvals,yvals[1]),col="grey45")
abline(h=-3); abline(h=0); abline(v=100);
text(15,-1.8,"high",cex=1.25 );
text(15,-2.1,"certainty",cex=1.25 );
text(15,-2.5,expression(P[e]<0.05),cex=1.25 );
text(48,-2,expression(paste("high uncertainty about",P[e])),
     cex=1.1,srt=-41, col="white")
text(70,-1,expression(P[e]>0.95),cex=1.25,col="black");
text(70,-0.6,"high certainty",cex=1.25,col="black");
title(expression(paste("(a) ", italic(n),"= 24, ", hat(mu),
                       "= -0.1, ", sigma[b],"=0.08")));








N=25; mu= -0.1; sigma2.b=0.25^2;

safe.limits=numeric(ngrid); dead.limits=numeric(ngrid)
minval=1e-12; maxval=1-minval;
for(j in 1:ngrid) {
  T=Tvals[j];
  safe.limits[j]=uniroot( function(x) {TMUpper(T,x)-0.05},
                          lower=minval,upper=maxval)$root;
  dead.limits[j]=uniroot( function(x) {TMLower(T,x)-0.95},
                          lower=minval,upper=maxval)$root;
}
matplot(Tvals,log10(cbind(safe.limits,dead.limits)),ylim=c(-
                                                             3,0),type="l",lty=1,col="gray75",xlim=c(1,100),
        xlab=xlabs,ylab=ylabs)
polygon(c(Tvals,rev(Tvals)),log10(c(safe.limits,rev(dead.limits))),
        col="gray85");
polygon(c(Tvals,rev(Tvals)),log10(c(dead.limits,rep(1,ngrid))),
        col="white");
logRelNe=seq(-3.125,-0.001,length=ngrid); RelNe=10^(logRelNe);
CIW=outer(Tvals,RelNe,CIWidth);
z=contourLines(Tvals,log10(RelNe),CIW,levels=0.8)
xvals=z[[1]]$x; yvals=z[[1]]$y;
polygon(c(xvals,100),c(yvals,yvals[1]),col="grey45")
abline(h=-3);abline(h=0); abline(v=100);
text(10,-2.1,"high",cex=1);
text(10,-2.4,"certainty",cex=1);
text(10,-2.7,expression(P[e]<0.05),cex=1);
text(60,-1.8,"high uncertainty",cex=1.25, col="white")
text(60,-2.1,expression(paste("about ",P[e])),cex=1.25, col="white")
#text(90,0.1,expression(P[e]%~~%1),cex=1.35);
title(expression(paste("(b) ", italic(n),"= 24, ",
                       hat(mu), "= -0.1, ", sigma[b],"=0.25")));
par(mfrow=c(1,1))
