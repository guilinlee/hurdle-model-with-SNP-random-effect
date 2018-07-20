
##snp for one dimensional
#b1 w starting value
##setwd("C:/Users/A0096342/Dropbox/r")
##setwd("/Users/liguilin/Dropbox/r")
############test function
testfun=function(b,snp.parms,FUN)
{ ra=range(b)
  x=seq((ra[1]-1),(ra[2]+1),0.001)
  y=FUN(x,snp.parms,F)
  plot(x,y,type="c",xlab="Truncated NB part b (K=2)",ylab="density");c=curve(FUN(x,snp.parms,F),add=TRUE)
  lines(density(b),col="red")
}

####### SNP for k=0
###paramter conditons r>0 mu=c(-Inf,Inf)
Nsnp0<-function(b,snp.parms)
{ r=snp.parms[1];mu=snp.parms[2]
  value=dnorm(b,mu,r,log=TRUE)
  return(-sum(value)) 
}
###optim function for k=0
#this function optimize the snp0 density and use the random sample of b 
#as a input and output the optimized value of r and mu and the -loglike value in order
optimfun0=function(sam.b)
{
  est=optim(c(0,0),Nsnp0,method="L-BFGS-B",
            lower=c(0.0001,-1000),upper=c(1000,1000),b=sam.b)
  result=list(par=est$par,lik=est$value)
}
######################################################### 
###snp for k=1
  # input: b : the untransformed random variable.
         ##parms: include three parameters: -pi/2<w1<pi/2; r>0; mu
         ##output: log density of random sample b
Nsnp1<-function(b,snp.parms)###loglikelihood
{ w1=snp.parms[1];r=snp.parms[2];mu=snp.parms[3]
  if((w1<=pi/2&&w1>=(-pi/2))&&r>0)
       
  {
    c1 = sin(w1); c2 = cos(w1)
    a0=c1;a1=c2
    z=(b-mu)/r
    pk=a0+a1*z
    value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
    return(-value)
  }
  else NA
}
###optim function for k=1
#this function optimize the snp1 density and use the random sample of b 
#as a input and output the optimized value of w1,r and mu and the -loglike value in order
optimfun1=function(sam.b)
{
  sw1 <- seq(-1.5, 1.5, 0.5)
  lkvalue=c(1:7);fitr=matrix(0,7,3)
  for (g in 1:7){
    gw1=sw1[g]; ga0 = sin(gw1); ga1 = cos(gw1)
    varb=var(sam.b);mb=mean(sam.b)
    expect_z=2*ga0*ga1
    expect_zsq=ga0^2+3*ga1^2
    var_z=expect_zsq-expect_z^2
    gr=sqrt(varb/var_z)
    gmu=mb-gr*expect_z 
    grid.op=optim(c(gw1,gr,gmu),Nsnp1,method="L-BFGS-B",
                 lower=c(-pi/2,0.000001,-10),upper=c(pi/2,10,10),b=sam.b)
    fitr[g,]=grid.op$par; lkvalue[g]=grid.op$value
  }
    min=which.min(lkvalue); 
    result=list(par=fitr[min,],lik=lkvalue[min])
}
####################################################### ########
################snp for k=2
# input: b : the untransformed random variable.
##parms: include three parameters: -pi/2<w1<pi/2;-pi/2<w1<pi/2; r>0; mu
##output: log density of random sample b
Nsnp2<-function(b,snp.parms)
{ w1=snp.parms[1]
  w2=snp.parms[2]
  r=snp.parms[3]
  mu=snp.parms[4]
  
    c1 = sin(w1); c2 = cos(w1)*sin(w2);c3=cos(w1)*cos(w2)
    a0=1.1944776*c1-0.2705981*c3; 
    a1=c2; 
    a2=-0.2705981*c1+0.6532815*c3
    z=(b-mu)/r
    pk=a0+a1*z+a2*z^2

    value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
    if(!is.finite(value)) 
      value=-1e+20
    return(-value)

}
####optim function
optimfun2=function(sam.b){
sw1 <- seq(-1.5, 1.5, 0.5)
sw2 <- seq(-1.5, 1.5, 0.5)
mb=mean(sam.b)
varb=var(sam.b)
gw<- expand.grid(sw1, sw2)
lkvalue=c(1:49);fitr=matrix(0,49,4);
for (g in 1:49){
  gw1=gw[[1]][g]
  gw2=gw[[2]][g]
  c1 = sin(gw1); c2 = cos(gw1)*sin(gw2);c3=cos(gw1)*cos(gw2)
  a0=1.1944776*c1-0.2705981*c3; a1=c2; a2=(-0.2705981*c1+0.6532815*c3)
  expect_z = 2*a0*a1+6*a1*a2;
  expect_zsq = a0^2+3*a1^2+15*a2^2+6*a0*a2;
  var_z=expect_zsq-expect_z^2
  gr=sqrt(varb/var_z)
  gmu=mb-gr*expect_z   
  fit=optim(c(gw1,gw2,gr,gmu),Nsnp2,b=sam.b,method="L-BFGS-B",
            lower=c(-pi/2,-pi/2,0.000001,-10),upper=c(pi/2,pi/2,5,10))
  fitr[g,]=fit$par; lkvalue[g]=fit$value
}
min=which.min(lkvalue); 
result=list(par=fitr[min,],Nlik=lkvalue[min])
}

#################################################################
###snp for k=3
# input: b : the untransformed random variable.
##parms: include three parameters: -pi/2<w1<pi/2;-pi/2<w2,w3<pi/2; r>0; mu
##output: log density of random sample b
Nsnp3=function(b,snp.parms)
{w1=snp.parms[1]
 w2=snp.parms[2]
 w3=snp.parms[3]
 r=snp.parms[4]
 mu=snp.parms[5]
    c1 = sin(w1); c2 = cos(w1)*sin(w2); c3 = cos(w1)*cos(w2)*sin(w3);
    c4 = cos(w1)*cos(w2)*cos(w3);
    a0=1.1944776*c1-0.2705981*c3;
    a1=1.5582767*c2-0.2679064*c4;
    a2=-0.2705981*c1+0.6532815*c3;
    a3=-0.2679064*c2+0.3080468*c4
    z=(b-mu)/r
    pk=(a0+a1*z+a2*z^2+a3*z^3)
    value=(sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r)))
 if(!is.finite(value)) 
   value=-1e+20
 return(-value)
}
######optim function
optimfun3=function(sam.b)
{sw1 = seq(-1.5, 1.5, 1.5)
 sw2 = seq(-1.5, 1.5, 1.5)
 sw3=seq(-1.5, 1.5, 1.5)
 mb=mean(sam.b)
 varb=var(sam.b)
 gw=expand.grid(sw1, sw2,sw3)
 lkvalue=c(1:27);fitr=matrix(0,27,5);
 for (g in 1:27){
  gw1=gw[g,1];gw2=gw[g,2];gw3=gw[g,3]
gc1 = sin(gw1); gc2 = cos(gw1)*sin(gw2); gc3 = cos(gw1)*cos(gw2)*sin(gw3);
gc4 = cos(gw1)*cos(gw2)*cos(gw3);
ga0=1.1944776*gc1-0.2705981*gc3;
ga1=1.5582767*gc2-0.2679064*gc4;
ga2=-0.2705981*gc1+0.6532815*gc3;
ga3=-0.2679064*gc2+0.3080468*gc4
expect_z = 2*ga0*ga1+6*ga1*ga2+6*ga0*ga3+30*ga2*ga3;
expect_zsq = ga0^2+3*ga1^2+15*ga2^2+105*ga3^2+6*ga0*ga2+30*ga1*ga3;
var_z=(expect_zsq-expect_z^2)
gr=sqrt(varb/var_z)
gmu=mb-gr*expect_z 
fit=optim(c(gw1,gw2,gw3,gr,gmu),Nsnp3,b=sam.b,method="L-BFGS-B",
          lower=c(-1.5,1.5,-1.5,0.000001,-10),upper=c(1.5,1.5,1.5,1,8))
fitr[g,]=fit$par; lkvalue[g]=fit$value
 }
min=which.min(lkvalue); 
result=list(par=fitr[min,],Nlik=lkvalue[min])
}
#############################################k=4
# input: b : the untransformed random variable.
##parms: include three parameters: -pi/2<w1<pi/2;-pi/2<w2,w3,w4<pi/2; r>0; mu
##output: log density of random sample b
Nsnp4=function(b,snp.parms)
  
{
  w1=snp.parms[1];w2=snp.parms[2];w3=snp.parms[3];w4=snp.parms[4];r=snp.parms[5];mu=snp.parms[6]
  if((w1<=pi/2&&w1>=(-pi/2))&&(w2<=pi/2&&w2>=(-pi/2))&&
       (w3<=pi/2&&w3>=(-pi/2))&&(w4<=pi/2&&w4>=(-pi/2))&&r>0)
   { c1 = sin(w1);
    c2 = cos(w1)*sin(w2); c3 = cos(w1)*cos(w2)*sin(w3);
    c4 = cos(w1)*cos(w2)*cos(w3)*sin(w4); c5 = cos(w1)*cos(w2)*cos(w3)*cos(w4)
    a0=1.28273545*c1-0.4779613*c3+0.03380517*c5; a1=1.5582767*c2-0.2679064*c4; 
    a2=-0.47796128*c1+1.3210538*c3-0.16238779*c5; a3=-0.2679064*c2+ 0.3080468*c4;
    a4=0.03380517*c1-0.1623878*c3+0.11897093*c5
    z=(b-mu)/r
    pk=a0+a1*z+a2*z^2+a3*z^3+a4*z^4
    value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
    return(-value)}
 else NA
}
######grid search
optimfun4=function(sam.b)
{
  sw1=seq(-1.5,1.5,1.5)
  sw2=seq(-1.5,1.5,1.5)
  sw3=seq(-1.5,1.5,1.5)
  sw4=seq(-1.5,1.5,1.5)
  gw=expand.grid(sw1,sw2,sw3,sw4)
  mb=mean(sam.b);varb=var(sam.b)
  lkvalue=c(1:81);fitr=matrix(0,81,6)
  for (g in 1:81){
    gw1=gw[g,1];gw2=gw[g,2];gw3=gw[g,3];gw4=gw[g,4];
    c1 = sin(gw1);
    c2 = cos(gw1)*sin(gw2); c3 = cos(gw1)*cos(gw2)*sin(gw3);
    c4 = cos(gw1)*cos(gw2)*cos(gw3)*sin(gw4); c5 = cos(gw1)*cos(gw2)*cos(gw3)*cos(gw4)
    a0=1.28273545*c1-0.4779613*c3+0.03380517*c5; a1=1.5582767*c2-0.2679064*c4; 
    a2=-0.47796128*c1+1.3210538*c3-0.16238779*c5; a3=-0.2679064*c2+ 0.3080468*c4;
    a4=0.03380517*c1-0.1623878*c3+0.11897093*c5
    expect_z = 2*a0*a1+6*a1*a2+6*a0*a3+30*a2*a3+30*a1*a4+210*a3*a4;
    expect_zsq = a0^2+3*a1^2+15*a2^2+105*a3^2+945*a4^2+6*a0*a2+30*a1*a3+30*a0*a4+210*a2*a4;
    var_z=(expect_zsq-expect_z^2)
    gr=sqrt(varb/var_z)
    gmu=mb-gr*expect_z 
    fit=optim(c(gw1,gw2,gw3,gw4,gr,gmu),Nsnp4,b=sam.b,method="L-BFGS-B",
              lower=c(-1.5,-1.5,-1.5,-1.5,0.000001,-10),upper=c(1.5,1.5,1.5,1.5,8,8))
    fitr[g,]=fit$par; lkvalue[g]=fit$value
  }
  min=which.min(lkvalue); 
  result=list(par=fitr[min,],Nlik=lkvalue[min])
}
    
#####################snpBIC function compute the snp density estimation result choosn from k=0,1,2,3,4b
#####################by the critria BIC
#############BIC = k*log(n)-2logL, where k is the number of parameters, and n is the sample size
snpBIC=function(BICb)
{nBICb=length(BICb);BICp=log(nBICb)
 snpk0=optimfun0(BICb);snpk1=optimfun1(BICb);
 snpk2=optimfun2(BICb);
 BIC0=BICp*2+2*snpk0$Nlik;BIC1=BICp*3+2*snpk1$Nlik;
 BIC2=BICp*4+2*snpk2$Nlik;
 BIC=c(BIC0,BIC1,BIC2)
 BICmin=which.min(BIC);K=BICmin-1
 BICoptim=get(paste0("snpk",K))
 result=list(par=BICoptim$par,Nlik=BICoptim$Nlik,k=K)
}

#########snp function in this document return the loglikelihood instead of negative lk
snp0<-function(b,snp.parms)
{ r=snp.parms[1];mu=snp.parms[2]
  if(r>0){
  
   value=dnorm(b,mu,r,log=TRUE)
    return(sum(value)) }
  else
    NA
  
}
##
snp1<-function(b,snp.parms)###loglikelihood
{ w1=snp.parms[1];r=snp.parms[2];mu=snp.parms[3]
  if((w1<=pi/2&&w1>=(-pi/2))&&r>0)
    
  {
    c1 = sin(w1); c2 = cos(w1)
    a0=c1;a1=c2
    z=(b-mu)/r
    pk=a0+a1*z
    value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
    return(value)
    
  }
  else NA
}
##
snp2<-function(b,snp.parms)
{ w1=snp.parms[1]
  w2=snp.parms[2]
  r=snp.parms[3]
  mu=snp.parms[4]
  if((w1<=pi/2&&w1>=(-pi/2))&&(w2<=pi/2&&w2>=(-pi/2))&&r>0)
    
  {
    c1 = sin(w1); c2 = cos(w1)*sin(w2);c3=cos(w1)*cos(w2)
    a0=1.1944776*c1-0.2705981*c3; 
    a1=c2; 
    a2=-0.2705981*c1+0.6532815*c3
    z=(b-mu)/r
    pk=a0+a1*z+a2*z^2
   
    value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
     return(value)
  }
  else NA
}
##
snp3=function(b,snp.parms)
{w1=snp.parms[1]
 w2=snp.parms[2]
 w3=snp.parms[3]
 r=snp.parms[4]
 mu=snp.parms[5]
 if((w1<=pi/2&&w1>=(-pi/2))&&(w2<=pi/2&&w2>=(-pi/2))&&
      (w3<=pi/2&&w3>=(-pi/2))&&r>0)
 {
   c1 = sin(w1); c2 = cos(w1)*sin(w2); c3 = cos(w1)*cos(w2)*sin(w3);
   c4 = cos(w1)*cos(w2)*cos(w3);
   a0=1.1944776*c1-0.2705981*c3;
   a1=1.5582767*c2-0.2679064*c4;
   a2=-0.2705981*c1+0.6532815*c3;
   a3=-0.2679064*c2+0.3080468*c4
   z=(b-mu)/r
   pk=(a0+a1*z+a2*z^2+a3*z^3)
   value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
    return(value)
   
 }
 else NA
}
##
snp4=function(b,snp.parms)
  
{w1=snp.parms[1];w2=snp.parms[2];w3=snp.parms[3];w4=snp.parms[4];r=snp.parms[5];mu=snp.parms[6]
 if((w1<=pi/2&&w1>=(-pi/2))&&(w2<=pi/2&&w2>=(-pi/2))&&
      (w3<=pi/2&&w3>=(-pi/2))&&(w4<=pi/2&&w4>=(-pi/2))&&r>0)
 {
   c1 = sin(w1);
   c2 = cos(w1)*sin(w2); c3 = cos(w1)*cos(w2)*sin(w3);
   c4 = cos(w1)*cos(w2)*cos(w3)*sin(w4); c5 = cos(w1)*cos(w2)*cos(w3)*cos(w4)
   a0=1.28273545*c1-0.4779613*c3+0.03380517*c5; a1=1.5582767*c2-0.2679064*c4; 
   a2=-0.47796128*c1+1.3210538*c3-0.16238779*c5; a3=-0.2679064*c2+ 0.3080468*c4;
   a4=0.03380517*c1-0.1623878*c3+0.11897093*c5
   z=(b-mu)/r
   pk=a0+a1*z+a2*z^2+a3*z^3+a4*z^4
   value=sum(log(pk^2)-0.5*(z^2)-1/2*log(2*pi)-log(r))
   return(value)
 }
 else NA
}

##############debug
## 
trace("Nsnp3",quote(if(any(is.nan(value))) {browser()}),at=18,print=F)
