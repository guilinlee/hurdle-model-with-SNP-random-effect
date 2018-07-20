snp0sdm=function(para){
  return(para)
}
snp1sdm=function(para){
  w1=para[1];r=para[2];mu=para[3]
  a0 = sin(w1); a1 = cos(w1)
  expect_z=2*a0*a1
  expect_zsq=a0^2+3*a1^2
  var_z=expect_zsq-expect_z^2
  varb=r^2*var_z
  mb=r*expect_z+mu
  varb=as.vector(varb)
  mb=mb[[1]]
  value=c(sqrt(varb),mb)
  return(value)}

snp2sdm=function(para){
  w1=para[1];w2=para[2];r=para[3];mu=para[4]
  c1 = sin(w1); c2 = cos(w1)*sin(w2);c3=cos(w1)*cos(w2)
  a0=1.1944776*c1-0.2705981*c3; a1=c2; a2=(-0.2705981*c1+0.6532815*c3)
  expect_z = 2*a0*a1+6*a1*a2;
  expect_zsq = a0^2+3*a1^2+15*a2^2+6*a0*a2;
  var_z=expect_zsq-expect_z^2
  varb=r^2*var_z
  mb=r*expect_z+mu
  varb=varb[[1]]
  mb=mb[[1]]
  value=c(sqrt(varb),mb)
  return(value)}


log.condi=function(data,condi.parms,pb,formu,condi.den.type)
{ 
  ######dropping factor levels in a subsetted data frame in R
  data$LOT=factor(data$LOT)
  Z=dummy(data$LOT)
  
  ##compute the value of the regre.parms*regre
  np=length(condi.parms)
  regvalue=rep(0,length(data$LOT));rf=np
  while (rf>0)
  {
    form=condi.parms[rf]*data[,formu[rf]]
    regvalue=regvalue+form
    rf=rf-1
  }
  ###
  ita=as.vector(regvalue+Z%*%pb)
  y=data$Mod
  switch(condi.den.type,
         "bernoulli"={
           value=y*ita-log(1+exp(ita))
         },
         "trun.Poi"={
           lambda=as.vector(exp(ita));
           value1=dpois(y,lambda=lambda,log=TRUE)
           value= value1-log(1-exp(-lambda))
         },
         
         "trun.NB"={ lambda=exp(ita)
                     loglik0=dnbinom(y, size = theta,mu =lambda, log = TRUE)
                     loglik=loglik0-log(1-(theta/(theta+lambda))^theta)
                     
         },
         stop("Enter something that switches me!")
  )
  #value gives the density variable for any dataset data, parameter beta1, and random variable,b1 
}
####density f(yi|beta) is the product likelihood for one lot
##f(yi)is the marginal density of y. this function compute the log density value 
##of marginal density yi
##condi.parms is the parameters for conditonal density of y|b, 
##pGsk is a density transform related parameters, snp.parms is 
##parameters for density of b . absbset is the quadrature points of b
####snpfunction should have the value snptype
###output is a vector with length LGh
log.fyi=function(data,condi.parms,pGsk,snp.parms,absbset,snpfunction)
{  Qua.mar.y=c(1:LGh)  
   for(l in 1:LGh){
     absb_l=absbset[l];weight_l=weight[l];absc_l=absc[l]
     log.condi_l=sum(log.condi(data,condi.parms,absb_l,formu.value,condi.den.type.value))
     ####here, formu.value and condi.den.type.value are decided before running the code
     ##this step compute the sum log value for conditonal density of y for one lot for the kth value of absb
     log.mar.y_l=log.condi_l+snpfunction(absb_l,snp.parms)+log(pGsk)+log(weight_l/dnorm(absc_l))
     ####sum value of log likelihood of marginal y for one lot
     Qua.mar.y[l]=log.mar.y_l
   }
   return(Qua.mar.y)
}
##########weight function in Q. a vector with length LGh
#wit=f(yi|z)f(z)/f(yi) after change variable b to u, where u is stardard normal , 
##now f(u)is replacedby f(z)*weight
wit=function(data,condi.parms,pGsk,snp.parms,absbset,snpfunction)
{ weit=c(1:LGh)
  logfw=log.fyi(data,condi.parms,pGsk,snp.parms,absbset,snpfunction)
  for(l in 1:LGh)
  { weit[l]=1/(sum(exp(logfw-logfw[l])))}
  return(weit)
}
###logbayes.b
#######this function return the -loglikelihood to be minimized to compute the bayes b
logbayes.b=function(data,condi.parms,pb,snp.parms,snpfunction)
{     log.condi.y=sum(log.condi(data,condi.parms,pb,formu.value,condi.den.type.value))
      logdenyb=log.condi.y+snpfunction(pb,snp.parms)###because the first density function return the value, but snp log density return a negative value
      if(!is.finite(logdenyb)) 
        logdenyb=-1e+20
      return(-logdenyb)
}
###logbayes.bapprox
#######this function return the -loglikelihood to be minimized to compute the bayes b
###this function approximate the bayes b by omit the snpfunction(pb,snp.parms,T) 
#part in the situation that we have enough numbers samples within one lot
logbayes.bapprox=function(data,condi.parms,pb)
{ log.condi.y=sum(log.condi(data,condi.parms,pb,formu.value,condi.den.type.value))
  if(!is.finite(log.condi.y)) 
    log.condi.y=-1e+20
  return(-log.condi.y)
}

###########function for the first part of M, which is the weighted log(yi|b)
Q1=function(condi.parms,absbtotal,rwitset,tGsk)
{
  valueQ1=0
  ##begin the i step computation
  for(i in 1:lotm)
  {   
    EMdata.i=get(paste0("EMdata",i));
    Q1l=c(1:LGh)###logf(yi|bil;bt1) a vector of length LGh
    for(l in 1:LGh){
      absb_l=absbtotal[i,l]
      logfyi.bt1_l=sum(log.condi(EMdata.i,condi.parms,absb_l,formu.value,condi.den.type.value))
      Q1l[l]=logfyi.bt1_l+log(weight[l]/dnorm(absc[l]))+log(tGsk)
    }
    Q1i=Q1l%*%rwitset[i,]
    valueQ1=valueQ1+Q1i;
  }
  if(!is.finite(valueQ1)) 
    valueQ1=-1e+20
  return(-valueQ1)
}

####computations for Q2

Q2=function(snp.parms,absbtotal,rwitset,snpfunction)
  
{ 
  Q2value=0
  for(i in 1:lotm)
  {   
    logb.t1=sapply(absbtotal[i,],snpfunction,snp.parms=snp.parms)
    Q2i=logb.t1%*%rwitset[i,]
    Q2value=Q2value+Q2i
  }
  value=-Q2value
  if(!is.finite(value)) 
    value=1e+20 
  return(value)
} 

Q2optimize=function(absbtotal,rwitset,snpfunctiontype){
  switch(snpfunctiontype,
         "0"={
           est=optim(c(0.00015,0),Q2,
                     lower=c(0.0001,-1000),upper=c(1000,1000),
                     method="L-BFGS-B",absbtotal=absbtotal,rwitset=rwitset,snpfunction=snp0)
           result=list(par=est$par,lik=est$value)
         },
         "1"={sw1 = seq(-1.5, 1.5, 0.5)
              lkvalue=c(1:7);fitr=matrix(0,7,3)
              for (g in 1:7){
                gw1=sw1[g]; gc1 = sin(gw1); gc2 = cos(gw1)
                ga0=gc1;ga1=gc2
                wein=sum(rwitset)
                mb=sum(absbtotal*rwitset)/wein
                varb=sum((absbtotal-mb)^2*rwitset)/wein;
                expect_z=2*ga0*ga1
                expect_zsq=ga0^2+3*ga1^2
                var_z=expect_zsq-expect_z^2
                gr=sqrt(varb/var_z)
                gmu=mb-gr*expect_z 
                grid.op=optim(c(gw1,gr,gmu),Q2,method="L-BFGS-B",
                              lower=c(-pi/2,0.000001,-10),upper=c(pi/2,10,10),
                              absbtotal=absbtotal,rwitset=rwitset,snpfunction=snp1)
                fitr[g,]=grid.op$par; lkvalue[g]=grid.op$value
              }
              min=which.min(lkvalue); 
              result=list(par=fitr[min,],lik=lkvalue[min])
         },
         
         "2"={ sw1 = seq(-1.5, 1.5, 0.5)
               sw2 = seq(-1.5, 1.5, 0.5)
               wein=sum(rwitset)
               mb=sum(absbtotal*rwitset)/wein
               varb=sum((absbtotal-mb)^2*rwitset)/wein;
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
                 fit=optim(c(gw1,gw2,gr,gmu),Q2,method="L-BFGS-B",
                           lower=c(-pi/2,-pi/2,0.000001,-10),upper=c(pi/2,pi/2,10,10),
                           absbtotal=absbtotal,rwitset=rwitset,snpfunction=snp2)
                 fitr[g,]=fit$par; lkvalue[g]=fit$value
               }
               min=which.min(lkvalue); 
               result=list(par=fitr[min,],lik=lkvalue[min])
         },
         "3"={sw1 <- seq(-1.5, 1.5, 1)
              sw2 <- seq(-1.5, 1.5, 1)
              sw3=seq(-1.5, 1.5, 1)
              wein=sum(rwitset)
              mb=sum(absbtotal*rwitset)/wein
              varb=sum((absbtotal-mb)^2*rwitset)/wein;
              gw<- expand.grid(sw1, sw2,sw3)
              lkvalue=c(1:64);fitr=matrix(0,64,5);
              for (g in 1:64){
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
                fit=optim(c(gw1,gw2,gw2,gr,gmu),Q2,method="L-BFGS-B",
                          lower=c(-pi/2,-pi/2,-pi/2,0.000001,-10),upper=c(pi/2,pi/2,pi/2,10,10),
                          absbtotal=absbtotal,rwitset=rwitset,snpfunction=snp3)
                fitr[g,]=fit$par; lkvalue[g]=fit$value
              }
              min=which.min(lkvalue); 
              result=list(par=fitr[min,],lik=lkvalue[min])
              
         },
         "4"={sw1=seq(-1.5,1.5,1)
              sw2=seq(-1.5,1.5,1)
              sw3=seq(-1.5,1.5,1)
              sw4=seq(-1.5,1.5,1)
              gw=expand.grid(sw1,sw2,sw3,sw4)
              
              wein=sum(rwitset)
              mb=sum(absbtotal*rwitset)/wein
              varb=sum((absbtotal-mb)^2*rwitset)/wein;
              
              lkvalue=c(1:256);fitr=matrix(0,256,6)
              for (g in 1:256){
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
                fit=optim(c(gw1,gw2,gw2,gw3,gr,gmu),Q2,method="L-BFGS-B",
                          lower=c(-pi/2,-pi/2,-pi/2,-pi/2,0.000001,-10),upper=c(pi/2,pi/2,pi/2,pi/2,8,8),
                          absbtotal=absbtotal,rwitset=rwitset,snpfunction=snp4)
                fitr[g,]=fit$par; lkvalue[g]=fit$value
              }
              min=which.min(lkvalue); 
              result=list(par=fitr[min,],lik=lkvalue[min])},
         stop("Enter something that switches me!")
  )
}

em.NB <- function(snp.inits,beta.inits,b.inits,Gsk.inits,maxit,tol)
  
{
  ####### Initial parameter estimates
  flag = 0
  tsnp=snp.inits;tbeta=beta.inits;tGsk=Gsk.inits;tb=b.inits
  
  ####### Iterate between expectation and maximization steps
  for(iter in 1:maxit)
    
  {    tpar=c(tbeta,tsnp)
       ###compute the absb needed
       absbtotal=matrix(0,lotm,LGh)
       for(i in 1:lotm){absbtotal[i,]=tb[i]+tGsk*absc}
       ###compute the weight to use
       rwit=matrix(0,lotm,LGh)
       for(i in 1:lotm)
       { EMdata.i=get(paste0("EMdata",i))
         rwit[i,]=wit(EMdata.i,tbeta,tGsk,tsnp,absbtotal[i,],snpftype)}
       #####M step for Q1
       Q1new=optim(par=tbeta,fn=Q1,method="L-BFGS-B",
                   lower=Q1lower,upper=Q1upper,absbtotal=absbtotal,rwitset=rwit,tGsk=tGsk)
       ##get the Q1 paramters new beta
       t1beta=Q1new$par
       
       ####M step for Q2
       Q2new=Q2optimize(absbtotal,rwit,K)
       t1snp=Q2new$par
       ###get the Q new snp parameters
       t1par=c(t1beta,t1snp)
       
       ########## Stop iteration if the difference between the current and 
       ###new estimates is less than a tolerance level
       if( all(abs(tpar - t1par) < tol) ){ flag = 1; break;}
       
       ######## Otherwise continue iteration
       else
       {
         tbeta = t1beta; tsnp=t1snp; 
         #########update b
         ####estimate bayes estimate of bi
         b.t1=c(1:lotm)
         
         for (i in 1:lotm)
         { EMdatai=get(paste0("EMdata",i))
           fit=optimize(logbayes.bapprox,c(-1000,1000),data=EMdatai,condi.parms=tbeta)
           b.t1[i]=fit$minimum }
         tb=b.t1;tGsk=sd(tb)
       }
        print(iter);
       print(abs(tpar-t1par));
       print(t1par);
            print(tb);
       print(sd(tb));
       print(mean(tb));
       print(Q1new$value);
       print(Q2new$lik);
       
  }
  
  if(!flag) warning("Didn't converge\n")
  result=list(Nloglike1=Q1new$value,Nloglike2=Q2new$lik,par=t1par,bayeb=tb)
  return(result)
}
em.NBps <- function(snp.inits,beta.inits,b.inits,Gsk.inits,maxit,tol)
  
{
  ####### Initial parameter estimates
  flag = 0
  tsnp=snp.inits;tbeta=beta.inits;tGsk=Gsk.inits;tb=b.inits
  
  ####### Iterate between expectation and maximization steps
  for(iter in 1:maxit)
    
  {    tpar=c(tbeta,tsnp)
       ###compute the absb needed
       absbtotal=matrix(0,lotm,LGh)
       for(i in 1:lotm){absbtotal[i,]=tb[i]+tGsk*absc}
       ###compute the weight to use
       rwit=matrix(0,lotm,LGh)
       for(i in 1:lotm)
       { EMdata.i=get(paste0("EMdata",i))
         rwit[i,]=wit(EMdata.i,tbeta,tGsk,tsnp,absbtotal[i,],snpftype)}
       #####M step for Q1
       Q1new=optim(par=tbeta,fn=Q1,method="L-BFGS-B",
                   lower=Q1lower,upper=Q1upper,absbtotal=absbtotal,rwitset=rwit,tGsk=tGsk)
       ##get the Q1 paramters new beta
       t1beta=Q1new$par
       
       ####M step for Q2
       Q2new=Q2optimize(absbtotal,rwit,K)
       t1snp=Q2new$par
       ###get the Q new snp parameters
       t1par=c(t1beta,t1snp)
       
       ########## Stop iteration if the difference between the current and 
       ###new estimates is less than a tolerance level
       if( all(abs(tpar - t1par) < tol) ){ flag = 1; break;}
       
       ######## Otherwise continue iteration
       else
       {
         tbeta = t1beta; tsnp=t1snp; 
         #########update b
         ####estimate bayes estimate of bi
       }
            print(iter);
       print(abs(tpar-t1par));
       print(t1par);
             print(tb);
       print(sd(tb));
       print(mean(tb));
       print(Q1new$value);
       print(Q2new$lik);
      
  }
  
  if(!flag) warning("Didn't converge\n")
  result=list(Nloglike1=Q1new$value,Nloglike2=Q2new$lik,par=t1par,bayeb=tb)
  return(result)
}
em.NB.approx <- function(snp.inits,beta.inits,b.inits,Gsk.inits,maxit,tol)
  
{
  ####### Initial parameter estimates
  flag = 0
  tsnp=snp.inits;tbeta=beta.inits;tGsk=Gsk.inits;tb=b.inits
  
  ####### Iterate between expectation and maximization steps
  for(iter in 1:maxit)
    
  {    tpar=c(tbeta,tsnp)
       ###compute the absb needed
       absbtotal=matrix(0,lotm,LGh)
       for(i in 1:lotm){absbtotal[i,]=tb[i]+tGsk*absc}
       ###compute the weight to use
       rwit=matrix(0,lotm,LGh)
       for(i in 1:lotm)
       { EMdata.i=get(paste0("EMdata",i))
         rwit[i,]=wit(EMdata.i,tbeta,tGsk,tsnp,absbtotal[i,],snpftype)}
       #####M step for Q1
       Q1new=optim(par=tbeta,fn=Q1,method="L-BFGS-B",
                   lower=Q1lower,upper=Q1upper,absbtotal=absbtotal,rwitset=rwit,tGsk=tGsk)
       ##get the Q1 paramters new beta
       t1beta=Q1new$par
       
       ####M step for Q2
       Q2new=Q2optimize(absbtotal,rwit,K)
       t1snp=Q2new$par
       ###get the Q new snp parameters
       t1par=c(t1beta,t1snp)
       
       ########## Stop iteration if the difference between the current and 
       ###new estimates is less than a tolerance level
       if( all(abs(tpar - t1par) < tol) ){ flag = 1; break;}
       
       ######## Otherwise continue iteration
       else
       {
         tbeta = t1beta; tsnp=t1snp; 
         #########update b
         ####estimate bayes estimate of bi
         b.t1=c(1:lotm)
         baysnp0=snpdm(tsnp)
        for (i in 1:lotm)
         { EMdatai=get(paste0("EMdata",i))
           fit=optimize(logbayes.b,c(-100,100),data=EMdatai,
                         condi.parms=tbeta,snp.parms=baysnp0,snpfunction=snp0)
           b.t1[i]=fit$minimum }
         tb=b.t1;tGsk=sd(tb)
       }
       print(iter);
       print(abs(tpar-t1par));
       print(t1par);
       print(tb);
       print(sd(tb));
       print(mean(tb));
       print(Q1new$value);
       print(Q2new$lik);
       
  }
  
  if(!flag) warning("Didn't converge\n")
  result=list(Nloglike1=Q1new$value,Nloglike2=Q2new$lik,par=t1par,bayeb=tb)
  return(result)
}
