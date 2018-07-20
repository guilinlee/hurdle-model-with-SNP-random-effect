
require(doSNOW)
require(foreach)
require(rlecuyer)
source('dummy.R')
source('snp onedim.R')
source('functions for EM-AGQ.R')
source('functions for EMGH.R')
source("simulation function2.R")

# 4.b Run the analysis on all availible Cores
cluster<-makeCluster(20, type = "SOCK")
registerDoSNOW(cluster)
###prepare quadrature pointsf
load("GH.Rdata");LGh=9
absc=GH[,1];weight=GH[,2];
W=c(-1,-0.5,0,0.5,1)
lotm=200;nid=3
foreach(simiter= 1:100) %dopar% {
  simresult=simfun.lognorm(0.3,lotm,nid,0.6,0.3,1,-2.066074)
  EMdata=simresult$data
  EMdata$LOT=factor(EMdata$LOT)
  for (i in 1:lotm){
    assign(paste0("EMdata",i),subset(EMdata,EMdata$LOT==i))}
  ###define the formu=X1+X2
  formu.value=2########the corresponding column of the design matrix, write as a vector
  condi.den.type.value="bernoulli"####choosen from "bernoulli","trun.Poi","trun.NB"
  #condi.den.type.value="trun.Poi"#
  ####m is the number of lots
  Q1lower=-10;Q1upper=10
  ######parameters
  b.in=simresult$samplep.b;Gsk.in=sd(b.in)
  beta.in=0
  ################ni que ding diyibu yao b.in de initial value?
  K="0";snpftype=snp0;snpdm=snp0sdm
  snp0.in=optimfun0(b.in)
  snp.in=snp0.in$par;#####w1,w2,r,mu
  ptm=proc.time()
  result0=em.NB.approx(snp.in,beta.in,b.in,Gsk.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  time0=ptm1-ptm
  Nlik0=result0$Nloglike1+result0$Nloglike2
  para0=c(time0,Nlik0,result0$par)
  
  ptm=proc.time()
  result0GH=em.NB.GH(snp.in,beta.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  time0GH=ptm1-ptm
  Nlik0GH=result0GH$Nloglike1+result0GH$Nloglike2
  para0GH=c(time0GH,Nlik0GH,result0GH$par)
  
  ###########################
  K="1";snpftype=snp1;snpdm=snp1sdm
  snp1.in=optimfun1(b.in)
  snp.in=snp1.in$par;#####w1,w2,r,mu
  ptm=proc.time()
  result1=em.NB.approx(snp.in,beta.in,b.in,Gsk.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  time1=ptm1-ptm
  Nlik1=result1$Nloglike1+result1$Nloglike2
  para1=c(time1,Nlik1,result1$par)
  ##
  ptm=proc.time()
  result1GH=em.NB.GH(snp.in,beta.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  time1GH=ptm1-ptm
  Nlik1GH=result1GH$Nloglike1+result1GH$Nloglike2
  para1GH=c(time1GH,Nlik1GH,result1GH$par)
  
  ############################
  K="2";snpftype=snp2;snpdm=snp2sdm
  snp2.in=optimfun2(b.in)
  snp.in=snp2.in$par;#####w1,w2,r,mu
  ptm=proc.time()
  result2=em.NB.approx(snp.in,beta.in,b.in,Gsk.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  time2=ptm1-ptm
  Nlik2=result2$Nloglike1+result2$Nloglike2
  para2=c(time2,Nlik2,result2$par)
  ###
  ptm=proc.time()
  result2GH=em.NB.GH(snp.in,beta.in,maxit=1000,tol=0.001)
  ptm1=proc.time()
  Nlik2GH=result2GH$Nloglike1+result2GH$Nloglike2
  time2GH=ptm1-ptm
  para2GH=c(time2GH,Nlik2GH,result2GH$par)
  #############
  para=c(para0,para0GH,para1,para1GH,para2,para2GH)
  save(para,file=paste0("gamma2",simiter,".Rdata"))
}

