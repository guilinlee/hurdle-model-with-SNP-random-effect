simfun.norm=function(a,lotsize,nid,r1,m1){
  ni=nid*length(W)
  samplep.b=rnorm(lotsize,m1,r1)
  X1=rep(W,each=nid);X1=rep(X1,lotsize)
  n=ni*lotsize;
  logitp=matrix(0,lotsize,length(W))
  for (i in 1:lotsize){
    for (j in 1:length(W))
    {logitp[i,j]<-a*W[j]+samplep.b[i]
    }
  }
  p=1/(exp(-logitp)+1)##this is the p 
  
  y=c(1:n)
  for (i in 1:lotsize)
  {
    for (j in 1:length(W))
    {
      y[((j-1)*nid+1+(i-1)*ni):(j*nid+(i-1)*ni)]=rbinom(nid,1,p[i,j])}
  }
  sim.data=data.frame(LOT=rep(1:lotsize,each=ni),X1=X1,Mod=y)
  return(list(data=sim.data,p=p,a=a,samplep.b=samplep.b))
}

simfun.lognorm=function(a,lotsize,nid,r1,m1,br,bmu){
  ni=nid*length(W)
  samplep.z=rlnorm(lotsize,m1,r1)
  samplep.b=br*samplep.z+bmu
  X1=rep(W,each=nid);X1=rep(X1,lotsize)
  n=ni*lotsize;
  logitp=matrix(0,lotsize,length(W))
  for (i in 1:lotsize){
    for (j in 1:length(W))
    {logitp[i,j]<-a*W[j]+samplep.b[i]
    }
  }
  p=1/(exp(-logitp)+1)##this is the p 
  
  y=c(1:n)
  for (i in 1:lotsize)
  {
    for (j in 1:length(W))
    {
      y[((j-1)*nid+1+(i-1)*ni):(j*nid+(i-1)*ni)]=rbinom(nid,1,p[i,j])}
  }
  sim.data=data.frame(LOT=rep(1:lotsize,each=ni),X1=X1,Mod=y)
  return(list(data=sim.data,p=p,a=a,samplep.b=samplep.b))
}







simfun.gamma=function(a,lotsize,nid,sh,sc,r1,mu1){
  ni=nid*length(W)
  sig=rgamma(lotsize,sh,scale=sc);
  samplep.b=r1*sig+mu1     
  X1=rep(W,each=nid);X1=rep(X1,lotsize)
  n=ni*lotsize;
  logitp=matrix(0,lotsize,length(W))
  for (i in 1:lotsize){
    for (j in 1:length(W))
    {logitp[i,j]<-a*W[j]+samplep.b[i]
    }
  }
  p=1/(exp(-logitp)+1)##this is the p 
  
  y=c(1:n)
  for (i in 1:lotsize)
  {
    for (j in 1:length(W))
    {
      y[((j-1)*nid+1+(i-1)*ni):(j*nid+(i-1)*ni)]=rbinom(nid,1,p[i,j])}
  }
  sim.data=data.frame(LOT=rep(1:lotsize,each=ni),X1=X1,Mod=y)
  return(list(data=sim.data,p=p,a=a,samplep.b=samplep.b))
}





simfun.mixnorm=function(a,lotsize,nid,mr,mm1,mm2,br,bmu,mix){  
  ni=nid*length(W)             
  #Sample N random uniforms U
  U =runif(lotsize)
  #Variable to store the samples from the mixture distribution                                             
  rand.samples = rep(NA,lotsize)
  #Sampling from the mixture
  for( k in 1:lotsize){
    if(U[k]<mix){
      rand.samples[k] = rnorm(1,mm1,mr)
      ##rand.samples[k] = rnorm(1,-0.3,0.17)
    }
    else{
      rand.samples[k] = rnorm(1,mm2,mr)
      ###rand.samples[k] = rnorm(1,0.35,0.17)}
  }}
  samplep.b=br*rand.samples+bmu  
  X1=rep(W,each=nid);X1=rep(X1,lotsize)
  n=ni*lotsize;
  logitp=matrix(0,lotsize,length(W))
  for (i in 1:lotsize){
    for (j in 1:length(W))
    {logitp[i,j]<-a*W[j]+samplep.b[i]
    }
  }
  p=1/(exp(-logitp)+1)##this is the p 
  
  y=c(1:n)
  for (i in 1:lotsize)
  {
    for (j in 1:length(W))
    {
      y[((j-1)*nid+1+(i-1)*ni):(j*nid+(i-1)*ni)]=rbinom(nid,1,p[i,j])}
  }
  sim.data=data.frame(LOT=rep(1:lotsize,each=ni),X1=X1,Mod=y)
  return(list(data=sim.data,p=p,a=a,samplep.b=samplep.b))
}