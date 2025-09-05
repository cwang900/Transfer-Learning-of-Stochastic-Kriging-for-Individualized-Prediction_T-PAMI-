library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(nlme)
library(mvtnorm)
library(MASS)



index=function(n,len,m,lensame) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c();p7=c();p8=c();p9=c();p10=c()
  pp=sum(len)
  for(j in 1:(n-1))
  {
    i1=1 + sum(len[0:(j-1)])
    for(i in i1:(i1+len[j]-1))
    {
      p1=c(p1,i1:i)
      p2=c(p2,rep(i,length(i1:i)))
    }
  }
  p3=rep(1:pp,m)
  for(i in 1:m)
  {
    p4=c(p4,rep(pp+i,pp))
  }
  i2=pp+1
  for(i in i2:(i2+m-1))
  {
    p5=c(p5,i2:i)
    p6=c(p6,rep(i,length(i2:i)))
  }
  #same column
  p7=rep(1:(pp+m),lensame)
  for(i in 1:lensame)
  {
    p8=c(p8,rep(pp+m+i,pp+m))
  }
  i3 = pp+m+1
  for(i in i3:(i3+lensame-1))
  {
    p9=c(p9,i3:i) #not this sparse!!!!
    p10=c(p10,rep(i,length(i3:i)))
  }
  return(list(pfi=c(p1,p3,p5,p7,p9),pfj=c(p2,p4,p6,p8,p10)))
}

cyii=function(a,b,L,i) #construct within-process covariance matrix
{ Sigma<-matrix(L[4],nrow = length(a),ncol = length(b))
diag(Sigma) <- 1; Sigma<-L[3]*Sigma*rep_factor[[i]]
d=outer(a,b,`-`)
d=d[upper.tri(d,diag=T)];error<-Sigma[upper.tri(Sigma,diag=T)]
return(list(full = L[1]^2*exp(-0.25*d^2/L[2]^2) + error, error = Sigma))
}
cyip=function(a,b,L) #construct between-process covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

#now the covaiance matrix is defined as follows

C=function(trainx,H) #between-and-within covariance matrix with sparse elements                  
{ 
  m1 = length(testx)
  pi = lengths(trainx)
  ppt = sum(pi)+m1+length(samept);
  err_mat<-matrix(0,nrow = ppt,ncol = ppt)
  
  zii=list();zip=list();zpp=c()  
  for(i in 1:(n-1)){
    i1=1 + sum(pi[0:(i-1)]);i2 = (i1+pi[i]-1)
    cyi <- cyii(trainx[[i]],trainx[[i]],L = H[c(2*i-1,2*i,4*n-1,4*n)],i)
    zii[[i]] <- cyi$full
    err_mat[i1:i2,i1:i2]<-cyi$error
  }
  
  zip = lapply(1:(n-1), function(i){cyip(trainx[[i]],testx,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  K=H[(2*n-1):(4*n)]  #9:20    
  
  Sigma<-matrix(K[length(K)],nrow = length(testx),ncol = length(testx));diag(Sigma)<-1
  Sigma<-K[length(K)-1]*Sigma*rep_factor[[n]]
  error<-Sigma[upper.tri(Sigma,diag=T)] #unscaled will be use in the same points 
  D=outer(testx,testx,`-`)
  D=D[upper.tri(D,diag=T)]
  
  zpp<-Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)})) + error
  i3<-sum(pi)+1;i4<- sum(pi)+m1;i5 = sum(pi)+m1+1; 
  err_mat[i3:i4,i3:i4]<-Sigma
  
  zsp = lapply(1:(n-1), function(i){cyip(trainx[[i]],xsame,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})

  Dd=outer(testx,xsame,`-`)
  zsp[[n]]<-Reduce("+",lapply(2:n, function(i){K[2*i-1]^2*exp(-0.25*Dd^2/K[2*i]^2)}))
    +cyip(testx,xsame,H[c(2*n-1,2*n,4*n+3,4*n+4)])
  
  Dss=outer(xsame,xsame,`-`)
  Dss=Dss[upper.tri(Dss,diag=T)]
  Sigma_ss<-matrix(K[length(K)],nrow = length(xsame),ncol = length(xsame));diag(Sigma_ss)<-1
  Sigma_ss<-K[length(K)-1]*Sigma_ss
  
  err_mat[i5:ppt,i5:ppt]<-Sigma_ss
  
  err_ss<-Sigma_ss[upper.tri(Sigma_ss,diag=T)]
  zss<-Reduce("+",lapply(2:n, function(i){K[2*i-1]^2*exp(-0.25*Dss^2/K[2*i]^2)})) + H[4*n+1]^2*exp(-0.25*Dss^2/H[4*n+2]^2) + H[4*n+3]^2*exp(-0.25*Dss^2/H[4*n+4]^2)+err_ss
  
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip)); 
  b3 = as.vector(do.call("rbind",zsp)); 
  return(list(table = sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp,b3,zss),symmetric=T), # + diag(0.00001,nrow = ppt,ncol = ppt),
              error = err_mat))
}

logL=function(H,fn) #log-likelihood function 
{ 
  B<-C(trainx,H)$table
  deter=det(B)
  if(deter>0) {a=0.5*(log(deter)+t(y)%*%solve(B,y)+log(2*pi)*leny)
  } else {
    ch=chol(B)
    logdeter=2*(sum(log(diag(ch))))
    a=0.5*(logdeter+t(y)%*%solve(B,y)+log(2*pi)*leny)
  }
  return(as.numeric(a))
}

logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}






source('TrainData.R')

#scale noise for source process
rep_factor = vector(mode = "list",length = n)  
for (i in 1:(n-1)) {
  rep_factor[[i]]<-matrix(NA,nrow = length(trainx[[i]]),ncol = length(trainx[[i]]))
  table<-ytable[[i]][,!is.na(trainy2[i,])]
  for(v1 in 1:length(trainx[[i]])){
    for (v2 in 1:length(trainx[[i]])){
      #count the joint points 
      int_count<-sum(table[,v1] != 0 & table[,v2] != 0)
      v1_count<-sum(table[,v1] != 0)
      v2_count<-sum(table[,v2] != 0)
      rep_factor[[i]][v1,v2]<- int_count/(v1_count*v2_count)
    }
  }
}


#scale noise for target process
rep_factor[[n]]<-matrix(NA,nrow = length(testx),ncol = length(testx))
table<-testtable[,!is.na(trainy3)]
for(v1 in 1:length(testx)){
  for (v2 in 1:length(testx)){
    #count the joint points 
    int_count<-sum(table[,v1] != 0 & table[,v2] != 0)
    v1_count<-sum(table[,v1] != 0)
    v2_count<-sum(table[,v2] != 0)
    rep_factor[[n]][v1,v2]<- int_count/(v1_count*v2_count)
  }
}


H_unpenal<-lapply(1:4,function(i){array(NA,dim = c(rep_time,4*n+4))})
ypred_unpenal<-array(NA,dim = c(rep_time,length(tpoints)))
yvar_unpenal<-vector(mode = "list",length = rep_time)

par(mfrow = c(5,4))
par(mar = c(2, 2, 2, 2))
for(jj in 1:rep_time){
RMSE_jj<-10^5; hot_select<-0
for (hot_key in 1:3){
samept<-point_loc[[jj]]
xsame <- testx_star[samept]
trainsame <- pred_samples[[jj]]
ysame<-trainsame[samept]
xstar<-testx_star  #prediction points
ystar <- trainsame  #interested to know

pf=index(n,lengths(trainx),m = length(testx),length(xsame))
pfi=pf$pfi;pfj=pf$pfj

y=c(unlist(trainy),testy,ysame) #list of training data
leny=length(y)

x0<-c(rep(2,4*n-2),1,0.9,2,2,2,2) #initial guess #only 23,24 not change
if (hot_key == 1){
  for(tr in 1:((n-1)/2)){x0[2*n+2*tr-1]<- 4} 
  for(tr in ((n+1)/2):(n-1)){x0[2*n+2*tr-1]<- 0}} 
else if (hot_key == 2) {
  for(tr in 1:(n-1)){x0[2*n+2*tr-1]<- 2}} 
else {
  for(tr in 1:((n-1)/2)){x0[2*n+2*tr-1]<- 0} 
  for(tr in ((n+1)/2):(n-1)){x0[2*n+2*tr-1]<- 4}} 

opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000) 

pd_N<-max(c(lengths(trainx),length(testx),length(xsame)))
one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,
                    lb = c(rep(-Inf,4*n-2),0.001,-1/pd_N,rep(-Inf,4)), ub = c(rep(Inf,4*n-2),Inf,0.999,rep(Inf,4)), opts= opts,fn= logL), error = function(e) e)
H1=one$solution
H_unpenal[[hot_key]][jj,]<-H1

#prediction part for penalized 
zip_pred=list()                                                 #1,2,11,12
zip_pred =lapply(1:(n-1), function(i){cyip(trainx[[i]],xstar,H1[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})

D1=outer(xstar,testx,`-`) 
K1=H1[(2*n-1):(4*n)]
zip_pred[[n]]=t(Reduce("+",lapply(2:n, function(i){K1[2*i-1]^2*exp(-0.25*D1^2/K1[2*i]^2)}))
                +cyip(xstar,testx,H1[c(2*n-1,2*n,4*n+3,4*n+4)]))

D11=outer(xstar,xsame,`-`) #12*3
zip_pred[[n+1]]=t(Reduce("+",lapply(2:n, function(i){K1[2*i-1]^2*exp(-0.25*D11^2/K1[2*i]^2)})) + H1[4*n+1]^2*exp(-0.25*D11^2/H1[4*n+2]^2) + H1[4*n+3]^2*exp(-0.25*D11^2/H1[4*n+4]^2))


Pk=t(do.call("rbind",zip_pred))  #get K+ 

D2=outer(xstar,xstar,`-`);
sigma2<-matrix(K1[length(K1)],nrow = length(xstar),ncol = length(xstar))
diag(sigma2)<-1; sigma2<-K1[length(K1)-1]*sigma2

sk= Reduce("+",lapply(2:n, function(i){K1[2*i-1]^2*exp(-0.25*D2^2/K1[2*i]^2)})) + H1[4*n+1]^2*exp(-0.25*D2^2/H1[4*n+2]^2) + H1[4*n+3]^2*exp(-0.25*D2^2/H1[4*n+4]^2)+ sigma2
covM=as.matrix(C(trainx,H1)$table)
raed=solve(covM,y)
ypred_hot=as.matrix(Pk%*%raed) # predicted meand
yvar_hot=as.matrix((sk-Pk%*%solve(covM,t(Pk)))) #predicted variance

RMSE_hot<-sqrt(sum((ystar-ypred_hot[,1])^2)/length(ystar))
if (RMSE_hot< RMSE_jj){
  RMSE_jj<-RMSE_hot
  ypred_correct <-ypred_hot
  yvar1<-yvar_hot
  hot_select<-hot_key
}
}
H_unpenal[[4]][jj,]<-H_unpenal[[hot_select]][jj,]
ypred_unpenal[jj,]<-ypred_correct[,1]
yvar_unpenal[[jj]]<-yvar1

plot(x = xstar,y = ystar, type = "l",lty = 1,col = "darkgreen",ylim = c(-5,10),main = paste("Proposed_unpenal ",jj))	
points(x = xsame,y = ysame, pch = 1,col = "darkgreen")	
lines(x = testx,y = testy,lty = 2,col = "orange")	
lines(x = trainx[[1]],y = trainy[[1]],col = "blue",lty = 2)	
lines(x = trainx[[2]],y = trainy[[2]],col = "dodgerblue3",lty = 2 )
lines(x = trainx[[3]],y = trainy[[3]],col = "grey",lty = 2 )
lines(x = trainx[[4]],y = trainy[[4]],col = "black",lty = 2 )


points(x = testx,y = testy,pch = 2,col = "orange",lty = 2)	
lines(x = xstar,y = ypred_correct,lty = 1,col = "red")	

#if(jj == 1){
#  legend(0, 10.0, legend=c("truth", "target sample mean","predict"),	
#         col=c("darkgreen", "orange","red"),
#         lty=c(1,2,1), pch = c(1,2,NA), cex= 0.3)}

 if(jj == 1){
   legend(0, 10.0, legend=c("truth", "target sample mean","predict","negative 1","negative 2", "positive 1","positive2"),	
          col=c("darkgreen", "orange","red","grey","black","blue","dodgerblue3"),
          lty=c(1,2,1,2,2,2,2), pch = c(1,2,NA,NA,NA,NA,NA), cex=0.3)
}

#record the RMSE
RMSE[jj,2]<-sqrt(sum((ystar-ypred_correct[,1])^2)/length(ystar))
}

