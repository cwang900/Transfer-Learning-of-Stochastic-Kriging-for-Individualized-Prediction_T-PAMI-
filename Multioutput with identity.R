library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(nlme)
library(mvtnorm)
library(MASS)

index=function(n,len,m) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3=c();p4=c();p5=c();p6=c()
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
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}

cyii=function(a,b,L,i) #main diagonal elements of covariance matrix
  #here assume L[3] is sigma, L[4] is rho
{ Sigma<-matrix(0,nrow = length(a),ncol = length(b))
diag(Sigma) = L[3]*diag(rep_factor[[i]])
d=outer(a,b,`-`)
d=d[upper.tri(d,diag=T)];error<-Sigma[upper.tri(Sigma,diag=T)]
return(list(full = L[1]^2*exp(-0.25*d^2/L[2]^2) + error,error = Sigma))
}

cyip=function(a,b,L) #5 #off diagonal elements of covariance matrix, which 30*30
  #L[1] is a_ii,L[2] is b_ii, L[3] is a_in, L[4] is b_in
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

C=function(trainx,H) #covariance matrix                            
{ 
  m1 = length(testx)
  pi = lengths(trainx)
  ppt = sum(pi)+m1
  err_mat<-matrix(0,nrow = ppt,ncol = ppt)
  
  zii=list();zip=list();zpp=c()  
  for(i in 1:(n-1)){
    i1=1 + sum(pi[0:(i-1)]);i2 = (i1+pi[i]-1)
    cyi <- cyii(trainx[[i]],trainx[[i]],L = H[c(2*i-1,2*i,4*n-1)],i)
    zii[[i]] <- cyi$full
    err_mat[i1:i2,i1:i2]<-cyi$error
  }
  
  zip = lapply(1:(n-1), function(i){cyip(trainx[[i]],testx,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})
  K=H[(2*n-1):(4*n-1)]      
  
  Sigma<-diag(K[length(K)],nrow = length(testx),ncol = length(testx))
  diag(Sigma)<-diag(Sigma)*diag(rep_factor[[n]]) #scale the error by repeititon
  error<-Sigma[upper.tri(Sigma,diag=T)] #unscaled will be use in the same points 
  D=outer(testx,testx,`-`)
  D=D[upper.tri(D,diag=T)]
  
  zpp<-Reduce("+",lapply(1:n, function(i){K[2*i-1]^2*exp(-0.25*D^2/K[2*i]^2)})) + error
  i3<-sum(pi)+1;i4<- sum(pi)+m1;
  err_mat[i3:i4,i3:i4]<-Sigma
  
  b1=unlist(zii);b2=as.vector(do.call("rbind",zip)); 
  return(list(table = sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,zpp),symmetric=T), error = err_mat))
}

logL=function(H,fn) #loglikelihood 
{ 
  B=C(trainx,H)$table #after update, it become none positive definite
  deter=det(B)
  if(deter>0) {a=0.5*(log(deter)+t(y)%*%solve(B,y)+log(2*pi)*leny)
  } else {
    ch=chol(B)
    logdeter=2*(sum(log(diag(ch))))
    a=0.5*(logdeter+t(y)%*%solve(B,y)+log(2*pi)*leny)
  }
  penal = 0
  for(i in 1:(n-1)){
    if (abs(H[2*n+2*i-1]<=eta)){
      ait = H[2*n+2*i-1]
      penal = penal + ait^2/(2*eta)
    }
    else {
      penal = penal + abs(H[2*n+2*i-1]) - eta/2
    }
  }
  return(as.numeric(a) + gamma*penal)
}

logL_grad=function(H,fn)
{
  return(nl.grad(H,fn))
}

source('TrainData.R')

#scale noise for source process
rep_factor = vector(mode = "list",length = n) 
for (i in 1:(n-1)) {	
  rep_factor[[i]]<-matrix(0,nrow = length(trainx[[i]]),ncol = length(trainx[[i]]))	
  table<-ytable[[i]][,!is.na(trainy2[i,])]	
  for(v1 in 1:length(trainx[[i]])){	
      v1_count<-sum(table[,v1] != 0)	
      rep_factor[[i]][v1,v1]<- 1/v1_count
    }	
}	

##----------------------------------------------------------------------------------
H_MGP<-lapply(1:4,function(i){array(NA,dim = c(rep_time,4*n-1))})
ypred_MGP<-array(NA,dim = c(rep_time,length(tpoints)))
yvar_MGP<-vector(mode = "list",length = rep_time)

yMGP_x<-vector(mode = "list",length = rep_time)
yMGP_y<-vector(mode = "list",length = rep_time)


par(mfrow = c(5,4))
par(mar = c(2, 2, 2, 2))
for(jj in 1:rep_time){
 RMSE_jj<-10^5; hot_select<-0
 for (hot_key in 1:3){
 samept<-point_loc[[jj]]	
 xsame <- testx_star[samept] 
 trainsame <- pred_samples[[jj]]	  
 ysame<-trainsame[samept]	
  
  #add new profile into trainy31
  trainy31[[num[n]+1]]<-rep(0,length = length(testx_star))
  trainy31[[num[n]+1]][samept]<-ysame
  
  #overwrite testx and testy 
  testtable<-matrix(unlist(trainy31),nrow = num[n]+1,ncol = length(testx_full),byrow = TRUE)
  testsum <- colSums(testtable) #sum on single points 
  counttest<-apply(testtable, 2, function(c)sum(c!=0))
  trainy3<-rep(NA,length = length(testx_full))
  for(j in 1:length(testx_full)){ #mean of the target profile
    if(counttest[j] == 0){
      trainy3[j]<-NA
    }
    trainy3[j]= testsum[j]/counttest[j]
  }
  
  testx<-testx_full[!is.na(trainy3)]
  testy<-trainy3[!is.na(trainy3)]
  yMGP_x[[jj]]<-testx
  yMGP_y[[jj]]<-testy
  
  #scale noise for target process
  rep_factor[[n]]<-matrix(0,nrow = length(testx),ncol = length(testx))	
  table<-testtable[,!is.na(trainy3)]	
  for(v1 in 1:length(testx)){	
    v1_count<-sum(table[,v1] != 0)	
    rep_factor[[n]][v1,v1]<- 1/v1_count
  }
  
  xstar<-testx_star  #prediction points	
  ystar <- trainsame  #interested to know	
  
  pf=index(n,lengths(trainx),m = length(testx))
  pfi=pf$pfi;pfj=pf$pfj
  
  gamma =  0.1 #Gamma_opt
  eta = 0.00001
  
y=c(unlist(trainy),testy) #list of training data
leny=length(y)

x0<-c(rep(2,4*n-2),1) #initial guess for length and scaling parameters

if (hot_key == 1){
  for(tr in 1:((n-1)/2)){x0[2*n+2*tr-1]<- 4} 
  for(tr in ((n+1)/2):(n-1)){x0[2*n+2*tr-1]<- 0}} 
else if (hot_key == 2) {
  for(tr in 1:(n-1)){x0[2*n+2*tr-1]<- 2}} 
else {
  for(tr in 1:((n-1)/2)){x0[2*n+2*tr-1]<- 0} 
  for(tr in ((n+1)/2):(n-1)){x0[2*n+2*tr-1]<- 4}} 

opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000) 

one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,
                    lb = c(rep(-Inf,4*n-2),0.001), ub = c(rep(Inf,4*n-2),Inf), opts= opts,fn= logL ), error = function(e) e)
H1=one$solution
H_MGP[[hot_key]][jj,]<-H1

#prediction part
zip_pred=list()
zip_pred =lapply(1:(n-1), function(i){cyip(trainx[[i]],xstar,H1[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i)])})

D1=outer(xstar,testx,`-`) 
K1=H1[(2*n-1):(4*n-1)] 
zip_pred[[n]]=t(Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D1^2/K1[(2*i)]^2)})))


Pk=t(do.call("rbind",zip_pred))  

D2=outer(xstar,xstar,`-`);
sigma2<-matrix(0,nrow = length(xstar),ncol = length(xstar))
diag(sigma2)<-K1[length(K1)]


sk=Reduce("+",lapply(1:n, function(i){K1[(2*i-1)]^2*exp(-0.25*D2^2/K1[(2*i)]^2)})) + sigma2
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
H_MGP[[4]][jj,]<-H_MGP[[hot_select]][jj,]
ypred_MGP[jj,]<-ypred_correct[,1]
yvar_MGP[[jj]]<-yvar1

plot(x = xstar,y = ystar, type = "l",lty = 1,col = "darkgreen",ylim = c(-5,10),main = paste("MI uncorrelated ",jj))	
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
 
RMSE[jj,4]<-sqrt(sum((ystar-ypred_correct[,1])^2)/length(ystar))
}
