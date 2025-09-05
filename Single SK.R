library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(nlme)
library(mvtnorm)
library(MASS)


cyip=function(a,b,L) #construct within-process covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}



index_sk=function(m,lensame) #creating index for sparse matrix elements
{
  p1=c();p2=c();p3 = c();p4 = c();p5 = c();p6 = c()
  
  #same column
  for(i in 1:m)
  {
    p1=c(p1,1:i)
    p2=c(p2,rep(i,length(1:i)))
  }
  i3 = m+1
  #same column
  p3=rep(1:m,lensame)
  for(i in 1:lensame)
  {
    p4=c(p4,rep(m+i,m))
  }
  for(i in i3:(i3+lensame-1))
  {
    p5=c(p5,i3:i) #not this sparse!!!!
    p6=c(p6,rep(i,length(i3:i)))
  }
  return(list(pfi=c(p1,p3,p5),pfj=c(p2,p4,p6)))
}

C_sk=function(H) #covariance matrix with sparse elements                      
{ 
  m1 = length(testx)
  K<-H     
  ppt <- m1+length(samept);
  err_mat<-matrix(0,nrow = ppt,ncol = ppt)
  #error of target matrix 
  Sigma<-matrix(K[length(K)],nrow = m1,ncol = m1)
  diag(Sigma)<-rep(1,length = m1)
  Sigma <-K[length(K)-1]*Sigma
  Sigma_t<-Sigma*rep_SK
  error_t<-Sigma_t[upper.tri(Sigma_t,diag=T)] #unscaled will be use in the same points 
  #First 30 by 30 matrix
  
  D=outer(testx,testx,`-`)
  D=D[upper.tri(D,diag=T)]
  zpp<-K[1]^2*exp(-0.25*D^2/K[2]^2) + error_t
  err_mat[1:m1,1:m1]<-Sigma_t
  
  #fill the same point on the side zsp, no error, margin from 
  zsp<-cyip(testx,xsame,H[c(1,2,5,6)])
  
  #incorporate same point from the line zss, with error 
  Dss=outer(xsame,xsame,`-`)
  Dss=Dss[upper.tri(Dss,diag=T)]
  Sigma_ss<-matrix(K[length(K)],nrow = length(xsame),ncol = length(xsame));diag(Sigma_ss)<-1
  Sigma_ss<-K[length(K)-1]*Sigma_ss;err_ss<-Sigma_ss[upper.tri(Sigma_ss,diag=T)]
  err_mat[(m1+1):ppt,(m1+1):ppt]<-Sigma_ss
  
  zss<-Reduce("+",lapply(1:2, function(i){K[2*i+1]^2*exp(-0.25*Dss^2/K[2*i+2]^2)})) + err_ss
  return(list(table = sparseMatrix(i=pfi,j=pfj,x=c(zpp,zsp,zss),symmetric=T),  #+ diag(0.0001,nrow = ppt,ncol = ppt),
              error = err_mat))
}


logL=function(H,fn) #log-likelihood function 
{ 
  #H_count<<-H_count+1
  B<<-C_sk(H)$table #after update, it become none positive definite
  #Hhist[[H_count]]<<-B
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

#scale noise for target process
rep_SK<-matrix(0,nrow = length(testx),ncol = length(testx))
table<-testtable[,!is.na(trainy3)]
for(v1 in 1:length(testx)){
  for (v2 in 1:length(testx)){
    #count the joint points 
    int_count<-sum(table[,v1] != 0 & table[,v2] != 0)
    v1_count<-sum(table[,v1] != 0)
    v2_count<-sum(table[,v2] != 0)
    rep_SK[v1,v2]<- int_count/(v1_count*v2_count)
  }
}


#recording length ans scaling parameters 
H_sk<-array(NA,dim = c(rep_time,8))
ypred_sk<-array(NA,dim = c(rep_time,length(tpoints)))
yvar_sk<-vector(mode = "list",length = rep_time)

par(mfrow = c(5,4))
par(mar = c(2, 2, 2, 2))
for(jj in 1:rep_time){
samept<-point_loc[[jj]]
xsame <- testx_star[samept]
trainsame <- pred_samples[[jj]]
ysame<-trainsame[samept]
xstar<-testx_star  #prediction points
ystar <- trainsame  #interested to know


pf=index_sk(length(testx),length(xsame))
pfi=pf$pfi;pfj=pf$pfj


y=c(testy,ysame) #list of training data
leny=length(y)
x0 = c(2,2,2,2,2,2,1,0.9) #initiate value for kernel parameters

opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000) 

pd_N<-max(c(length(testx),length(xsame)))
one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,
                    lb = c(rep(-Inf,6),0.001,0), ub = c(rep(Inf,6),Inf,0.999), opts= opts,fn= logL ), error = function(e) e)

H1=one$solution;
H0<-H1;
H_sk[jj,]<-H1

#prediction part
zip_pred=list()

zip_pred[[1]]<-cyip(testx,xstar,H0[c(1,2,5,6)])

K1=H0[3:8] #use to transfer
D2=outer(xstar,xsame,`-`)                                
zip_pred[[2]]<-t(Reduce("+",lapply(1:2, function(i){K1[2*i-1]^2*exp(-0.25*D2^2/K1[2*i]^2)}))) 
Pk=t(do.call("rbind",zip_pred)) 

D3=outer(xstar,xstar,`-`);
sigma2<-matrix(K1[length(K1)],nrow = length(xstar),ncol = length(xstar))
diag(sigma2)<-1; sigma2<-K1[length(K1)-1]*sigma2


sk=Reduce("+",lapply(1:2, function(i){K1[2*i-1]^2*exp(-0.25*D3^2/K1[2*i]^2)})) + sigma2
covM=as.matrix(C_sk(H0)$table)
raed=solve(covM,y)
ypred_correct=Pk%*%raed;
yvar1=(sk-Pk%*%solve(covM,t(Pk)))

ypred_sk[jj,]<-ypred_correct[,1]
yvar_sk[[jj]]<-yvar1

plot(x = xstar,y = ystar, type = "l",lty = 1,col = "darkgreen",ylim = c(-5,10), main = paste("Single SK  ",jj))
points(x = xsame,y = ysame, pch = 1,col = "darkgreen")
lines(x = testx,y = testy,lty = 2,col = "orange")
points(x = testx,y = testy,pch = 2,col = "orange")
lines(x = xstar,y = ypred_correct,lty = 1,col = "red")
if(jj == 1){
legend(0, 10.0, legend=c("truth", "target","predict"),
       col=c("darkgreen", "orange","red"), lty=c(1,2,1), pch = c(1,2,NA), cex=0.3)
}

#record RMSE
RMSE[jj,5]<-sqrt(sum((ystar-ypred_correct[,1])^2)/length(ystar))
}



