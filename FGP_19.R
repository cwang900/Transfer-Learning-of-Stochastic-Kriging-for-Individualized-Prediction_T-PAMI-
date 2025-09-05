library(Matrix)
library(nloptr)
library(minqa)
library(optimx)
library(rootSolve)
library(nlme)
library(mvtnorm)
library(MASS)

#I can create len = [40,40,40,40,20]

rep_time <- 100   #number of repetition H, for RMSE
n<-5              #number of processes: 4 + 1
num <- c(5,5,5,5,2) #Number of heterogeneous replication in each processes
RMSE<-array(NA,dim = c(rep_time,6))

index=function(n,len,m) #creating index for sparse matrix elements
{ #(5,[40,40,40,40],20) sample input
  #p1,p2 are i,j for block1  
  #p3,p4 are i,j for block 2
  #p5,p6 are i,j fpr block 3,that is the (20,20)
  
  #diagonal among sources
  p1=c();p2=c(); pw1 = c(); pw2 = c(); p3=c();p4=c();p5=c();p6=c()
  pp=sum(len) #160
  for(j in 1:(n-1)) #1:4 among 5
  {
    i1=1 + sum(len[0:(j-1)]) #1,41,81,121
    for(i in i1:(i1+len[j]-1)) #1:40, 41:80,81:120,121:160, i1 = 1, iy = 1:40
    {
      p1=c(p1,i1:i) #1, 1:2,1:3,...1:40, row
      p2=c(p2,rep(i,length(i1:i))) #1 with rep 1, 2 with rep 2, 3 with rep 3 column
    }
  }
  
  #between row and column
  for (j in 1:(n-2)){ #read first line to last line
    for (i in (1+j):(n-1)){  #read first column
      pr_cum = sum(len[0:(j-1)]) #row culmulative
      pp_cum = sum(len[1:(i-1)]) #column culmulative
      
      #combine row indices
      pr = rep((pr_cum+1):(pr_cum+len[j]),len[i])         
      pw1 = c(pw1,pr)
      #combine column indices      
      pc = c()
      for (ii in 1:len[i]){
        pc = c(pc,rep(pp_cum+ii,len[j])) #41 or 81 columwise read the block 
      }
      pw2 = c(pw2,pc) #this 
    }
  }
  
  #between source and target, target to target 
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
  return(list(pfi=c(p1,pw1,p3,p5),pfj=c(p2,pw2,p4,p6)))
}

#pf is checked, then next is to fill it up 

cyii=function(a,b,L,i) #construct within-process covariance matrix
{ Sigma<-matrix(0,nrow = length(a),ncol = length(b))
diag(Sigma) = L[3]*diag(rep_factor[[i]])
d=outer(a,b,`-`)
d=d[upper.tri(d,diag=T)];error<-Sigma[upper.tri(Sigma,diag=T)]
return(list(full = L[1]^2*exp(-0.25*d^2/L[2]^2) + error,error = Sigma))
}

cyip=function(a,b,L) #construct between-process covariance matrix
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}

C=function(trainx,H) #between-and-within covariance matrix with sparse elements                    
{ 
  m1 = length(testx)
  pi = lengths(trainx)
  ppt = sum(pi)+m1
  err_mat<-matrix(0,nrow = ppt,ncol = ppt)
  
  zii=list();zss = list();zip=list();zpp=c()  
  for(i in 1:(n-1)){#1:4
    i1=1 + sum(pi[0:(i-1)]);i2 = (i1+pi[i]-1)
    cyi <- cyii(trainx[[i]],trainx[[i]],L = H[c(2*i-1,2*i,2*n+1)],i)
    zii[[i]] <- cyi$full
    err_mat[i1:i2,i1:i2]<-cyi$error
  }
  
  #add the full matrix, 4*4 contains 3+2+1 extra matrix
  read = 1
  for (j in 1:(n-2)){  #i,j in the order of reading 
    for (i in (1+j):(n-1)){
      cyss <- cyip(trainx[[j]],trainx[[i]],L = H[c(2*j-1,2*j,2*i-1,2*i)]) #totall 11 paramters when n = 5
      zss[[read]] = cyss
      read = read + 1
    }
  }

  zip = lapply(1:(n-1), function(i){cyip(trainx[[i]],testx,H[c(2*i-1,2*i,2*n-1,2*n)])})
  #K=H[(2*n-1):(4*n-1)]      
  
  Sigma<-diag(H[2*n+1],nrow = length(testx),ncol = length(testx))
  diag(Sigma)<-diag(Sigma)*diag(rep_factor[[n]]) 
  error<-Sigma[upper.tri(Sigma,diag=T)] 
  D=outer(testx,testx,`-`)
  D=D[upper.tri(D,diag=T)]
  
  zpp<- H[2*n-1]^2*exp(-0.25*D^2/H[2*n]^2) + error
  i3<-sum(pi)+1;i4<- sum(pi)+m1;
  err_mat[i3:i4,i3:i4]<-Sigma
  
  b1=unlist(zii); b2 = unlist(zss); b3=as.vector(do.call("rbind",zip)); 
  return(list(table = sparseMatrix(i=pfi,j=pfj,x=c(b1,b2,b3,zpp),symmetric=T), error = err_mat))
}




logL=function(H,fn) #log-likelihood 
{ 
  
  B=C(trainx,H)$table
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
  rep_factor[[i]]<-matrix(0,nrow = length(trainx[[i]]),ncol = length(trainx[[i]]))	
  table<-ytable[[i]][,!is.na(trainy2[i,])]	
  for(v1 in 1:length(trainx[[i]])){	
    v1_count<-sum(table[,v1] != 0)	
    rep_factor[[i]][v1,v1]<- 1/v1_count
  }	
}	

##----------------------------------------------------------------------------------
H_FGP<-array(NA,dim = c(rep_time,2*n+1))
ypred_FGP<-array(NA,dim = c(rep_time,length(tpoints)))
yvar_FGP<-vector(mode = "list",length = rep_time)

yFGP_x<-vector(mode = "list",length = rep_time) #since it includes the extra replication, so it varies 
yFGP_y<-vector(mode = "list",length = rep_time)

par(mfrow = c(4,4))
par(mar = c(2, 2, 2, 2))
for(jj in 1:100){
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
    yFGP_x[[jj]]<-testx
    yFGP_y[[jj]]<-testy
    
    #scale noise for target process
    rep_factor[[n]]<-matrix(0,nrow = length(testx),ncol = length(testx))	
    table<-testtable[,!is.na(trainy3)]	
    for(v1 in 1:length(testx)){	
      v1_count<-sum(table[,v1] != 0)	
      rep_factor[[n]][v1,v1]<- 1/v1_count
    }
    
    xstar<-testx_star  #predicted input location
    ystar <- trainsame  #underlying truth of the individual unit of interest
    
    pf=index(n,lengths(trainx),m = length(testx))
    pfi=pf$pfi;pfj=pf$pfj
    
    y=c(unlist(trainy),testy) #list of training data
    leny=length(y)
    
    x0<-c(rep(1,2*n),1) #initial guess for length and scaling parameters
    
    opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000) 
    #MLE estimation
    one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,
                        lb = c(rep(-Inf,2*n),0.001), ub = c(rep(Inf,2*n),Inf), opts= opts,fn= logL), error = function(e) e)
    H1=one$solution
    H_FGP[jj,]<-H1
    
    #Prediction for the individual units of interest
    zip_pred=list()
    zip_pred =lapply(1:(n-1), function(i){cyip(trainx[[i]],xstar,H1[c(2*i-1,2*i,2*n-1,2*n)])})
    
    D1=outer(xstar,testx,`-`) 
    zip_pred[[n]]=t(H1[(2*n-1)]^2*exp(-0.25*D1^2/H1[(2*n)]^2))
    
    
    Pk=t(do.call("rbind",zip_pred))  
    
    D2=outer(xstar,xstar,`-`);
    sigma2<-matrix(0,nrow = length(xstar),ncol = length(xstar))
    diag(sigma2)<-H1[2*n+1]
    
    
    sk= H1[(2*n-1)]^2*exp(-0.25*D2^2/H1[(2*n)]^2) + sigma2
    covM=as.matrix(C(trainx,H1)$table)
    raed=solve(covM,y)
    ypred_F=as.matrix(Pk%*%raed) # predicted meand
    yvar_F=as.matrix((sk-Pk%*%solve(covM,t(Pk)))) #predicted variance
    
    ypred_FGP[jj,] = ypred_F
    yvar_FGP[[jj]] = yvar_F
    
    plot(x = xstar,y = ystar, type = "l",lty = 1,col = "darkgreen",ylim = c(-5,10), main = paste("Full MGP",jj))
    points(x = xsame,y = ysame, pch = 16,col = "darkgreen")
    lines(x = xstar,y = ypred_F,lty = 1,col = "red")
    if(jj == 1){
      legend(0, 10.0, legend=c("truth", "predict"),
             col=c("darkgreen","red"), lty=c(1,1), pch = c(1,NA), cex=0.6)
    }
    #record the RMSE
    RMSE[jj,6]<-sqrt(sum((ystar-ypred_F[,1])^2)/length(ystar))
}