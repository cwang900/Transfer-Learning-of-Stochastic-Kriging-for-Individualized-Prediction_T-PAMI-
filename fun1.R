fun_random=function(x,i,lose = TRUE) { 
  #defining all the source processes
  rho <- 0.9; n <- length(x); 
  Sigma <- matrix(rho, n, n)
  diag(Sigma) <- 1
  yError<- mvrnorm(1, mu = rep(0,n), Sigma)   #Error Structure 
  if     (i==1){
    output <- 2*cos(5*x) + yError} 
  else if(i==2){ 
    output <- 2*sin(5*x) + yError} 
  else if(i==3){
    output <- 1/2*exp(x) + yError
  }
  else{
    output <- 1/5*(x+3)*(x+4) + yError
  }
  complete<-output
  if(lose == TRUE){
    r1<-sample(10:15,1)
    lose_pt<-sample(1:n,r1,replace = FALSE)
    output[lose_pt]<-0
  }
  return(list(output = output,complete = complete))
}

library(MASS)

fun2_random=function(x,lose = TRUE){ 
  ##defining the target processes
  rho <- 0.9; n <- length(x); 
  Sigma <- matrix(rho, n, n)
  diag(Sigma) <- 1
  yError<- mvrnorm(1, mu = rep(0,n), Sigma)
  output <- cos(5*x) + sin(5*x) + yError 
  complete<-output
  if(lose == TRUE){
    r2<-sample(10:15,1)
    #cause random sampling among the input locations
    pt_loc<-sample(1:40,r2,replace = FALSE) 
    output[-pt_loc]<-0
  }
  return(list(output = output,complete = complete))
}



