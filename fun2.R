fun_random=function(x,i,lose = TRUE) {
  #defining all the source processes
  rho <- 0.9; n <- length(x);
  Sigma <- matrix(rho, n, n) #Equi-correlation Structure 
  diag(Sigma) <- 1
  yError<- mvrnorm(1, mu = rep(0,n), Sigma)
  if     (i==1){
    output <- 0.3*(x-3)^3 + yError}
  else if(i==2){
    output <- 0.3*x^2+2*sin(2*x) + yError}
  else if (i == 3){
    output <- (x-2)^2+ yError
  }
  else if(i == 4){
    output <- (x-1)*(x-2)*(x-4) + yError
  }
  else{
    print("only four source is designed")
  }
  complete<-output
  #cause segmentally collected data among the input locations
  if(lose == TRUE){
    r<-sample(seq(0,5,by = 1),1,replace = FALSE) 
    if(r!= 0){ output[(n-r+1):n]<-0}
    r1<-sample(25:30,1)
    pt_loc<-sample(1:n,r1,replace = FALSE) 
    output[-pt_loc]<-0
  } 
  return(list(output = output,complete = complete))
}

fun2_random=function(x,lose = TRUE){
  #defining the target processes
  rho <- 0.9; n <- length(x); 
  Sigma <- matrix(rho, n, n)
  diag(Sigma) <- 1
  yError<- mvrnorm(1, mu = rep(0,n), Sigma)
  output <- 0.2*(x-3)^3+0.15*x^2+sin(2*x) + yError
  complete<-output
  #cause segmentally collected data among the input locations
  if(lose == TRUE){
    r<-sample(seq(20,25,by = 1),1,replace = FALSE)
    output[(n-r+1):n]<-0
    r1<-sample(10:15,1)
    pt_loc<-sample(1:(n-r),r1,replace = FALSE) 
    output[-pt_loc]<-0
  }
  return(list(output = output,complete = complete))
}