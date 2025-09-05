library(MASS)

set.seed(300)
source('fun1.R')
num <- c(5,5,5,5,2) #number of heterogeneous replication in each processes

tpoints<-seq(0,pi,length.out = 40) #input locations 
trainy1=lapply(1:(n-1),function(i){lapply(1:num[i],function(j){fun_random(tpoints,i)$output})}) #Initiate replications

#aggregate mean and average noise for all sources 
trainy2=matrix(NA, nrow=(n-1),ncol=length(tpoints)) 
ytable<-vector(mode = "list",length = (n-1))
for(i in 1:(n-1)){
  ytable[[i]]<-matrix(unlist(trainy1[[i]]),nrow = num[i],ncol = length(tpoints),byrow = TRUE)
  ptsum <- colSums(ytable[[i]]) #sum on single points 
  countpt<-apply(ytable[[i]], 2, function(c)sum(c!=0))
  for(j in 1:length(tpoints)){
    if (countpt[j] == 0){
      trainy2[i,j]<-NA
    }
    trainy2[i,j]<- ptsum[j]/countpt[j]
  }
}

#aggregate mean and average noise for the target process
testx_full <- seq(0,pi,length = 40)
trainy31=lapply(1:num[n],function(j){fun2_random(testx_full)$output}) 
testtable<-matrix(unlist(trainy31),nrow = num[n],ncol = length(testx_full),byrow = TRUE)
testsum <- colSums(testtable)
counttest<-apply(testtable, 2, function(c)sum(c!=0))
trainy3<-rep(NA,length = length(testx_full))
for(j in 1:length(testx_full)){ #mean of the target profile
  if(counttest[j] == 0){
    trainy3[j]<-NA
  }
  trainy3[j]= testsum[j]/counttest[j]
}

#filter out the null entries after taking the average
trainx<-lapply(1:(n-1),function(i){tpoints[!is.na(trainy2[i,])]})
trainy<-lapply(1:(n-1),function(i){trainy2[i,][!is.na(trainy2[i,])]})
testx<-testx_full[!is.na(trainy3)]
testy<-trainy3[!is.na(trainy3)]

testx_star<-seq(0,pi,length = 40)
pred_samples<-lapply(1:rep_time, function(j){fun2_random(testx_star,lose = FALSE)$output}) #this is a y, random sampled 
point_loc<-vector(mode = "list",length = rep_time)
for(ii in 1:rep_time){
  r2<-1:40
  point_loc[[ii]]<-sort(sample(r2[!is.na(trainy3)],10,replace = FALSE))
}


#The count matrix need to be positive definite matrix, (easy to achieve, as long as not too sparse)

#The (1,pho) structure need to have requirements as well 

#the by Schur product theorem, Hadamard product (matrices) gives a positive definite matrix 

# x0_list<-vector(mode = "list",length = ) #5 is 5 boost set
# for(i in 1:5){
#   x0_list[[i]]<-c(rep(i,4*n-2),1,0.9)
# }
# rm(n,i)

# (N-1)*(eig_1) + (eig_2) = N
# eig_1 = 1-p 
#therfore
# eig_1 = 1-p > 0, then p < 1
# eig_2 = N-(N-1)*eig_1 = N-(N-1)(1-p) = N - (N-PN-1+P) = P(N-1)+1> 0, then P> (-1)/(N-1) 

