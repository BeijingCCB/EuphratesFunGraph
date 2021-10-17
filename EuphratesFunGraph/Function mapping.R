library(mvtnorm)
library(pbapply)

get_miu3 =function(B,t){B[1]/(1+exp((4*B[2]*(B[3]-t)/B[1])+2))}
SAD1_get_matrix = function(par, times = t, options=list()) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi<- par[1]   
  v2 <- par[2]   
  tmp <- (1-phi^2)   
  sigma <- array(1, dim=c(n,n))   
  for(i in 1:n)   
  {     
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp     
    sigma[i:n,i] <- sigma[i,i:n]   
  }   
  sigma <- sigma * abs(v2)   
  return(sigma); 
}
get_u = function(A,t){get_miu3(A[1:3],t)-get_miu3(A[4:6],t)}
get_sig= function(par,times = t, options=list()){SAD1_get_matrix(par = c(par[1],par[2]),times = t) + SAD1_get_matrix(par = c(par[3],par[4]),times = t) 
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,-1],2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))
}
H0 = function(yt,t,par){
  miu=get_u(par[1:6],t)
  sigma=get_sig(par = par[7:10],times=t)
  L0 = c()
  L0 = sum(dmvnorm(yt,miu,sigma,log = TRUE))
  return(-L0)
}
H1 = function(yt,t,m0,m1,m2,par,hypar){
  
  
  
  
  
  p1.0 <- yt[c(m0),]
  p1.1 <- yt[c(m1),]
  p1.2 <- yt[c(m2),]
  miu0 = get_u(par[1:6],t)
  miu1 = get_u(par[7:12],t)
  miu2 = get_u(par[13:18],t)
  sigma = get_sig(par = par[19:22],times=t)
  
  L1.0 = sum(dmvnorm(p1.0,miu0,sigma,log = TRUE))
  L1.1=sum(dmvnorm(p1.1,miu1,sigma,log = TRUE))
  L1.2=sum(dmvnorm(p1.2,miu2,sigma,log = TRUE))
  L1= L1.1+L1.0+L1.2
  return(-L1)
}
H11 = function(yt,t,m0,m1,par){
  
  
  
  
  
  p1.0 <- yt[c(m0),]
  p1.1 <- yt[c(m1),]
  
  miu0 = get_u(par[1:6],t)
  miu1 = get_u(par[7:12],t)
  sigma = get_sig(par = par[13:16],times=t)
  
  L1.1 = sum(dmvnorm(p1.0,miu0,sigma,log = TRUE))
  L1.0=sum(dmvnorm(p1.1,miu1,sigma,log = TRUE))
  L1= L1.1+L1.0
  return(-L1)
}


t <- seq(20,120,20)
parl0 <- optim(par = c(get_initial_par(cknew,t),get_initial_par(saltnew,t),0.95,0.018,0.95,0.018),H0,yt = (cknew[,-1]-saltnew[,-1]),t = t)

parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)


parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)
parl0 <- optim(par = parl0$par,H0,yt =  (cknew[,-1]-saltnew[,-1]), t = t)




optim_diff <- function(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,pheno_diff,t,m0,m1,m2){
  
  itime <- 100
  itimes <- 1
  if(length(m2)==0){
    par0 <-as.numeric(c(parl0$par[1:6],parl0$par))
    repeat{
      a <- optim(par = par0 ,H11,yt = pheno_diff, t = t,m0=m0,m1=m1)
      
      b <- optim(a$par,H11,yt = pheno_diff, t = t,m0=m0,m1=m1)
      
      
      # cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$value-b$value) < 0.01 )||itimes == itime){ #itimes瓒??鍛????瓒???
        break
      }else{
        par0 <- b$par
      }
    }
    b
  }else{
    par0 <-as.numeric(c(parl0$par[1:6],parl0$par[1:6],parl0$par))
    repeat{
      a <- optim(par = par0 ,H1,yt = pheno_diff, t = t,m0=m0,m1=m1,m2=m2)
      
      b <- optim(a$par,H1,yt = pheno_diff, t = t,m0=m0,m1=m1,m2=m2)
      
      
      # cat("Logistic_diff",itimes,b$value,'\n')
      
      itimes <- itimes + 1
      
      if(all( abs(a$value-b$value) < 0.01 )||itimes == itime){ #itimes瓒??鍛????瓒???
        break
      }else{
        par0 <- b$par
      }
    }
    b
  }
  
  
}

get_lr <- function(SNP,control,stress){
  vnp=SNP
  
  m0 <- which( vnp == 0)
  m1 <- which( vnp == 1)
  m2 <- which( vnp == 2)
  pheno_ck <- control[,-1]
  pheno_salt <- stress[,-1]
  pheno_ck0 <- pheno_ck[m0,]
  pheno_ck1 <- pheno_ck[m1,]
  pheno_salt0 <- pheno_salt[m0,]
  pheno_salt1 <- pheno_salt[m1,]
  
  parl1 <- optim_diff(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,(control[,-1]-stress[,-1]),t,m0,m1,m2)
  LR_plastic3 <- c(-2*(-parl0$value+parl1$value),parl1$par)
}


LR_plastic <- apply(hyhsnp,1,get_lr,control=cknew,stress=saltnew)



lr <- c()
for (i in 1:1000) {
  m <- sample(1:76,76)
  
  cknew <- cknew[m,]
  saltnew <- saltnew[m,]
  
  LR_plasticper <- apply(hyhsnp,1,get_lr,control=cknew,stress=saltnew)
  LR_plasticper2 <- matrix(NA,ncol=length(LR_plasticper),nrow=2)
  for (i in 1:length(LR_plasticper)) {
    LR_plasticper2[1,i] <- LR_plasticper[[i]][[1]]
  }
  lr[i] <- max(LR_plasticper2[1,])
  
}
Threshold_value <- sort(lr,decreasing = T)[50]
