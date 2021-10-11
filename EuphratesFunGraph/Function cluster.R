requiredPackages = c("mvtnorm","reshape2","pbapply",'parallel','orthopolynom')
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}


ykmeans <- function(data,k,iter.max=10){
  
  l <- nrow(data)
  mdist <- dist(data,method = "euclidean")
  mdist <- as.numeric(mdist[1:l-1])
  smdist <- sort(mdist,decreasing = T)
  wcenter <- c()
  for (i in 1:(k+20)) {
    wcenter[i] <-  which(mdist==smdist[i])
    
  }
  wwcenter <- as.numeric(names(table(wcenter)))[1:(k-1)]
  wcenter <- wwcenter+1
  center <- data[c(1,wcenter),]
  ycl <- kmeans(data,centers = center,iter.max = iter.max)
  return(ycl)
}




get_init_par <- function(data,k){
  
  f3  <- function(y) {
    
    lop <- legendre.polynomials(n=3,normalized = F)
    lop_matrix <- as.matrix(as.data.frame(polynomial.values(polynomials = lop,x=scaleX(seq(20,120,4),u=-1,v=1))))
    colnames(lop_matrix) <- paste0("L",0:3)
    lop_fit <- lm(y~lop_matrix[,2:4]) 
    return(as.numeric(lop_fit$coefficients))
  }
  
  
  
  
  
  #????k-means???ó?ʼ????--------------------------------------
  
  init_cluster <- ykmeans(data,k,iter.max = 10)
  pro <- table(init_cluster$cluster)/nrow(data) #????ģ?͸??ʸ???
  
  cuM <- init_cluster$centers
  cusd <- diag(cov(df))
  
  init_curve_para <-  t(apply(cuM, 1, f3))
  init_sd_para <- c(mean(cusd),0.4) #SAD1
  init_pro <- pro
  
  #????mclust
  #init_cluster <- mclust::densityMclust(df,k)
  #cuM <- t(init_cluster$parameters$mean)
  #init_pro <- init_cluster$parameters$pro
  #init_sd_para <- c(mean(cusd),0.5)
  
  return_object <- list(init_sd_para,init_curve_para,init_pro)
  names(return_object)<-c("init_sd_par","init_curve","init_pro")
  return(return_object)
}


get_cluster <- function(data,k,input){
  
  requiredPackages = c("mvtnorm","reshape2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  Delta <- 10000; iter <- 0; itermax <- 2000;
  SAD1_get_matrix <- function(par,data){
    p <-  ncol(data)
    v2 <- par[1]
    phi <- par[2]
    tmp <- (1-phi^2)
    sigma <- array(dim=c(p,p))
    for(i in 1:p){
      sigma[i,i:p] <- phi^( c(i:p) - i ) * (1-phi^(2*i ))/tmp
      sigma[i:p,i] <- sigma[i,i:p]}
    sigma <- sigma*abs(v2)
    return(sigma)
  } 
  AR1_get_matrix <- function(par,data){
    n <- ncol(data)
    rho <- par[2]
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    return(par[1]*rho^exponent)
  }
  Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
  {
    u <- -1;
    v <- 1;
    if (is.null(tmin)) tmin<-min(t);
    if (is.null(tmax)) tmax<-max(t);
    ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
    np.order <- length(mu)-1;
    L <- mu[1] + ti*mu[2];
    if (np.order>=2)
      L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
    if (np.order>=3)
      L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
    if (np.order>=4)
      L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
    if (np.order>=5)
      L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
    if (np.order>=6)
      L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
    if (np.order>=7)
      L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
    if (np.order>=8)
      L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
    if (np.order>=9)
      L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
    if (np.order>=10)
      L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
    if (np.order>=11)
    {
      for(r in 11:(np.order))
      {
        kk <- ifelse(r%%2==0, r/2, (r-1)/2);
        for (k in c(0:kk) )
        {
          L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
        }
      }
    }
    return(L);
  }
  mle <- function(par,data,prob){
    par1 <- par[1:2]
    par2 <- matrix(par[-c(1:2)],nrow = k,ncol = 2)
    temp_S <- sapply(1:k, function(c) dmvnorm(data,
                                              Legendre.model(t=seq(20,120,4),mu=par2[c,]),
                                              SAD1_get_matrix(par1,data))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  while ( Delta > 1000 && iter <= itermax ) {
    # initiation
    if(iter == 0){
      init_sd_para <- input[[1]]
      init_curve_para <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_sd_para,init_curve_para)
    LL_mem <- mle(old_par,data,pro)
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             Legendre.model(t=seq(20,120,4),mu=init_curve_para[c,]),
                                             SAD1_get_matrix(init_sd_para,data))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_sd_para <- new_par$par[1:2]
    init_curve_para <- matrix(new_par$par[-c(1:2)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    #if (Delta > 500)
    #break
    cat('\n',"iter=",iter,"LL=",L_Value,"delta",Delta,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  
  cluster <- apply(omega,1,which.max)
  #clustered_df <- data.frame(cbind(row.names(data),df_fitted,cluster))
  clustered_df <- data.frame(cbind(row.names(data),df,cluster))
  colnames(clustered_df) <- c("row.names(data)",age,"cluster")
  long_df <- melt(clustered_df,id.vars=c("row.names(data)","cluster"))
  long_df[,4] <- as.numeric(long_df[,4])
  long_df[,3] <- as.numeric(long_df[,3])
  colnames(long_df) <- c("gene","cluster","time","fpkm")
  
  clustered_df <- clustered_df[,-1]
  return_object <- list(init_sd_para,init_curve_para,pro,LL_mem,BIC,clustered_df)
  names(return_object)<-c("sd_par", "curve_par", "pro", "LL", "BIC", "clustered_data")
  
  return(return_object)
  
}

t <- seq(20,120,4)
age <- t
df <- hy_genetic_effect
input=get_init_par(df,2)
cluster_par <- get_cluster(df,2,input) 
