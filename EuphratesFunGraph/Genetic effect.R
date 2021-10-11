Logistic_diff = function(A,t){get_miu3(A[1:3],t)-get_miu3(A[4:6],t)}
rrloss <- function(yy,par,t){
  sum((yy-Logistic_diff(par,t))^2)
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,-1],2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))
}


get_VG <- function(marker_data,t,LR){
  
  T_marker <-marker_data#marker_data[,colnames(marker_data)%in%c(pheno_diff[,1])]
  
  diff_vg <- c() 
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(which(T_marker[a,]==1))
    aa <- as.numeric(which(T_marker[a,]==0))
    Aa <- as.numeric(which(T_marker[a,]==2))
    all <- as.numeric(which(T_marker[a,]!=9))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A鍩哄洜棰戠巼
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a鍩哄洜棰戠巼
    
    mean_AA <- Logistic_diff(LR[[a]][8:13],t)
    mean_aa <- Logistic_diff(LR[[a]][2:7],t)
    
    #mean_AA <- Logistic_diff(FunMap_par[,a][8:13],t)
    #mean_aa <- Logistic_diff(FunMap_par[,a][2:7],t)
    AE <- (mean_AA - mean_aa)
    
    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{
      mean_Aa <- Logistic_diff(LR[[a]][14:19],t)
      #mean_Aa <- Logistic_diff(FunMap_par[,a][14:19],t)
      AE <- (mean_AA - mean_aa)/2
      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE
      
    }
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  rownames(diff_vg) <- c(1:dim(marker_data)[1])
  return(sqrt(diff_vg)) 
}






hy_genetic_effect1 <- get_VG(hyhsnp,t=seq(20,120,4),LR_plastic)


