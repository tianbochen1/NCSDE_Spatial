#######################################################################     METHOD     ###########################################################################

####################           INITIALIZATION          ################
initial=function(Xs,d=4,K=100,pen="diff",m=2,dim=length(Xs)) {
	len <- length(Xs)
	if (is.vector(Xs[[1]]) || is.ts(Xs[[1]])) {
		freq <- spec.pgram(Xs[[1]],plot=FALSE)$freq#[-c(1:(length(Xs[[1]])/100),2450:2500)];
		I <- matrix(0,nr=length(freq),nc=len)
		for (i in 1:len) I[,i] <- spec.pgram(Xs[[i]],plot=FALSE)$spec#[-c(1:(length(Xs[[1]])/100),2450:2500)]
		B <- fda::bsplineS(freq,breaks=seq(0,.5,length.out=K-d+2),norder=d); #C <- qr.Q(qr(rep(1,K)),TRUE)[,2:K]; B <- B%*%C
		if (pen=="diff") {if (m!=0) {L <- diff(diag(K),differences=m)} else {L <- diag(K)}; P <- t(L)%*%L} else 
			{P <- fda::bsplinepen(fda::create.bspline.basis(nbasis=K,norder=d))}; #P <- t(C)%*%P%*%C
	} else {
		freq <- get_Fourier_freqs(nrow(Xs[[1]]),ncol(Xs[[1]]))/2/pi; indx <- 1:nrow(freq)
		#indx <- (which(apply(abs(freq), 1, sum) == 0)+1):nrow(freq); indx <- indx[!(indx%in%unique(which(freq==0,arr.ind=1))[,1])]; freq <- freq[indx,];
		I <- matrix(0,nr=nrow(freq),nc=len)
		for (i in 1:len) I[,i] <- as.vector(fftshift(Re(periodogram2(Xs[[i]]))))[indx]
		B1 <- fda::bsplineS(freq[,1],breaks=seq(min(freq[,1]),max(freq[,1]),length.out=K-d+2),norder=d);
		B2 <- fda::bsplineS(freq[,2],breaks=seq(min(freq[,2]),max(freq[,2]),length.out=K-d+2),norder=d);
		B <- mgcv::tensor.prod.model.matrix(X=list(B1,B2))
		if (pen=="diff") {if (m!=0) {L <- diff(diag(K),differences=m)} else {L <- diag(K)}; P <- D.tensor(L1=L,L2=L)} else 
			{P1 <- P2 <- fda::bsplinepen(fda::create.bspline.basis(nbasis=K,norder=d)); P <- P1%x%diag(nrow(P2)) + diag(nrow(P1))%x%P2}; #P <- t(C)%*%P%*%C
	}; if (sum(I < .Machine$double.eps)) I[I < .Machine$double.eps] <- .Machine$double.eps
	smooth <- solve(t(B)%*%B)%*%t(B)
	Theta.tA <- (smooth%*%log(I))
	svd.Theta.tA <- svd(Theta.tA);
	k <- min(sum(svd.Theta.tA$d/svd.Theta.tA$d[1]>.Machine$double.eps),dim);
	Theta <- as.matrix(svd.Theta.tA$u[,1:k]); A <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))
	scree <- round(apply(A^2,2,sum)/sum(A^2)*100,1); A <- A[,1:k];
	A <- A%*%diag(sign(Theta[1,]),ncol(Theta)); Theta <- Theta%*%diag(sign(Theta[1,]),ncol(Theta)); 
	return(list(freq=freq, I=I, B=B, Theta=Theta, A=A, P=P, smooth=smooth, Theta.tA=Theta.tA, scree=scree, svd.Theta.tA=svd.Theta.tA))
}

#### I add KK: number of common basis
initial2=function(Xs,d=4,K=100,pen="diff",m=2,dim=length(Xs),KK) {
  len <- length(Xs)
  if (is.vector(Xs[[1]]) || is.ts(Xs[[1]])) {
    freq <- spec.pgram(Xs[[1]],plot=FALSE)$freq#[-c(1:(length(Xs[[1]])/100),2450:2500)];
    I <- matrix(0,nr=length(freq),nc=len)
    for (i in 1:len) I[,i] <- spec.pgram(Xs[[i]],plot=FALSE)$spec#[-c(1:(length(Xs[[1]])/100),2450:2500)]
    B <- fda::bsplineS(freq,breaks=seq(0,.5,length.out=K-d+2),norder=d); #C <- qr.Q(qr(rep(1,K)),TRUE)[,2:K]; B <- B%*%C
    if (pen=="diff") {if (m!=0) {L <- diff(diag(K),differences=m)} else {L <- diag(K)}; P <- t(L)%*%L} else 
    {P <- fda::bsplinepen(fda::create.bspline.basis(nbasis=K,norder=d))}; #P <- t(C)%*%P%*%C
  } else {
    freq <- get_Fourier_freqs(nrow(Xs[[1]]),ncol(Xs[[1]]))/2/pi; indx <- 1:nrow(freq)
    #indx <- (which(apply(abs(freq), 1, sum) == 0)+1):nrow(freq); indx <- indx[!(indx%in%unique(which(freq==0,arr.ind=1))[,1])]; freq <- freq[indx,];
    I <- matrix(0,nr=nrow(freq),nc=len)
    for (i in 1:len) I[,i] <- as.vector(fftshift(Re(periodogram2(Xs[[i]]))))[indx]
    B1 <- fda::bsplineS(freq[,1],breaks=seq(min(freq[,1]),max(freq[,1]),length.out=K-d+2),norder=d);
    B2 <- fda::bsplineS(freq[,2],breaks=seq(min(freq[,2]),max(freq[,2]),length.out=K-d+2),norder=d);
    B <- mgcv::tensor.prod.model.matrix(X=list(B1,B2))
    if (pen=="diff") {if (m!=0) {L <- diff(diag(K),differences=m)} else {L <- diag(K)}; P <- D.tensor(L1=L,L2=L)} else 
    {P1 <- P2 <- fda::bsplinepen(fda::create.bspline.basis(nbasis=K,norder=d)); P <- P1%x%diag(nrow(P2)) + diag(nrow(P1))%x%P2}; #P <- t(C)%*%P%*%C
  }; if (sum(I < .Machine$double.eps)) I[I < .Machine$double.eps] <- .Machine$double.eps
  smooth <- solve(t(B)%*%B)%*%t(B)
  Theta.tA <- (smooth%*%log(I))
  svd.Theta.tA <- svd(Theta.tA);
  k <- min(sum(svd.Theta.tA$d/svd.Theta.tA$d[1]>.Machine$double.eps),KK);
  Theta <- as.matrix(svd.Theta.tA$u[,1:k]); A <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))
  scree <- round(apply(A^2,2,sum)/sum(A^2)*100,1); A <- A[,1:k];
  A <- A%*%diag(sign(Theta[1,]),ncol(Theta)); Theta <- Theta%*%diag(sign(Theta[1,]),ncol(Theta)); 
  return(list(freq=freq, I=I, B=B, Theta=Theta, A=A, P=P, smooth=smooth, Theta.tA=Theta.tA, scree=scree, svd.Theta.tA=svd.Theta.tA,sv=svd.Theta.tA$d))
}# return singular values for future use


init.Theta.A.i <- function(i=1, dim) {
	svd.Theta.tA <- svd(init$Theta.tA[,i]);
	k <- min(sum(svd.Theta.tA$d/svd.Theta.tA$d[1]>.Machine$double.eps),dim);
	Theta <- as.matrix(svd.Theta.tA$u[,1:k]); A <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))
	A <- A[,1:k]; A <- A%*%diag(sign(Theta[1,]),ncol(Theta)); Theta <- Theta%*%diag(sign(Theta[1,]),ncol(Theta)); 
	return(list(Theta=Theta, A=A, svd.Theta.tA=svd.Theta.tA))
}

#######################         Evaluating the Penalized log-Likelihood     ###############
plogLike=function(Theta=init$Theta,A=init$A,lambda=rep(0,ncol(init$Theta))) {
	W <- list(); logLike <- 0;
	for (i in 1:nrow(A)) {
		W[[i]]=init$B%*%Theta%*%A[i,];
		logLike <- logLike - (sum((W[[i]]+(init$I[,i])*exp(-W[[i]])))) 
	}
	plogLike <- -2*logLike + sum(diag(lambda,ncol(Theta))%*%diag(t(Theta)%*%init$P%*%Theta))
	return(plogLike)
}
logLike=function(Theta=init$Theta,A=init$A) {
  W <- list(); logLike <- 0;
  for (i in 1:nrow(A)) {
    W[[i]]=init$B%*%Theta%*%A[i,];
    logLike <- logLike - (sum((W[[i]]+(init$I[,i])*exp(-W[[i]])))) 
  }
  loglike=-2*logLike
  return(loglike)
} 


#ploglike with penalty 2
plogLike22=function(Theta=init$Theta,A=init$A,lambda=rep(0,ncol(init$Theta)),sv=init$sv,dima=dima,lambda2) {
  W <- list(); logLike <- 0;
  for (i in 1:nrow(A)) {
    W[[i]]=init$B%*%Theta%*%A[i,];
    logLike <- logLike - (sum((W[[i]]+(init$I[,i])*exp(-W[[i]])))) 
  }
  L1 <- diff(diag(dima[1]),differences=2)
  L2 <- diff(diag(dima[2]),differences=2)
  PP <- D.tensor(L1,L2)
  K=ncol(A)
  asmo=rep(0,K)
  for(i in 1:K){ 
    asmo[i]=sv[i]/sum(sv[1:K])*t(A[,i])%*%PP%*%A[,i]
  }
  plogLike <- -2*logLike + sum(diag(lambda,ncol(Theta))%*%diag(t(Theta)%*%init$P%*%Theta))+sum(asmo)*lambda2
  return(plogLike)
}

#######################           Iterations for a given tuning parameter          #############

iterate=function(lambda=0, eps=0.01, n.iter=250, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=TRUE) {
	Theta <- list(); A <- list(); l <- 1; m <- 2
	if (!is.null(which.spec)) {
		init.i <- init.Theta.A.i(i=which.spec, dim=ncol(init$Theta));
		init$A <- init.i$A; init$Theta <- init.i$Theta; init$I <- as.matrix(init$I[,which.spec])
	}
	if (is.double(lambda)) lambda <- rep(lambda, ncol(init$Theta)); Lambda <- lambda;
	diff <- eps + nrow(init$A); Theta[[l]] <- init$Theta; A[[l]] <- init$A; AIC <- plogLike(Theta=Theta[[l]],A=A[[l]],lambda=0);
	while (abs(diff) > eps && l < n.iter) {
	   ll <- 0; repit <- TRUE; 
		A.chng <- matrix(0,nrow=nrow(A[[l]]),ncol=ncol(A[[l]])); A[[(l+1)]] <- A.chng;
		Theta.chng <- matrix(0,nrow=nrow(Theta[[l]]),ncol=ncol(Theta[[l]])); Theta[[(l+1)]] <- Theta.chng;
		temp1 <- array(0,dim=c(nrow(init$P),nrow(init$P),ncol(Theta[[l]]))); tmp1i <- temp1;
		temp2 <- array(0,dim=c(nrow(init$P),ncol(Theta[[l]])));
		for (i in 1:nrow(A[[l]])) {
			temp3 <- -(apply(init$B,2,sum)-t(init$B)%*%((init$I[,i])*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))
			temp4 <- -t(init$B*as.vector(init$I[,i]*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))%*%init$B;
			A.chng[i,] <- solve( t(Theta[[l]])%*%temp4%*%Theta[[l]])%*%(t(Theta[[l]])%*%temp3);
			for (k in 1:ncol(Theta[[l]])) {
				temp1[,,k] <- temp1[,,k] + A[[l]][i,k]^2*temp4
				temp2[,k] <- temp2[,k] + A[[l]][i,k]*temp3
				if (i==nrow(A[[l]])) {
					if (qr(temp1[,,k] - lambda[k]*init$P)$rank < nrow(temp1[,,k])) {print("choose smaller lambda"); break()}
					tmp1i[,,k] <- solve(temp1[,,k] - lambda[k]*init$P)
					Theta.chng[,k] <- tmp1i[,,k]%*%(temp2[,k]-lambda[k]*init$P%*%Theta[[l]][,k])
				}
			}
		}
		while (repit) {
			tau <- (0.5)^ll; 
			for (i in 1:nrow(A[[l]])) {A[[(l+1)]][i,] <- A[[l]][i,] - tau*A.chng[i,]}
			for (k in 1:ncol(Theta[[l]])) {Theta[[(l+1)]][,k] <- Theta[[l]][,k] - tau*Theta.chng[,k]}
			if (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) <= plogLike(Theta[[l]],A[[l]],lambda)) {repit <- FALSE} else {ll <- ll + 1}
		}
		if (trace.graph) {
			ts.plot(exp(init$B%*%Theta[[l]]%*%t(A[[l]])), gpars=list(xaxt="n",xlab=""), main=paste("Step:",l,"- lambda =",round(lambda,3)), col=ind, log="y")#, ylim=c(0,20))
			if (exists("paras")) plot.fw(paras,unique(sort(ind)),rng=range(exp(init$B%*%Theta[[l]]%*%t(A[[l]]))))
			axis(1, at = 1:length(init$freq), labels = init$freq)
		}
		diff <- plogLike(Theta[[l]],A[[l]],lambda)
		if  (update.lam==TRUE) {
			DF <- NULL; for (k in 1:ncol(Theta[[l]])) {DF[k] <- sum(diag(tmp1i[,,k]%*%temp1[,,k]))}; 
			if (combine.lam) {lambda <- rep(1/(sum(diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(sum(DF)-(m-1))), ncol(init$Theta))}
			else {lambda <- 1/((diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(DF-(m-1)))}; Lambda <- rbind(Lambda,lambda) # "ncol(Theta[[l]])" or "1"
			AIC[(l+1)] <- plogLike(Theta=Theta[[(l+1)]],A=A[[(l+1)]],lambda=0) + 2*sum(DF);
		}; diff <- diff - plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda); l <- l+1;
		# cat(paste("\n l=",(l)," half=",ll," pLogL=",round(plogLike(Theta[[(l)]],A[[(l)]],lambda))," diff=",round(diff,-log10(eps)),sep=""),"\n");# browser()
		# if  (update.lam==TRUE) {cat(paste("AIC=",AIC[(l)]),"\n");}
	}
	  svd.Theta.tA <- svd(Theta[[l]]%*%t(A[[l]]))
	  Theta[[l]] <- as.matrix(svd.Theta.tA$u[,1:ncol(Theta[[l]])]); 
	  A[[l]] <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))[,1:ncol(Theta[[l]])];
	  A[[l]] <- A[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]])); 
	  Theta[[l]] <- Theta[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]]));
		if  (update.lam!=TRUE) {
			df <- 0; for (k in 1:ncol(Theta[[l]])) {df <- df + sum(diag(tmp1i[,,k]%*%temp1[,,k]))}
			AIC <- plogLike(Theta=Theta[[l]],A=A[[(l)]],lambda=0) + 2*df
		}
	return (list(Theta=Theta,A=A,lambda=Lambda,plogLike=plogLike(Theta[[l]],A[[l]],lambda),AIC=AIC))
}


##### spatially smooth the A in each iteration. lambda_2 is chosen by roughness-lambda2 plot.
iterate2=function(lambda=0, eps=0.01, n.iter=250, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=TRUE,sv=init$sv,dima,lambda2,bw) {
  Theta <- list(); A <- list(); l <- 1; m <- 2
  if (!is.null(which.spec)) {
    init.i <- init.Theta.A.i(i=which.spec, dim=ncol(init$Theta)); ind <- ind[which.spec]
    init$A <- init.i$A; init$Theta <- init.i$Theta; init$I <- as.matrix(init$I[,which.spec])
  }
  
  
  if (is.double(lambda)) lambda <- rep(lambda, ncol(init$Theta)); Lambda <- lambda;
  diff <- eps + nrow(init$A); Theta[[l]] <- init$Theta; A[[l]] <- init$A; AIC <- plogLike2(Theta=Theta[[l]],A=A[[l]],lambda=0,sv=init$sv,dima,lambda2);
  while (abs(diff) > eps && l < n.iter) {
    ll <- 0; repit <- TRUE; 
    A.chng <- matrix(0,nrow=nrow(A[[l]]),ncol=ncol(A[[l]])); A[[(l+1)]] <- A.chng;
    Theta.chng <- matrix(0,nrow=nrow(Theta[[l]]),ncol=ncol(Theta[[l]])); Theta[[(l+1)]] <- Theta.chng;
    temp1 <- array(0,dim=c(nrow(init$P),nrow(init$P),ncol(Theta[[l]]))); tmp1i <- temp1;
    temp2 <- array(0,dim=c(nrow(init$P),ncol(Theta[[l]])));
    for (i in 1:nrow(A[[l]])) {
      temp3 <- -(apply(init$B,2,sum)-t(init$B)%*%((init$I[,i])*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))
      temp4 <- -t(init$B*as.vector(init$I[,i]*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))%*%init$B;
      A.chng[i,] <- solve(t(Theta[[l]])%*%temp4%*%Theta[[l]])%*%(t(Theta[[l]])%*%temp3);
      for (k in 1:ncol(Theta[[l]])) {
        temp1[,,k] <- temp1[,,k] + A[[l]][i,k]^2*temp4
        temp2[,k] <- temp2[,k] + A[[l]][i,k]*temp3
        if (i==nrow(A[[l]])) {
          if (qr(temp1[,,k] - lambda[k]*init$P)$rank < nrow(temp1[,,k])) {print("choose smaller lambda"); break()}
          tmp1i[,,k] <- solve(temp1[,,k] - lambda[k]*init$P)
          Theta.chng[,k] <- tmp1i[,,k]%*%(temp2[,k]-lambda[k]*init$P%*%Theta[[l]][,k])
        }
      }
    }
    while (repit) {
      tau <- (0.5)^ll; 
      for (i in 1:nrow(A[[l]])) {A[[(l+1)]][i,] <- A[[l]][i,] - tau*A.chng[i,]}
      for (k in 1:ncol(Theta[[l]])) {Theta[[(l+1)]][,k] <- Theta[[l]][,k] - tau*Theta.chng[,k]}
      if (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) <= plogLike(Theta[[l]],A[[l]],lambda)) {repit <- FALSE} else {ll <- ll + 1}
    }
    for(i in 1:ncol(A[[l]])){A[[(l+1)]][,i]=as.vector(image.smooth(matrix(A[[l]][,i],dima[1],dima[2]),theta = lambda2*(1/2)^(l-1) )$z)}##Spatially smooth columns of A
    
    if (trace.graph) {
      ts.plot(exp(init$B%*%Theta[[l]]%*%t(A[[l]])), gpars=list(xaxt="n",xlab=""), main=paste("Step:",l,"- lambda =",round(lambda,3)), col=ind, log="y")#, ylim=c(0,20))
      if (exists("paras")) plot.fw(paras,unique(sort(ind)),rng=range(exp(init$B%*%Theta[[l]]%*%t(A[[l]]))))
      axis(1, at = 1:length(init$freq), labels = init$freq)
    }
    diff <- plogLike2(Theta[[l]],A[[l]],lambda,sv=init$sv,dima,lambda2)
    if  (update.lam==TRUE) {
      DF <- NULL; for (k in 1:ncol(Theta[[l]])) {DF[k] <- sum(diag(tmp1i[,,k]%*%temp1[,,k]))}; 
      if (combine.lam) {lambda <- rep(1/(sum(diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(sum(DF)-(m-1))), ncol(init$Theta))}
      else {lambda <- 1/((diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(DF-(m-1)))}; Lambda <- rbind(Lambda,lambda) # "ncol(Theta[[l]])" or "1"
      AIC[(l+1)] <- plogLike2(Theta=Theta[[(l+1)]],A=A[[(l+1)]],lambda=0,sv=init$sv,dima,lambda2) + 2*sum(DF);
    }; diff <- diff - plogLike2(Theta[[(l+1)]],A[[(l+1)]],lambda,sv=init$sv,dima,lambda2); l <- l+1;
    cat(paste("\n l=",(l)," half=",ll," pLogL=",round(plogLike2(Theta[[(l)]],A[[(l)]],lambda,sv=init$sv,dima,lambda2))," diff=",round(diff,-log10(eps)),sep=""),"\n");# browser()
    if  (update.lam==TRUE) {cat(paste("AIC=",AIC[(l)]),"\n");}
   
    
     }
  svd.Theta.tA <- svd(Theta[[l]]%*%t(A[[l]]))
  Theta[[l]] <- as.matrix(svd.Theta.tA$u[,1:ncol(Theta[[l]])]); 
  A[[l]] <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))[,1:ncol(Theta[[l]])];
  A[[l]] <- A[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]])); 
  Theta[[l]] <- Theta[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]]));
  if  (update.lam!=TRUE) {
    df <- 0; for (k in 1:ncol(Theta[[l]])) {df <- df + sum(diag(tmp1i[,,k]%*%temp1[,,k]))}
    AIC <- plogLike2(Theta=Theta[[l]],A=A[[(l)]],lambda=0,sv=init$sv,dima,lambda2) + 2*df
  }
  return (list(Theta=Theta,A=A,lambda=Lambda,plogLike=plogLike2(Theta[[l]],A[[l]],lambda,sv=init$sv,dima,lambda2),AIC=AIC,sv=init$sv,dima,lambda2,log=logLike(Theta=Theta[[l]],A=A[[l]]),plog=plogLike(Theta=Theta[[l]],A=A[[l]],lambda)))
}

#######################################################################     END OF METHOD     ###################################################################


#######################           Iterations for a given tuning parameter          #############
iter=function(lambda=0, eps=0.01, n.iter=250, update.lam=FALSE) {
	Theta <- list(); A <- list(); #if (is.double(lambda)) lambda <- rep(lambda, ncol(init$Theta)); 
	diff <- eps + nrow(init$A); l <- 1; Theta[[l]] <- init$Theta; A[[l]] <- init$A; Lambda <- lambda;
	while (abs(diff) > eps && l < n.iter) {
	   ll <- 0; repit <- TRUE; 
		A.chng <- matrix(0,nrow=nrow(A[[l]]),ncol=ncol(A[[l]])); A[[(l+1)]] <- A.chng;
		Theta.chng <- matrix(0,nrow=nrow(Theta[[l]]),ncol=ncol(Theta[[l]])); Theta[[(l+1)]] <- Theta.chng;
		temp1 <- array(0,dim=c(nrow(init$P),nrow(init$P),ncol(Theta[[l]]))); tmp1i <- temp1;
		temp2 <- array(0,dim=c(nrow(init$P),ncol(Theta[[l]])));
		for (i in 1:nrow(A[[l]])) {
			temp3 <- -(apply(init$B,2,sum)-t(init$B)%*%((init$I[,i])*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))
			temp4 <- -t(init$B*as.vector(init$I[,i]*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))%*%init$B;
			A.chng[i,] <- solve(t(Theta[[l]])%*%(temp4-lambda*init$P)%*%Theta[[l]])%*%(t(Theta[[l]])%*%(temp3-lambda*init$P%*%Theta[[l]]%*%A[[l]][i,]));
			for (k in 1:ncol(Theta[[l]])) {
				temp1[,,k] <- temp1[,,k] + A[[l]][i,k]^2*temp4
				temp2[,k] <- temp2[,k] + A[[l]][i,k]*temp3
				if (i==nrow(A[[l]])) {
					if (qr(temp1[,,k] - lambda*apply(ans$A^2,2,sum)[k]*init$P)$rank < nrow(temp1[,,k])) {print("choose smaller lambda"); break()}
					tmp1i[,,k] <- solve(temp1[,,k] - lambda*apply(ans$A^2,2,sum)[k]*init$P)
					Theta.chng[,k] <- tmp1i[,,k]%*%(temp2[,k]-lambda*init$P%*%Theta[[l]]%*%apply(diag(A[[l]][,k],nrow(A[[l]]))%*%A[[l]],2,sum))
				}
			}
		}
		while (repit) {
			tau <- (0.5)^ll; 
			for (i in 1:nrow(A[[l]])) {A[[(l+1)]][i,] <- A[[l]][i,] - tau*A.chng[i,]}
			for (k in 1:ncol(Theta[[l]])) {Theta[[(l+1)]][,k] <- Theta[[l]][,k] - tau*Theta.chng[,k]}
			if (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) <= plogLike(Theta[[l]],A[[l]],lambda)) {repit <- FALSE} else {ll <- ll + 1}
		}
	   ts.plot(exp(init$B%*%Theta[[l]]%*%t(A[[l]])),main=paste("Step:",l,"lam:",lambda), col=ind, ylim=c(0,20))		
		diff <- (plogLike(Theta[[(l+1)]],A[[(l+1)]],lambda) - plogLike(Theta[[l]],A[[l]],lambda));
		cat(paste("\n l=",l," half=",ll," pLogL=",round(plogLike(Theta[[l+1]],A[[l+1]],lambda))," diff=",round(diff,-log10(eps)),sep=""),"\n");
		if  (update.lam==TRUE) {
			DF <- NULL; for (k in 1:ncol(Theta[[l]])) {DF[k] <- sum(diag(tmp1i[,,k]%*%temp1[,,k]))}; 
			lambda <- 1/(sum(diag(A[[(l+1)]]%*%t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]%*%t(A[[(l+1)]])))/(sum(DF)-nrow(A[[l]])+1)); Lambda <- rbind(Lambda,lambda)
		}; l <- l+1;
#	} 
	  svd.Theta.tA <- svd(Theta[[l]]%*%t(A[[l]]))
	  Theta[[l]] <- as.matrix(svd.Theta.tA$u[,1:ncol(Theta[[l]])]); 
	  A[[l]] <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))[,1:ncol(Theta[[l]])];
	  A[[l]] <- A[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]])); 
	  Theta[[l]] <- Theta[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]]));
};	df <- 0; for (k in 1:ncol(Theta[[l]])) {df <- df + sum(diag(tmp1i[,,k]%*%temp1[,,k]))}
	AIC <- plogLike(Theta=Theta[[l]],A=A[[(l)]],lambda=0) + 2*df
	return (list(Theta=Theta[[l]],A=A[[l]],lambda=Lambda,plogLike=plogLike(Theta[[l]],A[[l]],lambda),AIC=AIC))
}
#######################################################################     END OF METHOD     ###################################################################

D.tensor <- function(L1,L2) {
  d1 <- dim(L1)[2];  d2 <- dim(L2)[2]
  D1 <- t(L1)%*%L1;  D2 <- t(L2)%*%L2
  D <- kronecker(diag(1, d2),D1)+kronecker(D2,diag(1,d1))
  D
}

exp.Cov <- function (x, rho=1, sigma=1) {
	sigma^2*exp(-abs(x)/rho)
}

gauss.Cov <- function (x, rho=1, sigma=1) {
	sigma^2*exp(-abs(x)^2/(2*rho^2))
}

matern.Cov <- function (x, nu=1/2,rho=1,sigma=1) {
	(sigma^2*2^(1-nu)/gamma(nu))*(sqrt(2*nu)*abs(x)/rho)^nu*besselK(sqrt(2*nu)*abs(x)/rho, nu)
}

gauss.Spec <- function (w, rho=1, sigma=1) {
	d <- ifelse(is.matrix(w), ncol(w), 1);
	sigma^2*(2*pi*rho^2)^(d/2)*exp(-2*pi^2*rho^2*abs(w)^2)
}

matern.Spec <- function(w,nu,rho=1,sigma=1) {
	d <- ifelse(is.matrix(w), ncol(w), 1);
	sigma^2*((2*sqrt(pi))^d*gamma(nu+d/2)*(2*nu/rho^2)^nu/gamma(nu))*(2*nu/rho^2+4*pi^2*abs(w)^2)^(-nu-d/2)
}

plot.fw <- function(paras, ind=ind, rng=NULL,lwd=2,col=4) {
	f.w <- matrix(0,nrow=length(ind),ncol=length(init$freq));
	for (i in 1:length(ind)) { 
		f.w[i,] = (2*pi)*LSTS::fdensity(ar = paras[[ind[i]]], lambda = init$freq*2*pi, sd=1) 
		par(new = TRUE); if (is.null(rng)) rng <- range(f.w[i,]);
		plot(f.w[i,], axes = FALSE, xlab = '', ylab = '', log="y", lwd=lwd, col=col, type="l", ylim=rng)
	}
	return(t(f.w))
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

soilmean=function(x){
  a=x[,,1]
  for (i in 2:24){
    a=a+x[,,i]
  }
  a=a/24
  return(a)
}


bs=function(Xs,f,n,d,s){
  freq <- get_Fourier_freqs(nrow(Xs),ncol(Xs))/2/pi
  freq=freq[1:f,]
  B1 <- fda::bsplineS(freq[,1],breaks=seq(min(freq[,1]),max(freq[,1]),length.out=n),norder=d);
  B2 <- fda::bsplineS(freq[,2],breaks=seq(min(freq[,2]),max(freq[,2]),length.out=n),norder=d);
  B <- mgcv::tensor.prod.model.matrix(X=list(B1,B2))
  L1 <- diff(diag(n+d-2),differences=2)
  L2 <- diff(diag(n+d-2),differences=2)
  PP <- D.tensor(L1,L2)
  return(B%*%solve(t(B)%*%B+s*PP   )%*%t(B))
}

spa_ind=function(x,u,v){
  row=floor((x-1)/u)+1
  if (x%%u==0) col=u
  if (x%%u!=0) col=x%%u
  return(c(row,col))
}

spa_inv=function(x,y,u,v){
  return((x-1)*u+y)
}
  
spa_rough=function(A,dima){
  K=ncol(A)
  asmo=rep(0,K)
  for(p in 1:K){
    for(i in 1:dima[1]*dima[2]){
      row= spa_ind(i,dima[1],dima[2])[1]
      col= spa_ind(i,dima[1],dima[2])[2]
      if( row>1 & row<dima[2] & col>1 & col< dima[1] )
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/4*(A[spa_inv(row+1,col,dima[1],dima[2]),p]+A[spa_inv(row-1,col,dima[1],dima[2]),p]+
                                  +A[spa_inv(row,col+1,dima[1],dima[2]),p] +A[spa_inv(row,col-1,dima[1],dima[2]),p]))^2}
      if(row==1  & col==1) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/2*(A[spa_inv(row+1,col,dima[1],dima[2]),p]+A[spa_inv(row,col+1,dima[1],dima[2]),p]))^2}
      if(row==1  & col==dima[1]) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/2*(A[spa_inv(row,col-1,dima[1],dima[2]),p]+A[spa_inv(row+1,col,dima[1],dima[2]),p]))^2}
      if(row==dima[2]  & col==1) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/2*(A[spa_inv(row-1,col,dima[1],dima[2]),p]+A[spa_inv(row,col+1,dima[1],dima[2]),p]))^2}
      if(row==dima[2]  & col==dima[1]) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/2*(A[spa_inv(row-1,col,dima[1],dima[2]),p]+A[spa_inv(row,col-1,dima[1],dima[2]),p]))^2}
      if(row==1 & col>1 & col< dima[1])
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/3*(A[spa_inv(row,col-1,dima[1],dima[2]),p]+A[spa_inv(row,col+1,dima[1],dima[2]),p]
                                + +A[spa_inv(row+1,col,dima[1],dima[2]),p]))^2}
      if(row>1 & row<dima[2] & col== dima[1])
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/3*(A[spa_inv(row-1,col,dima[1],dima[2]),p]+A[spa_inv(row+1,col,dima[1],dima[2]),p]
                                + +A[spa_inv(row,col-1,dima[1],dima[2]),p]))^2}
      if(row==dima[2] & col>1 & col< dima[1]) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/3*(A[spa_inv(row,col-1,dima[1],dima[2]),p]+A[spa_inv(row,col+1,dima[1],dima[2]),p]
                                + +A[spa_inv(row-1,col,dima[1],dima[2]),p]))^2}
      if(row>1 & row<dima[2] & col==1) 
      {asmo[p]=asmo[p]+(A[spa_inv(row,col,dima[1],dima[2]),p]-
                          +1/3*(A[spa_inv(row-1,col,dima[1],dima[2]),p]+A[spa_inv(row+1,col,dima[1],dima[2]),p]
                                + +A[spa_inv(row,col+1,dima[1],dima[2]),p]))^2}
      
    }}
  return(asmo)
}


iterate13=function(lambda=0,lambda2=0, eps=0.01, n.iter=50, update.lam=TRUE, combine.lam=FALSE, which.spec=NULL, trace.graph=TRUE,dima=dima){
  #dima: dimension of A
  Theta <- list(); A <- list(); l <- 1; m <- 2
  if (!is.null(which.spec)) {
    init.i <- init.Theta.A.i(i=which.spec, dim=ncol(init$Theta)); ind <- ind[which.spec]
    init$A <- init.i$A; init$Theta <- init.i$Theta; init$I <- as.matrix(init$I[,which.spec])
  }
  if (is.double(lambda)) lambda <- rep(lambda, ncol(init$Theta)); Lambda <- lambda;
  diff <- eps + nrow(init$A); Theta[[l]] <- init$Theta; A[[l]] <- init$A;
  lambda2=rep(lambda2,ncol(A[[l]]))
  AIC <- plogLike2(Theta=Theta[[l]],A=A[[l]],lambda=0,init$sv,dima=dima,lambda2)
  lam2=lambda2
  while (abs(diff) > eps && l < n.iter) {
    ll <- 0; repit <- TRUE; 
    A.chng <- matrix(0,nrow=nrow(A[[l]]),ncol=ncol(A[[l]])); A[[(l+1)]] <- A.chng;
    
    Theta.chng <- matrix(0,nrow=nrow(Theta[[l]]),ncol=ncol(Theta[[l]])); Theta[[(l+1)]] <- Theta.chng;
    temp1 <- array(0,dim=c(nrow(init$P),nrow(init$P),ncol(Theta[[l]]))); tmp1i <- temp1;
    temp2 <- array(0,dim=c(nrow(init$P),ncol(Theta[[l]])));
    df_2=matrix(0,nrow(A[[l]]),ncol(A[[l]]))
    
    for (i in 1:nrow(A[[l]])) {
      temp3 <- -(apply(init$B,2,sum)-t(init$B)%*%((init$I[,i])*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))
      temp4 <- -t(init$B*as.vector(init$I[,i]*exp(-init$B%*%Theta[[l]]%*%t(A[[l]])[,i])))%*%init$B;
      
      for(p in 1:ncol(A[[l]])){
        temp5 = -t(init$B*as.vector(init$I[,i]*exp(-init$B%*%Theta[[l]][,p]%*%t(A[[l]])[p,i])))%*%init$B
        df_2[i,p]=solve(t(Theta[[l]][,p])%*%temp5%*%Theta[[l]][,p]-lambda2[p])%*%(t(Theta[[l]][,p])%*%temp5%*%Theta[[l]][,p])}
      row= spa_ind(i,dima[1],dima[2])[1]
      col= spa_ind(i,dima[1],dima[2])[2]
      #########PEN_2: 6 CASES
      
      ##middle cases
      if( row>2 & row<dima[2]-1 & col>2 & col< dima[1]-1 ){
        for(p in 1:ncol(A[[l]]))
        {A.chng[i,p]<-solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-10/4*lambda2[p])%*%(t(Theta[[l]][,p])%*%
        + temp3-(10/4*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-1*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]
        + +A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]
        + + A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p])  
        + +2/8*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p] 
        + + A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])
        + +2/16*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]
        + +A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p] +A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p])   )*lambda2[p])
        }}
      
      # 4 corners a11,a1n an1 ann
      if(row==1  & col==1){
        for(p in 1:ncol(A[[l]]))
        {A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-22/9*lambda2[p])%*%(t(Theta[[l]][,p])%*%
         + temp3-(22/9*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]
         + -10/6*(A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p] )
         + +4/9*(A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p] )  
         + +2/9*(A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p]))*lambda2[p])}}
      
      if(row==dima[2]  & col==1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-22/9*lambda2[p])%*%(t(Theta[[l]][,p])%*%
          + temp3-(22/9*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]
          + -10/6*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p] )
          + +4/9*(A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p] )  
          + +2/9*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p]))*lambda2[p])}}
      
      if(row==1  & col==dima[1]){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-22/9*lambda2[p])%*%(t(Theta[[l]][,p])%*%
          + temp3-(22/9*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]
          + -10/6*(A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p] )
          + +4/9*(A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p] )  
          + +2/9*(A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]))*lambda2[p])}}
      
      if(row==dima[2]  & col==dima[1]){
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-22/9*lambda2[p])%*%(t(Theta[[l]][,p])%*%
          + temp3-(22/9*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]
          + -10/6*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p] )
          + +4/9*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p] )  
          + +2/9*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]))*lambda2[p])}}
      ###a12 a21... 
	    
      if(row==1 & col==2){ 
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p])%*%(t(Theta[[l]][,p])%*%
                                                                                                   + temp3-(205/72*A[[l]][spa_inv(1,2,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(1,1,dima[1],dima[2]),p]
                                                                                                            + -4/3*A[[l]][spa_inv(1,3,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(1,4,dima[1],dima[2]),p]  
                                                                                                            +  +10/16*A[[l]][spa_inv(2,1,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(2,2,dima[1],dima[2]),p]
                                                                                                            +  +50/144*A[[l]][spa_inv(2,3,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(3,2,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==2 & col==1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]  )%*%(t(Theta[[l]][,p])%*%
                                                                                                     + temp3-(205/72*A[[l]][spa_inv(2,1,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(1,1,dima[1],dima[2]),p]
                                                                                                              + -4/3*A[[l]][spa_inv(3,1,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(4,1,dima[1],dima[2]),p]  
                                                                                                              +  +10/16*A[[l]][spa_inv(1,2,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(2,2,dima[1],dima[2]),p]
                                                                                                              +  +50/144*A[[l]][spa_inv(3,2,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(2,3,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==1 & col==dima[1]-1){
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]  )%*%(t(Theta[[l]][,p])%*%
                                                                                                     + temp3-(205/72*A[[l]][spa_inv(1,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(1,col+1,dima[1],dima[2]),p]
                                                                                                              + -4/3*A[[l]][spa_inv(1,col-1,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(1,col-2,dima[1],dima[2]),p]  
                                                                                                              +  +10/16*A[[l]][spa_inv(1,col+1,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(2,dima[1]-1,dima[1],dima[2]),p]
                                                                                                              +  +50/144*A[[l]][spa_inv(row+1,dima[1]-2,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(3,dima[1]-1,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==2 & col==dima[1]) {
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]  )%*%(t(Theta[[l]][,p])%*%
                                                                                                     + temp3-(205/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(1,col,dima[1],dima[2]),p]
                                                                                                              + -4/3*A[[l]][spa_inv(3,col,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]  
                                                                                                              +  +10/16*A[[l]][spa_inv(1,dima[1]-1,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(2,dima[1]-1,dima[1],dima[2]),p]
                                                                                                              +  +50/144*A[[l]][spa_inv(3,dima[1]-1,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(2,dima[1]-2,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==dima[2] & col==2){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p] )%*%(t(Theta[[l]][,p])%*%
                                                                                                    + temp3-(205/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(dima[2],1,dima[1],dima[2]),p]
                                                                                                             + -4/3*A[[l]][spa_inv(dima[2],3,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(dima[2],4,dima[1],dima[2]),p]  
                                                                                                             +  +10/16*A[[l]][spa_inv(dima[2]-1,1,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(dima[2]-1,2,dima[1],dima[2]),p]
                                                                                                             +  +50/144*A[[l]][spa_inv(dima[2]-1,3,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(dima[2]-2,2,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==dima[2]-1 & col==1){
        for(p in 1:ncol(A[[l]])){  #a21
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(205/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(dima[2],1,dima[1],dima[2]),p]
                                                                                                               + -4/3*A[[l]][spa_inv(dima[2]-2,1,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(dima[2]-3,1,dima[1],dima[2]),p]  
                                                                                                               +  +10/16*A[[l]][spa_inv(dima[2],2,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(dima[2]-1,2,dima[1],dima[2]),p]
                                                                                                               +  +50/144*A[[l]][spa_inv(dima[2]-2,2,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(dima[2]-1,3,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==dima[2] & col==dima[1]-1) {
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(205/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(dima[2],dima[1],dima[1],dima[2]),p]
                                                                                                               + -4/3*A[[l]][spa_inv(dima[2],dima[1]-2,dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(dima[2],dima[1]-3,dima[1],dima[2]),p]  
                                                                                                               +  +10/16*A[[l]][spa_inv(dima[2]-1,dima[1],dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(dima[2]-1,dima[1]-1,dima[1],dima[2]),p]
                                                                                                               +  +50/144*A[[l]][spa_inv(dima[2]-1,dima[1]-2,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(dima[2]-2,dima[1]-1,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==dima[2]-1 & col==dima[1]) {
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-205/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(205/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]-10/6*A[[l]][spa_inv(dima[2],dima[1],dima[1],dima[2]),p]
                                                                                                               + -4/3*A[[l]][spa_inv(dima[2]-2,dima[1],dima[1],dima[2]),p]+2/9*A[[l]][spa_inv(dima[2]-3,dima[1],dima[1],dima[2]),p]  
                                                                                                               +  +10/16*A[[l]][spa_inv(dima[2],dima[1]-1,dima[1],dima[2]),p]-14/12*A[[l]][spa_inv(dima[2]-1,dima[1]-1,dima[1],dima[2]),p]
                                                                                                               +  +50/144*A[[l]][spa_inv(dima[2]-2,dima[1]-1,dima[1],dima[2]),p]+2/16*A[[l]][spa_inv(dima[2]-1,dima[1]-2,dima[1],dima[2]),p])*lambda2[p])}}
      ######margins#######
      
      if(row==1 & col>2 & col< dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-185/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(185/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+2/9*(A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]
                                                                                                                                                                      + +A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p] )-4/3*(A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]
                                                                                                                                                                                                                             + +A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) +50/144*(A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]
                                                                                                                                                                                                                                                                                       + +A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]
                                                                                                               + +2/16*A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p])*lambda2[p])}}
      if(row>2 & row<dima[2]-1 & col== dima[1]){
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-185/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(185/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+2/9*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]
                                                                                                                                                                      + +A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p] )-4/3*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]
                                                                                                                                                                                                                             + +A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]) +50/144*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]
                                                                                                                                                                                                                                                                                       + +A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]
                                                                                                               + +2/16*A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p])*lambda2[p])}}
      if(row==dima[2] & col>2 & col< dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-185/72*lambda2[p])%*%(t(Theta[[l]][,p])%*%
                                                                                                   +temp3-(185/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+2/9*(A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]
                                                                                                                                                                  + +A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p] )-4/3*(A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]
                                                                                                                                                                                                                         + +A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) +50/144*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]
                                                                                                                                                                                                                                                                                   + + A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]
                                                                                                           + +2/16*A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p])*lambda2[p])}}
      if(row>2 & row<dima[2]-1 & col== 1){
        for(p in 1:ncol(A[[l]])){  
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-185/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      + temp3-(185/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+2/9*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]
                                                                                                                                                                      + +A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p] )-4/3*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]
                                                                                                                                                                                                                             + +A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]) +50/144*(A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]
                                                                                                                                                                                                                                                                                       + + A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]
                                                                                                               + +2/16*A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p])*lambda2[p])}}
      ###### a22#####################
      if(row==2 & col==2){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-194/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +temp3-(194/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+4/9*A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]
                                                                                                              + -14/12*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]) 
                                                                                                              + +50/144*(A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p])
                                                                                                              +  -2/2*(A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) 
                                                                                                              + +2/16*(A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p])
                                                                                                              + +2/8*A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p]  )*lambda2[p])}}
      if(row==2 & col==dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-194/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +temp3-(194/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+4/9*A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]
                                                                                                              + -14/12*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) 
                                                                                                              + +50/144*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])
                                                                                                              +  -2/2*(A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]) 
                                                                                                              + +2/16*(A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p])
                                                                                                              + +2/8*A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]  )*lambda2[p])}}
      if(row==dima[2]-1 & col==2){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-194/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +temp3-(194/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+4/9*A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]
                                                                                                              + -14/12*(A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]) 
                                                                                                              + +50/144*(A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p])
                                                                                                              +  -2/2*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) 
                                                                                                              + +2/16*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p])
                                                                                                              + +2/8*A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]  )*lambda2[p])}}
      if(row==dima[2]-1 & col==dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-194/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +temp3-(194/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+4/9*A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p]
                                                                                                              + -14/12*(A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]) 
                                                                                                              + +50/144*(A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p])
                                                                                                              +  -2/2*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]) 
                                                                                                              + +2/16*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p])
                                                                                                              +  +2/8*A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]  )*lambda2[p])}}
      #### inner margins 
      if(row==2 & col>2 & col<dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-187/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +  temp3-(187/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+25/72*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p] 
                                                                                                                                                                         +  +A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]  
                                                                                                                +  +2/16*(A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p]+A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p])  
                                                                                                                +  +2/8*(A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p]) 
                                                                                                                +  -2/2*(A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]) )*lambda2[p ])}}
      if(row==dima[2]-1 & col>2 & col<dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-187/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +  temp3-(187/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+25/72*(A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p] 
                                                                                                                                                                         +  +A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]  
                                                                                                                +  +2/16*(A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p]+A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p])  
                                                                                                                +  +2/8*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]) 
                                                                                                                +  -2/2*(A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]) )*lambda2[p])}}
      if(row>2 & row<dima[2]-1 & col==2){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-187/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +  temp3-(187/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+25/72*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p] 
                                                                                                                                                                         +  +A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]  
                                                                                                                +  +2/16*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+2,dima[1],dima[2]),p])  
                                                                                                                +  +2/8*(A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p]) 
                                                                                                                +  -2/2*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]) )*lambda2[p])}}
      if(row>2 & row<dima[2]-1 & col==dima[1]-1){
        for(p in 1:ncol(A[[l]])){
          A.chng[i,p] <- solve( t(Theta[[l]][,p])%*%temp4%*%Theta[[l]][,p]-187/72*lambda2[p]   )%*%(t(Theta[[l]][,p])%*%
                                                                                                      +  temp3-(187/72*A[[l]][spa_inv(row,col,dima[1],dima[2]),p]+25/72*(A[[l]][spa_inv(row-1,col+1,dima[1],dima[2]),p] 
                                                                                                                                                                         +  +A[[l]][spa_inv(row+1,col+1,dima[1],dima[2]),p])-14/12*A[[l]][spa_inv(row,col+1,dima[1],dima[2]),p]  
                                                                                                                +  +2/16*(A[[l]][spa_inv(row-2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row+2,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-2,dima[1],dima[2]),p])  
                                                                                                                +  +2/8*(A[[l]][spa_inv(row-1,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col-1,dima[1],dima[2]),p]) 
                                                                                                                +  -2/2*(A[[l]][spa_inv(row-1,col,dima[1],dima[2]),p]+A[[l]][spa_inv(row,col-1,dima[1],dima[2]),p]+A[[l]][spa_inv(row+1,col,dima[1],dima[2]),p]) )*lambda2[p])}}
      
      
      for (k in 1:ncol(Theta[[l]])) {
        temp1[,,k] <- temp1[,,k] + A[[l]][i,k]^2*temp4
        temp2[,k] <- temp2[,k] + A[[l]][i,k]*temp3
        if (i==nrow(A[[l]])) {
          if (qr(temp1[,,k] - lambda[k]*init$P)$rank < nrow(temp1[,,k])) {print("choose smaller lambda"); break()}
          tmp1i[,,k] <- solve(temp1[,,k] - lambda[k]*init$P)
          Theta.chng[,k] <- tmp1i[,,k]%*%(temp2[,k]-lambda[k]*init$P%*%Theta[[l]][,k])
        }
      }
    }
    
    while (repit) {
      tau <- (0.5)^ll; 
      for (i in 1:nrow(A[[l]])) {A[[(l+1)]][i,] <- A[[l]][i,] - tau*A.chng[i,]}
      for (k in 1:ncol(Theta[[l]])) {Theta[[(l+1)]][,k] <- Theta[[l]][,k] - tau*Theta.chng[,k]}
      if (plogLike2(Theta[[(l+1)]],A[[(l+1)]],lambda,init$sv,dima=dima,lambda2) <= plogLike2(Theta[[l]],A[[l]],lambda,init$sv,dima=dima,lambda2)) {repit <- FALSE} else {ll <- ll + 1}
    }
    if (trace.graph) {
      ts.plot(exp(init$B%*%Theta[[l]]%*%t(A[[l]])), gpars=list(xaxt="n",xlab=""), main=paste("Step:",l,"- lambda =",round(lambda,3)), col=ind, log="y")#, ylim=c(0,20))
      if (exists("paras")) plot.fw(paras,unique(sort(ind)),rng=range(exp(init$B%*%Theta[[l]]%*%t(A[[l]]))))
      axis(1, at = 1:length(init$freq), labels = init$freq)
    }
    diff <- plogLike2(Theta[[l]],A[[l]],lambda,init$sv,dima=dima,lambda2)
    if  (update.lam==TRUE) {
      DF <- NULL; for (k in 1:ncol(Theta[[l]])) {DF[k] <- sum(diag(tmp1i[,,k]%*%temp1[,,k]))}; 
      if (combine.lam) {lambda <- rep(1/(sum(diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(sum(DF)-(m-1))), ncol(init$Theta))}
      else {lambda <- 1/((diag(t(Theta[[(l+1)]])%*%init$P%*%Theta[[(l+1)]]))/(DF-(m-1)))}; Lambda <- rbind(Lambda,lambda) # "ncol(Theta[[l]])" or "1"
      
    } 
    
    
    if  (update.lam==TRUE) {
      DF2=rep(0,ncol(A[[l]]))
      for(k in 1:ncol(A[[l]])){
        DF2[k]=sum(df_2[,k])
        lambda2[k]=DF2[k]/spa_rough(A[[l]],dima)[k]
      }
      AIC[(l+1)] <- plogLike2(Theta=Theta[[(l+1)]],A=A[[(l+1)]],lambda=0,init$sv,dima=dima,lambda2) + 2*sum(DF)+2*sum(DF2)
      lam2=rbind(lam2,lambda2)}
    
    diff <- diff - plogLike2(Theta[[(l+1)]],A[[(l+1)]],lambda,init$sv,dima=dima,lambda2); l <- l+1;
    cat(paste("\n l=",(l)," half=",ll," pLogL=",round(plogLike2(Theta[[(l)]],A[[(l)]],lambda,init$sv,dima=dima,lambda2))," diff=",round(diff,-log10(eps)),sep=""),"\n");# browser()
    # svd.Theta.tA <- svd(Theta[[l]]%*%t(A[[l]]))
    # Theta[[l]] <- as.matrix(svd.Theta.tA$u[,1:ncol(Theta[[l]])]); 
    # A[[l]] <- as.matrix(svd.Theta.tA$v%*%diag(svd.Theta.tA$d,ncol(svd.Theta.tA$v)))[,1:ncol(Theta[[l]])];
    # A[[l]] <- A[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]])); 
    # Theta[[l]] <- Theta[[l]]%*%diag(sign(Theta[[l]][1,]),ncol(Theta[[l]]));
    if  (update.lam==TRUE) {cat(paste("AIC=",AIC[(l)]),"\n");}
  }
  return (list(Theta=Theta,A=A,lambda=Lambda,plogLike=plogLike2(Theta[[l]],A[[l]],lambda,init$sv,dima=dima,lambda2),AIC=AIC,lam2=lam2))
}
########penalized likelihood with 2 penalties
plogLike2=function(Theta=init$Theta,A=init$A,lambda=rep(0,ncol(init$Theta)),sv=init$sv,dima=dima,lambda2) {
  W <- list(); logLike <- 0;
  for (i in 1:nrow(A)) {
    W[[i]]=init$B%*%Theta%*%A[i,];
    logLike <- logLike - (sum((W[[i]]+(init$I[,i])*exp(-W[[i]])))) 
  }
  asmo=spa_rough(A,dima)
  
  plogLike <- -2*logLike + sum(diag(lambda,ncol(Theta))%*%diag(t(Theta)%*%init$P%*%Theta))+t(asmo)%*%lambda2
  return(plogLike)
}

vario =function(data){
  x = c()
  for(i in 1:40){
    x = append(x, rep(i,40))
  }
  y = rep(1:40,40)
  z = as.vector(t(data))
  raw.dat = data.frame(x=x,y=y,z=z)
  g <- gstat(formula=z~1, locations=~x+y, data=raw.dat)
  raw.vgm <- variogram(g)
  return(raw.vgm)
}


getobj = function(data){
  x = c()
  for(i in 1:40){
    x = append(x, rep(i,40))
  }
  y = rep(1:40,40)
  z = as.vector(t(data))
  raw.dat = data.frame(x=x,y=y,z=z)
  return(raw.data)
}


