"closedp.Mtb" <- function(X,dfreq=FALSE)
{

        X<-as.matrix(X)
        t <- ifelse(dfreq,dim(X)[2]-1,dim(X)[2])

    #####################################################################################################################################
    # Validation des arguments fournis en entrée
    
    # Argument dfreq
    if(!is.logical(dfreq)||!isTRUE(all.equal(length(dfreq),1))) stop("'dfreq' must be a logical object of length 1")
    
    # Argument X
    if (dfreq)
    {
        if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("Every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,t+1]%%1)!=0)) stop("The last column of 'X' must contain capture histories frequencies, therefore integers")
    } else {
        if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
    
    #####################################################################################################################################

        Y <- histfreq.t(X,dfreq=dfreq)
        histpos <- histpos.t(t)
        nbcap <- apply(histpos, 1, sum)
              
        pnames <- vector("character",t)
        for (i in 1:t) { pnames[i]<-paste("p",i,sep="") }
        cnames <- vector("character",t-1)
        for (i in 2:t) { cnames[i-1]<-paste("c",i,sep="") }


        # modele Mtb
        I1 <- matrix(c(0,0,1),ncol=1)  # matrices d'indicatices de non capture (wi=0) avant la première capture, i=1,...,t-1
        for(i in (3:t))   I1 <- cbind( c(rep(0, 2^(i-1)), rep(1,((2^(i-1))-1))), rbind(matrix(rep(0,(2^(i-1))*(i-2)),nrow=2^(i-1)),I1))       
        I2 <- matrix(rep(0,(2^t-1)*t),ncol=t)  # matrice d'indicatrices de la première capture
        poscapt1<-rep(1,2^(t-1))
        for(i in 2:t) poscapt1 <- c(poscapt1,rep(i,2^(t-i)))
        for(i in 1:(2^t-1)) I2[i,poscapt1[i]] <- 1        
        I3 <- histpos.t(t)[,-1]  # matrice d'indicatrices de capture (wi=1) après la première capture
        for(i in (2^(t-1)+1):(2^t-1)) I3[i,1:(poscapt1[i]-1)]<-rep(0,length(I3[i,1:(poscapt1[i]-1)]))
        I4 <- 1-histpos.t(t)[,-1]  # matrice d'indicatrices de non-capture (wi=0) après la première capture
        for(i in (2^(t-1)+1):(2^t-1)) I4[i,1:(poscapt1[i]-1)]<-rep(0,length(I4[i,1:(poscapt1[i]-1)]))
        # valeurs initiales des paramètres (N->n, tous les p -> 0.5 et ci=pi)
        init_loglinparam.Mtb <- as.vector(c(log(sum(Y)),rep(0,t+1)))
        # Fonction de log-vraisemblance
        ll_loglinparam.Mtb <- function(par,Y,I1,I2,I3,I4)
        {
                logmu <- par[1] + I1%*%log(1/(1+exp(par[2:t]))) + I2%*%log(exp(par[2:(t+1)])/(1+exp(par[2:(t+1)]))) + I3%*%log(exp(par[3:(t+1)]+par[t+2])/(1+exp(par[3:(t+1)]+par[t+2]))) + I4%*%log(1/(1+exp(par[3:(t+1)]+par[t+2])))
                loglikelihood <- t(Y)%*%logmu - sum(exp(logmu)) - sum(log(factorial(Y))) 
                loglikelihood
        }        
        # Optimisation de la fonction de log-vraisemblance en fonction des paramètres loglinéaires
        maxim <- optim(par=init_loglinparam.Mtb,fn=ll_loglinparam.Mtb,Y=Y,I1=I1,I2=I2,I3=I3,I4=I4,control=list(fnscale=-1,maxit=1000),method ="L-BFGS-B",hessian=TRUE)
        # Statistiques d'ajustement du modèle
        logmu <- maxim$par[1] + I1%*%log(1/(1+exp(maxim$par[2:t]))) + I2%*%log(exp(maxim$par[2:(t+1)])/(1+exp(maxim$par[2:(t+1)]))) + I3%*%log(exp(maxim$par[3:(t+1)]+maxim$par[t+2])/(1+exp(maxim$par[3:(t+1)]+maxim$par[t+2]))) + I4%*%log(1/(1+exp(maxim$par[3:(t+1)]+maxim$par[t+2])))
        se <- sqrt(exp(2*maxim$par[1])*solve(-maxim$hessian)[1,1]-exp(maxim$par[1]))
        dev <- 2*sum(Y*(pmax(log(Y),0)-logmu))
        df <- 2^t-1-(t+2)
        AIC <- -2*maxim$value+2*(t+2)   
        # Préparation des sorties
        Mtb <- matrix(c(exp(maxim$par[1]),se,dev,df,AIC),nrow=1)
        # Autres paramètres
        parMtb<-matrix(c(exp(maxim$par[1]),exp(maxim$par[2:(t+1)])/(1+exp(maxim$par[2:(t+1)])),exp(maxim$par[3:(t+1)]+maxim$par[t+2])/(1+exp(maxim$par[3:(t+1)]+maxim$par[t+2]))),nrow=1)
        colnames(parMtb) <- c("N",pnames,cnames)

    
        # Préparation des sorties
        tableau <- matrix(Mtb,nrow=1)
        dimnames(tableau) <- list("Mtb",c("abundance","stderr","deviance","df","AIC"))
        ans <- list(n=sum(Y),results=tableau,parMtb=parMtb)
        class(ans) <- "closedp.Mtb"
        ans

}


print.closedp.Mtb <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Abundance estimation and model fit:\n")
        tableau <- x$results
        tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
        tableau[,4] <- round(tableau[,4],0)
        tableau[,c(3,5)] <- round(tableau[,c(3,5)],3)       
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE)
        cat("\nWarning: Model Mtb is not log-linear.\nThe abundance estimation for this model can be unstable.")
        cat("\n\n")
        invisible(x)
}
