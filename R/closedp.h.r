"closedp.h" <- function(X, dfreq=FALSE, m="Mh", h="Poisson", a=2)
{

    #####################################################################################################################################
    # Validation des arguments fournis en entrée
    
    # Argument dfreq
    if(!is.logical(dfreq)||length(dfreq)!=1) stop("'dfreq' must be a logical object of length 1")

        X <- as.matrix(X)
        t <- if(dfreq) dim(X)[2]-1 else dim(X)[2]

    # Argument X
    X <- as.matrix(X)
    if (dfreq)
    {
        if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("Every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,t+1]%%1)!=0)) stop("The last column of 'X' must contain capture histories frequencies, therefore integers")
    } else {
        if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
    
    # Argument m
    if(length(m)!=1) stop("'m' must be of length 1")
    if(m!="Mh"&&m!="Mth") stop("'m' can only take the values 'Mh' and 'Mth'")
    
    # Argument h
    if(length(h)!=1) stop("'h' must be of length 1")
    if(!is.function(h)&&h!="Poisson") stop("'h' must be a function or a charater strings taking the value 'Poisson'")
    
    # Argument a
    if(length(a)!=1) stop("'a' must be of length 1")
    if (!is.numeric(a)) stop("'a' must be a numeric value")
    
    #####################################################################################################################################

        histpos <- histpos.t(t)
        Y <- histfreq.t(X,dfreq=dfreq)
        nbcap <- apply(histpos, 1, sum)

        mX2 <- if (is.function(h)) h(nbcap) else if (h=="Poisson") a^nbcap - 1 
        mX <- if (m=="Mh") cbind(nbcap,mX2) else if(m=="Mth") cbind(histpos,mX2)
        anaM <- glm(Y~mX,family=poisson)
        NM <- sum(na.rm=TRUE,Y)+exp(anaM$coef[1]) # calcul de la taille de la population N
        varcovM <- summary(anaM)$cov.unscaled
        erreurtypeM <- sqrt(exp(anaM$coef[1])+(exp(2*anaM$coef[1]))*varcovM[1,1])
        M <- matrix(c(NM,erreurtypeM,anaM$dev,anaM$df.residual,anaM$aic),nrow=1)


        # Préparation des sorties
        closedp.call<-match.call()
        modelname <- if (is.function(h)) paste(m,closedp.call$h) else if(h=="Poisson") paste(m,paste(h,a,sep=""))
        dimnames(M) <- list(modelname,c("abundance","stderr","deviance","df","AIC"))
        ans <- list(n=sum(na.rm=TRUE,Y),results=M,glm=anaM)
        class(ans) <- "closedp.custom"
        ans

}
