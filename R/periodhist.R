"periodhist" <- function(X,dfreq=FALSE,vt)
{

    #####################################################################################################################################
    # Validation des arguments fournis en entrée
    
    # Argument dfreq
    if(!is.logical(dfreq)||length(dfreq)!=1) stop("'dfreq' must be a logical object of length 1")

        X <- as.matrix(X)
        t <- if(dfreq) dim(X)[2]-1 else dim(X)[2]
    
    # Argument vt
    if(t!=sum(na.rm=TRUE,vt)) stop("The number of columns in 'X' is not equal to the total number of capture occasions (sum of the 'vt' components)")
    if (any((vt %% 1)!=0)) stop("The 'vt' components must be integers")

    # Argument X
    if (dfreq)
    {
        if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("Every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,t+1] %% 1)!=0)) stop("The last column of 'X' must contain capture histories frequencies, therefore integers")
    } else {
        if (any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
       
    #####################################################################################################################################

        I <- length(vt)
        Xinter <- if(dfreq) matrix(0,dim(X)[1],I+1) else matrix(0,dim(X)[1],I)
        for (i in (1:I))
        {
                if (i==1) { Xs <- matrix(X[,c(1:vt[i])],nrow=dim(X)[1]) } else
                { Xs <- matrix(X[,c((sum(na.rm=TRUE,vt[1:(i-1)])+1):sum(na.rm=TRUE,vt[1:i]))],nrow=dim(X)[1]) }
                Xinter[,i] <- apply(Xs,1,max)
        }
        if (dfreq) Xinter[,I+1] <- X[,sum(na.rm=TRUE,vt)+1]
        Y <- histfreq.t(Xinter,dfreq=dfreq)
        Xinterfreq <- cbind(histpos.t(I),Y)
        return(Xinterfreq)
}
