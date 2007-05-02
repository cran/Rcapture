"histfreq.t" <- function(X,dfreq=FALSE)
{

        X <- as.matrix(X)
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

        nl <- dim(X)[1]
        vecli <- rep(0,nl)      # vecteur du numero de l'historique (selon l'ordre de histpos.t) de chaque individu
        for (i in (1:nl))
        {
                vecli[i] <- 1+sum((1-X[i,1:t])*(2^{(t-1):0}))
        }
        X <- X[vecli<2^t,]      # pour mettre de cote les lignes de X comprenant uniquement des zeros s'il y en a
        vecli <- vecli[vecli<2^t]    
        nl <- length(vecli)
        Y <- rep(0,2^t-1)       # on cree un vecteur de 0 de la taille du nbre des differents historiques de capture possibles
        for (i in (1:nl))
        {
            Y[vecli[i]] <- Y[vecli[i]]+ifelse(dfreq,X[i,t+1],1)
        }
        return(Y)
}
