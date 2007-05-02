"histfreq.0" <-
function(X,dfreq=FALSE,vt)
{

    #####################################################################################################################################
    # Validation des arguments fournis en entrée
    
    # Argument dfreq
    if(!is.logical(dfreq)||!isTRUE(all.equal(length(dfreq),1))) stop("'dfreq' must be a logical object of length 1")
    
    #Argument X
    X <- as.matrix(X)
    if (dfreq)
    {
        if (any(X[,1:sum(vt)]!=1&X[,1:sum(vt)]!=0)) stop("Every columns of 'X' but the last one must contain only zeros and ones")
        if (any((X[,sum(vt)+1] %% 1)!=0)) stop("The last column of 'X' must contain capture histories frequencies, therefore integers")
    } else {
        if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
    }
    
    # Argument vt
    if(!isTRUE(all.equal(ifelse(dfreq,dim(X)[2]-1,dim(X)[2]),sum(vt))))
        stop("The number of columns in 'X' is not equal to the total number of capture occasions (sum of the 'vt' components)")
    if (any((vt %% 1)!=0)) stop("The 'vt' components must be integers")
    
    #####################################################################################################################################

        I <- length(vt) # nombre de periodes primaires
        Xh<-apply(X[,1:vt[1]],1,sum)
        for (i in 2:I)
        {
                Xh<-cbind(Xh,apply(X[,c((sum(vt[1:(i-1)])+1):sum(vt[1:i]))],1,sum))
        }
        nl<-dim(Xh)[1]
        vecli<-rep(0,nl)    # vecteur du numero de l'historique (selon l'ordre de histpos.0) de chaque individu
        cvt<-c(cumprod((vt+1)[I:2])[(I-1):1],1)
        for (i in (1:nl))
        {
                vecli[i]<-1+sum((vt-Xh[i,])*cvt)
        }
        X <- X[vecli<prod(vt+1),]      # pour mettre de cote les lignes de X comprenant uniquement des zeros s'il y en a
        vecli<-vecli[vecli<prod(vt+1)]
        nl<-length(vecli)
        Y <- rep(0,prod(vt+1)-1)       # on cree un vecteur de 0 de la taille du nbre  des differents historiques de capture possibles
        for (i in (1:nl))
        {
                Y[vecli[i]]<-Y[vecli[i]]+ifelse(dfreq,X[i,sum(vt)+1],1)
        }
        return(Y)
}
