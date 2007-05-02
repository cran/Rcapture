"histpos.0" <- function(vt)
{

    #####################################################################################################################################
    # Validation de l'argument fourni en entrée
    if (any((vt %% 1)!=0)) stop("The 'vt' components must be integers")
    #####################################################################################################################################

        vt<-vt+1
        I<-length(vt)
        vec <- (vt[I] - 1):0
        mat<-cbind(rep(((vt[I-1]-1):0),rep(length(vec),vt[I-1])),rep(vec,vt[I-1]))
        if (isTRUE(I>2))
        {
            for (i in (2:(I-1)))
            {
                    nl<-dim(mat)[1]
                    mvec<-rep((vt[I-i]-1):0,rep(nl,vt[I-i]))
                    for (j in (1:i))
                    {
                            mvec<-cbind(mvec,rep(mat[,j],vt[I-i]))
                    }
                    mat<-mvec
            }
        }
        colnames(mat)<-NULL
        mat[-dim(mat)[1],]
}
