"histpos.t" <- function(t)
{
    #####################################################################################################################################
    # Validation de l'argument fourni en entr�e
    if (length(t)!=1||any((t %% 1)!=0)) stop("'t' must be an integer of length 1")
    #####################################################################################################################################

    X0 <- cbind(c(1, 1, 0), c(1, 0, 1))
    if (t==2) { # il y n'y a que 2 occasions de capture 
        X0 
    } else {   # il y a plus de 2 occasions de capture
        for(i in (3:t)) {
              X1 <- cbind( c(rep(1, 2^(i-1)), rep(0,((2^(i-1))-1))), rbind(X0,rep(0,(i-1)),X0))
              X0 <- X1
        }
        X1
    }
}


"histpos.0" <- function(vt)
{
    #####################################################################################################################################
    # Validation de l'argument fourni en entr�e
    if (any((vt %% 1)!=0)) stop("the 'vt' components must be integers")
    #####################################################################################################################################

    I<-length(vt)
    if (I==1) {
         mat <- matrix(vt:1,ncol=1)
    } else {
         vt<-vt+1
         vec <- (vt[I] - 1):0
         mat<-cbind(rep(((vt[I-1]-1):0),rep(length(vec),vt[I-1])),rep(vec,vt[I-1]))
         if (I>2) {
              for (i in (2:(I-1))) {
                   nl<-dim(mat)[1]
                   mvec<-rep((vt[I-i]-1):0,rep(nl,vt[I-i]))
                   for (j in (1:i)) {
                           mvec<-cbind(mvec,rep(mat[,j],vt[I-i]))
                   }
                   mat<-mvec
              }     
         }
         mat <- mat[-dim(mat)[1],]
         colnames(mat)<-NULL
    }
    return(mat)    
}