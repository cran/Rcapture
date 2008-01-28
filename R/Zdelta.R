"Zdelta" <-
function (Xdelta)
{

        Xdelta <- as.matrix(Xdelta)
        I <-dim(Xdelta)[2]
        Z <- matrix(0,dim(Xdelta)[1],2*(I-1))
        Z[,1] <- (1-Xdelta[,1])
        if(I>2)
        {
            i <- 2
            for  (j in (2:(I-1)))
            {
                    Z[,j] <- Z[,(j-1)]*(1-Xdelta[,i])
                    i<- i+1
            }
        }
        Z[,I] <- (1-Xdelta[,I])
        if(I>2)
        {
            i <- 1
            for  (j in ((I+1):(2*(I-1))))
            {
                    Z[,j] <- Z[,(j-1)]*(1-Xdelta[,(I-i)])
                    i <- i+1
            }
        }
        Z
}
