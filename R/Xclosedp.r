"Xclosedp" <-
function(t,m,h,a)
{

        histpos <- histpos.t(t)

        if (m=="none") # no model
        {
            mXp <- as.matrix(apply(histpos,1,max))
        } else if (m=="M0") # modele M0
                {
                    mXp <- as.matrix(apply(histpos, 1, sum))
                } else if (m=="Mt") # modele Mt
                        {
                            mXp <- histpos
                        } else if (m=="Mh") # modele Mh
                                {
                                    mXp1 <- apply(histpos,1,sum) 
                                    if (h=="Chao")
                                    {
                                        mXp2 <- matrix(0,2^t-1,t-2)
                                        for (j in (3:t)) { mXp2[,j-2]<-pmax(mXp1-j+1,0) }
                                    } else
                                    if (h=="Poisson") mXp2 <- a^mXp1 - 1 else
                                    if (h=="Darroch") mXp2 <- (mXp1^2)/2 else
                                    if(is.function(h)) mXp2 <- h(mXp1)
                                    mXp <- cbind(mXp1,mXp2)
                                } else if (m=="Mth") # modele Mth
                                        {
                                            mXp1 <- apply(histpos,1,sum)
                                            if (h=="Chao")
                                            {
                                                mXp2 <- matrix(0,2^t-1,t-2)
                                                for (j in (3:t)) { mXp2[,j-2]<-pmax(mXp1-j+1,0) }
                                            } else
                                            if (h=="Poisson") mXp2 <- a^mXp1 - 1 else
                                            if (h=="Darroch") mXp2 <- (mXp1^2)/2 else
                                            if(is.function(h)) mXp2 <- h(mXp1)
                                            mXp <- cbind(histpos,mXp2)
                                        }
                            
          list(mat=mXp,nbparam=dim(mXp)[2])
}
