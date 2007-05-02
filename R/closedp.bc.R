"closedp.bc" <- function(X,dfreq=FALSE)
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

        mX <- histpos.t(t)
        Y.t <- histfreq.t(X,dfreq=dfreq)
        base.freq <- descriptive(X,dfreq=dfreq)$base.freq
        Y.0 <- base.freq[,1]
        Y.b <- base.freq[,2]
        nbcap <- apply(mX, 1, sum)


        # Afin de supprimer la génération de warnings.
        options(warn=-1)  #Sans cette option, une grande quantité de warnings sont générés car les valeurs de Y ne sont pas entières.
    
    
        # modele M0
        delta<-rep(0,length(Y.0))
        delta[1] <- -1/2
        delta[2] <- 1
        Yd <- pmax(0,Y.0 + delta)
        mX0 <- 1:t
        cons <- log(choose(t, 1:t))
        anaM0 <- glm(Yd~ offset(cons) + mX0,family=poisson)
        NM0 <- sum(Y.0)+exp(anaM0$coef[1]) # calcul de la taille de la population N
        varcovM0 <- summary(anaM0)$cov.unscaled
        erreurtypeM0 <- sqrt(exp(anaM0$coef[1])+(exp(2*anaM0$coef[1]))*varcovM0[1,1])
        M0 <- c(NM0,erreurtypeM0)


        # modele Mt
        delta<-rep(0,length(Y.t))
        delta[nbcap==1] <- (t-2)/(2*t)
        delta[nbcap==2] <- 2/(t*(t-1))
        Yd <- pmax(0,Y.t + delta)
        anaMt <- glm(Yd~mX,family=poisson)
        NMt <- sum(Y.t)+exp(anaMt$coef[1]) # calcul de la taille de la population N
        varcovMt <- summary(anaMt)$cov.unscaled
        erreurtypeMt <- sqrt(exp(anaMt$coef[1])+(exp(2*anaMt$coef[1]))*varcovMt[1,1])
        Mt <- c(NMt,erreurtypeMt)


        # modele Mh Chao
        # Calcul exact :
        NMhC <- sum(Y.0) + ((t - 1) * Y.0[1] * (Y.0[1] - 1))/(2 * t * (Y.0[2] + 1))
        erreurtypeMhC <- sqrt(((t - 1) * Y.0[1] * (Y.0[1] - 1))/(2 * t * (Y.0[2] + 1)) +
        ((t - 1)^2 * Y.0[1] * (Y.0[1] - 1) * (Y.0[1]^2 + 4 * Y.0[1] * Y.0[2] + 3 * Y.0[1] - 6 * Y.0[2] - 6))/(4 * t^2 * (Y.0[2] +1)^2 * (Y.0[2] + 2)))
        MhC <- c(NMhC,erreurtypeMhC)


        # modele Mh Poisson
        if(isTRUE(t>2))
        {
                delta<-rep(0,length(Y.0))
                delta[1] <- -3/4
                delta[2] <- 3/2
                delta[3] <- 1/4
                Yd <- pmax(0,Y.0 + delta)
                mXh <- cbind(1:t,2^(1:t)-1)
                anaMhP <- glm(Yd~ offset(cons) + mXh,family=poisson)
                NMhP <- sum(Y.0)+exp(anaMhP$coef[1]) # calcul de la taille de la population N
                varcovMhP <- summary(anaMhP)$cov.unscaled
                erreurtypeMhP <- sqrt(exp(anaMhP$coef[1])+(exp(2*anaMhP$coef[1]))*varcovMhP[1,1])
                MhP <- c(NMhP,erreurtypeMhP)
        } else {
            MhP <- c(NA,NA)
        }


        # modele Mh Darroch
        delta<-rep(0,length(Y.0))
        delta[1] <- -1
        delta[2] <- 2
        Yd <- pmax(0,Y.0 + delta)
        mXh <- cbind(1:t,((1:t)^2)/2)
        anaMhD <- glm(Yd~ offset(cons) + mXh,family=poisson)
        NMhD <- sum(Y.0)+exp(anaMhD$coef[1]) # calcul de la taille de la population N
        varcovMhD <- summary(anaMhD)$cov.unscaled
        erreurtypeMhD <- sqrt(exp(anaMhD$coef[1])+(exp(2*anaMhD$coef[1]))*varcovMhD[1,1])
        MhD <- c(NMhD,erreurtypeMhD)


        # modele Mth Chao
        if(isTRUE(t>2))
        {
            delta<-rep(0,length(Y.t))
            delta[nbcap==1] <- (t-2)/(2*t)
            delta[nbcap==2] <- 2/(t*(t-1))
            Yd <- pmax(0,Y.t + delta)
            mX4 <- matrix(0, 2^t - 1, t - 2)
            for(i in (3:t)) { mX4[(nbcap == i), (i - 2)] <- 1 }
            anaMthC <- glm(Yd ~ cbind(mX, mX4), family = poisson)
            NMthC <- sum(Y.t)+exp(anaMthC$coef[1]) # calcul de la taille de la population N
            varcovMthC <- summary(anaMthC)$cov.unscaled
            erreurtypeMthC <- sqrt(exp(anaMthC$coef[1])+(exp(2*anaMthC$coef[1]))*varcovMthC[1,1])
            MthC <- c(NMthC,erreurtypeMthC)
        } else {
            MthC <- c(NA,NA)
        }


        # modele Mth Poisson
        if(isTRUE(t>2))
        {
                delta<-rep(0,length(Y.t))
                delta[nbcap==1] <- (2*t-5)/(4*t)
                delta[nbcap==2] <- 3/(t*(t-1))
                delta[nbcap==3] <- 3/(2*t*(t-1)*(t-2))
                Yd <- pmax(0,Y.t + delta)
                mX1 <- 2^nbcap - 1  # variable d interaction
                mX3 <- cbind(mX,mX1) # fusion de ces matrices
                anaMthP <- glm(Yd~mX3,family=poisson)
                NMthP <- sum(Y.t)+exp(anaMthP$coef[1]) # calcul de la taille de la population N
                varcovMthP <- summary(anaMthP)$cov.unscaled
                erreurtypeMthP <- sqrt(exp(anaMthP$coef[1])+(exp(2*anaMthP$coef[1]))*varcovMthP[1,1])
                MthP <- c(NMthP,erreurtypeMthP)
        } else {
            MthP <- c(NA,NA)
        }


        # modele Mth Darroch
        delta<-rep(0,length(Y.t))
        delta[nbcap==1] <- (t-3)/(2*t)
        delta[nbcap==2] <- 4/(t*(t-1))
        Yd <- pmax(0,Y.t + delta)
        mX1 <- (nbcap^2)/2  # variable d interaction
        mX3 <- cbind(mX,mX1) # fusion de ces matrices
        anaMthD <- glm(Yd~mX3,family=poisson)
        NMthD <- sum(Y.t)+exp(anaMthD$coef[1]) # calcul de la taille de la population N
        varcovMthD <- summary(anaMthD)$cov.unscaled
        erreurtypeMthD <- sqrt(exp(anaMthD$coef[1])+(exp(2*anaMthD$coef[1]))*varcovMthD[1,1])
        MthD <- c(NMthD,erreurtypeMthD)


        # modele Mb
        if(isTRUE(t>2))
        {
            delta<-rep(0,length(Y.b))
            delta[1] <- +2
            delta[2] <- -1
            delta[3] <- -1
            Yd <- pmax(0,Y.b + delta)
            mXMb <- c(0:(t-1)) # vecteur de 0 a t-1 , de la taille du nbre d occasions de captures
            anaMb <- glm(Yd~mXMb,family=poisson)
            NMb <- exp(anaMb$coef[1])/(1-exp(anaMb$coef[2])) # calcul de la taille de la population N
            varcovMb <- summary(anaMb)$cov.unscaled 
            v <- NMb*c(1,exp(anaMb$coef[2])/(1-exp(anaMb$coef[2])))    # v est un vecteur cree pour faciliter le calcul de l erreur type
            erreurtypeMb <- sqrt((t(v)%*%varcovMb%*%v) - NMb) # calcul de l'erreur type 
            Mb <- c(NMb,erreurtypeMb)
        } else { 
            Mb <- c(NA,NA)  
        }
        


        # modele Mbh
        if(isTRUE(t>3))
        {
            nX<- X[X[,1]==0,-1]    # matrice des historiques observes avec les suppressions des individus non representatifs
            YMbh <- Y.b[-1] # vecteur du nombre d unites captures pour la premiere fois a l occasion j
            delta<-rep(0,length(YMbh))
            delta[1] <- +2
            delta[2] <- -1
            delta[3] <- -1
            Yd <- pmax(0,YMbh + delta)
            mXMbh <- c(0:(t-2))
            anaMbh <- glm(Yd~mXMbh,family=poisson)
            NMbh <- Y.b[1]+exp(anaMbh$coef[1])/(1-exp(anaMbh$coef[2]))
            varcovMbh <- summary(anaMbh)$cov.unscaled
            v <- (NMbh-Y.b[1]) * c(1, exp(anaMbh$coef[2])/(1 - exp(anaMbh$coef[2])))
            erreurtypeMbh <- sqrt((t(v) %*% varcovMbh %*% v) - NMbh+Y.b[1])
            Mbh <- c(NMbh,erreurtypeMbh)
        } else { 
            Mbh <- c(NA,NA)  
        }


        # Afin de permettre de nouveau la génération de warnings
        options(warn=0)


    # Préparation des sorties
    tableau <- rbind(M0,Mt,MhC,MhP,MhD,MthC,MthP,MthD,Mb,Mbh)
    dimnames(tableau) <- list(c("M0","Mt","Mh Chao","Mh Poisson2","Mh Darroch","Mth Chao","Mth Poisson2","Mth Darroch","Mb","Mbh"),c("abundance","stderr"))
    ans <- list(n=sum(Y.t),results=tableau)
    class(ans) <- "closedp.bc"
    ans
}


print.closedp.bc <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Abundance estimations with bias correction:\n")
        tableau <- round(x$results,1)
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE)
        cat("\n")
        invisible(x)
}
