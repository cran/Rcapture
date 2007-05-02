"closedp" <- function(X,dfreq=FALSE,neg=TRUE)
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
    
    # Argument neg
    if(!is.logical(neg)||!isTRUE(all.equal(length(neg),1))) stop("'neg' must be a logical object of length 1")
    
    #####################################################################################################################################

        Y <- histfreq.t(X,dfreq=dfreq)
        
        # Création de matrices qui seront utiles pour construire les différentes matrices X des modèles
        histpos <- histpos.t(t)
        nbcap <- apply(histpos, 1, sum)
        nbcap_av <- rep(0,2^(t-1))  # nbre d'occasions de capture avant la premiere capture
        for(i in (t-1):1) { nbcap_av<-c(nbcap_av,rep(t-i,2^(i-1))) }
        nbcap_ap <- (nbcap-1) # nbre de capture apres la premiere capture
        inv_c1<-1-histpos[,1]
        if (t>2)
        {
            mXchao <- matrix(0,2^t-1,t-2)
            for (i in 3:t) { mXchao[,i-2]<-pmax(nbcap-i+1,0) }
        }
        mXP <- 2^nbcap-1
        mXD <- nbcap^2/2
        
        # Création de vecteurs de noms qui permettent de rendre plus clair les objects glm fournis en sortie
        betanames <- vector("character",t)
        for (i in 1:t) { betanames[i]<-paste("beta",i,sep="") }
        etanames <- vector("character",t-2)
        for (i in 3:t) { etanames[i-2]<-paste("eta",i,sep="") }
        
        pnames <- vector("character",t)
        for (i in 1:t) { pnames[i]<-paste("p",i,sep="") }
        cnames <- vector("character",t-1)
        for (i in 2:t) { cnames[i-1]<-paste("c",i,sep="") }

        # Identifiant du moment de la première capture pour chaque historique : intervient dans le calcul des p
        ifirstcap <- 0
        for (i in 1:t) { ifirstcap <- c(ifirstcap,rep(i,2^(t-i))) }
        ifirstcap <- ifirstcap[-1]


        # modele M0
        mXM0. <- matrix(nbcap,ncol=1)
        colnames(mXM0.) <- "beta"
        anaM0 <- glm(Y~mXM0.,family=poisson)
        NM0 <- sum(Y)+exp(anaM0$coef[1]) # calcul de la taille de la population N
        varcovM0 <- summary(anaM0)$cov.unscaled
        erreurtypeM0 <- sqrt(exp(anaM0$coef[1])+(exp(2*anaM0$coef[1]))*varcovM0[1,1])
        M0 <- c(NM0,erreurtypeM0,anaM0$dev,anaM0$df.residual,anaM0$aic)
        # Autres paramètres
        parM0<-matrix(c(NM0,exp(anaM0$coef[2])/(1+exp(anaM0$coef[2]))),nrow=1)
        colnames(parM0) <- c("N","p")


        # modele Mt
        mXMt. <- histpos
        colnames(mXMt.) <- betanames
        anaMt <- glm(Y~mXMt.,family=poisson)
        NMt <- sum(Y)+exp(anaMt$coef[1]) # calcul de la taille de la population N
        varcovMt <- summary(anaMt)$cov.unscaled
        erreurtypeMt <- sqrt(exp(anaMt$coef[1])+(exp(2*anaMt$coef[1]))*varcovMt[1,1])
        Mt <- c(NMt,erreurtypeMt,anaMt$dev,anaMt$df.residual,anaMt$aic)
        # Autres paramètres
        parMt<-matrix(c(NMt,exp(anaMt$coef[2:(t+1)])/(1+exp(anaMt$coef[2:(t+1)]))),nrow=1)
        colnames(parMt) <- c("N",pnames)


        # modele Mh Chao
        if (t>2)
        {
            mXMhC. <- cbind(nbcap, mXchao)
            colnames(mXMhC.) <- c("beta",etanames)
            anaMhC <- glm(Y ~ mXMhC., family = poisson)
            # Réajustement du modèle en enlevant les eta négatifs
            ppositions <- 0
            if(neg)
            {
                param <- anaMhC$coef
                indic <- as.vector(c(0,0,ifelse(param[-(1:2)]<0,1,0)))
                while(isTRUE(sum(indic)>0)) # Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
                {
                    # Détermination de la position du premier eta négatif
                    pos <- 1
                    while(isTRUE(all.equal(indic[pos],0))) pos <- pos + 1
                    ppositions <- c(ppositions,pos)
                    # Retrait de la bonne colonne de mX et réajustement du modèle
                    mXMhC. <- mXMhC.[,-(pos-sum(ppositions<pos))]
                    anaMhC <- glm(Y~mXMhC.,family=poisson)
                    # Ajout de zéros dans le vecteur des paramètres loglinéaires
                    positions <- sort(ppositions[-1])                
                    param <- c(anaMhC$coef[1:(positions[1]-1)],0)
                    if(isTRUE(length(positions)>1))
                    {
                        for ( i in 2:length(positions))
                        {
                            if(isTRUE(all.equal(positions[i],positions[i-1]+1))) {
                                param <- c(param,0)
                            } else {
                                param <- c(param,anaMhC$coef[(positions[i-1]-i+2):(positions[i]-i)],0)
                            }
                        }
                    }
                    param <- c(param,anaMhC$coef[(positions[length(positions)]-length(positions)+1):length(anaMhC$coef)])
                    indic <- as.vector(c(0,0,ifelse(param[-(1:2)]<0,1,0)))
                }
            }
            posMhC <- sort(ppositions[-1])
            # Estimation de l'abondance
            NMhC <- sum(Y)+exp(anaMhC$coef[1]) # calcul de la taille de la population N
            varcovMhC <- summary(anaMhC)$cov.unscaled
            erreurtypeMhC <- sqrt(exp(anaMhC$coef[1])+(exp(2*anaMhC$coef[1]))*varcovMhC[1,1])
            MhC <- c(NMhC,erreurtypeMhC,anaMhC$dev,anaMhC$df.residual,anaMhC$aic)
            # Autres paramètres
            pMhC <- sum(anaMhC$fitted.values[ifirstcap==1])/NMhC
            parMhC<-matrix(c(NMhC,pMhC),nrow=1)
            colnames(parMhC) <- c("N","p")
        } else {
            anaMhC <- NULL
            MhC <- rep(NA,5)
            posMhC <- NULL
            parMhC <- rep(NA,2)
        }


        # modele Mh Poisson
        mXMhP. <- cbind(nbcap,mXP)
        colnames(mXMhP.) <- c("beta","tau")
        anaMhP <- glm(Y~mXMhP.,family=poisson)
        NMhP <- sum(Y)+exp(anaMhP$coef[1]) # calcul de la taille de la population N
        varcovMhP <- summary(anaMhP)$cov.unscaled
        erreurtypeMhP <- sqrt(exp(anaMhP$coef[1])+(exp(2*anaMhP$coef[1]))*varcovMhP[1,1])
        MhP <- c(NMhP,erreurtypeMhP,anaMhP$dev,anaMhP$df.residual,anaMhP$aic)
        # Autres paramètres
        pMhP <- sum(anaMhP$fitted.values[ifirstcap==1])/NMhP
        parMhP<-matrix(c(NMhP,pMhP),nrow=1)
        colnames(parMhP) <- c("N","p")


        # modele Mh Darroch
        mXMhD. <- cbind(nbcap,mXD)
        colnames(mXMhD.) <- c("beta","tau")
        anaMhD <- glm(Y~mXMhD.,family=poisson)
        NMhD <- sum(Y)+exp(anaMhD$coef[1]) # calcul de la taille de la population N
        varcovMhD <- summary(anaMhD)$cov.unscaled
        erreurtypeMhD <- sqrt(exp(anaMhD$coef[1])+(exp(2*anaMhD$coef[1]))*varcovMhD[1,1])
        MhD <- c(NMhD,erreurtypeMhD,anaMhD$dev,anaMhD$df.residual,anaMhD$aic)
        # Autres paramètres
        pMhD <- sum(anaMhD$fitted.values[ifirstcap==1])/NMhD
        parMhD<-matrix(c(NMhD,pMhD),nrow=1)
        colnames(parMhD) <- c("N","p")


        # modele Mth Chao
        if (t>2)
        {
            mXMthC. <- cbind(histpos, mXchao)
            colnames(mXMthC.) <- c(betanames,etanames)
            anaMthC <- glm(Y ~ mXMthC., family = poisson)
            # Réajustement du modèle en enlevant les eta négatifs
            ppositions <- 0
            if(neg)
            {
                param <- anaMthC$coef
                indic <- as.vector(c(rep(0,t+1),ifelse(param[-(1:(t+1))]<0,1,0)))
                while(isTRUE(sum(indic)>0)) # Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
                {
                    # Détermination de la position du premier eta négatif
                    pos <- 1
                    while(isTRUE(all.equal(indic[pos],0))) pos <- pos + 1
                    ppositions <- c(ppositions,pos)
                    # Retrait de la bonne colonne de mX et réajustement du modèle
                    mXMthC. <- mXMthC.[,-(pos-sum(ppositions<pos))]
                    anaMthC <- glm(Y~mXMthC.,family=poisson)
                    # Ajout de zéros dans le vecteur des paramètres loglinéaires
                    positions <- sort(ppositions[-1])                
                    param <- c(anaMthC$coef[1:(positions[1]-1)],0)
                    if(isTRUE(length(positions)>1))
                    {
                        for ( i in 2:length(positions))
                        {
                            if(isTRUE(all.equal(positions[i],positions[i-1]+1))) {
                                param <- c(param,0)
                            } else {
                                param <- c(param,anaMthC$coef[(positions[i-1]-i+2):(positions[i]-i)],0)
                            }
                        }
                    }
                    param <- c(param,anaMthC$coef[(positions[length(positions)]-length(positions)+1):length(anaMthC$coef)])
                    indic <- as.vector(c(rep(0,t+1),ifelse(param[-(1:(t+1))]<0,1,0)))
                }
            } 
            posMthC <- sort(ppositions[-1])
            # Estimation de l'abondance
            NMthC <- sum(Y)+exp(anaMthC$coef[1])
            varcovMthC <- summary(anaMthC)$cov.unscaled
            erreurtypeMthC <- sqrt(exp(anaMthC$coef[1])+(exp(2*anaMthC$coef[1]))*varcovMthC[1,1])
            MthC <- c(NMthC,erreurtypeMthC,anaMthC$dev,anaMthC$df.residual,anaMthC$aic)
            # Autres paramètres
            upredMthC <- rep(0,t)
            for ( i in 1:t ) { upredMthC[i] <- sum(anaMhC$fitted.values[ifirstcap==i]) }
            denoPMthC <- NMthC
            for ( i in 2:t ) { denoPMthC <- c(denoPMthC,NMthC-sum(upredMthC[1:(i-1)])) }          
            parMthC<-matrix(c(NMthC,upredMthC/denoPMthC),nrow=1)
            colnames(parMthC) <- c("N",pnames)
        } else {
            anaMthC <- NULL
            MthC <- rep(NA,5)
            posMthC <- NULL
            parMthC <- rep(NA,2)
        }


        # modele Mth Poisson
        mXMthP. <- cbind(histpos,mXP)
        colnames(mXMthP.) <- c(betanames,"tau")
        anaMthP <- glm(Y~mXMthP.,family=poisson)
        NMthP <- sum(Y)+exp(anaMthP$coef[1]) # calcul de la taille de la population N
        varcovMthP <- summary(anaMthP)$cov.unscaled
        erreurtypeMthP <- sqrt(exp(anaMthP$coef[1])+(exp(2*anaMthP$coef[1]))*varcovMthP[1,1])
        MthP <- c(NMthP,erreurtypeMthP,anaMthP$dev,anaMthP$df.residual,anaMthP$aic)
        # Autres paramètres
        upredMthP <- rep(0,t)
        for ( i in 1:t ) { upredMthP[i] <- sum(anaMhP$fitted.values[ifirstcap==i]) }
        denoPMthP <- NMthP
        for ( i in 2:t ) { denoPMthP <- c(denoPMthP,NMthP-sum(upredMthP[1:(i-1)])) }          
        parMthP<-matrix(c(NMthP,upredMthP/denoPMthP),nrow=1)
        colnames(parMthP) <- c("N",pnames)


        # modele Mth Darroch
        mXMthD. <- cbind(histpos,mXD)
        colnames(mXMthD.) <- c(betanames,"tau")
        anaMthD <- glm(Y~mXMthD.,family=poisson)
        NMthD <- sum(Y)+exp(anaMthD$coef[1]) # calcul de la taille de la population N
        varcovMthD <- summary(anaMthD)$cov.unscaled
        erreurtypeMthD <- sqrt(exp(anaMthD$coef[1])+(exp(2*anaMthD$coef[1]))*varcovMthD[1,1])
        MthD <- c(NMthD,erreurtypeMthD,anaMthD$dev,anaMthD$df.residual,anaMthD$aic)
        # Autres paramètres
        upredMthD <- rep(0,t)
        for ( i in 1:t ) { upredMthD[i] <- sum(anaMhD$fitted.values[ifirstcap==i]) }
        denoPMthD <- NMthD
        for ( i in 2:t ) { denoPMthD <- c(denoPMthD,NMthD-sum(upredMthD[1:(i-1)])) }          
        parMthD<-matrix(c(NMthD,upredMthD/denoPMthD),nrow=1)
        colnames(parMthD) <- c("N",pnames)


        # modele Mb
        mXMb. <- cbind(nbcap_av,nbcap_ap)
        colnames(mXMb.) <- c("beta1","beta2")
        anaMb <- glm(Y~mXMb.,family=poisson)
        NMb <-(exp(anaMb$coef[1])*(1+exp(anaMb$coef[3]))^t)/(1+exp(anaMb$coef[3])-exp(anaMb$coef[2])) # calcul de la taille de la population N
        varcovMb <- summary(anaMb)$cov.unscaled 
        v1<-1 
        v2<-exp(anaMb$coef[2])/(1+exp(anaMb$coef[3])-exp(anaMb$coef[2])) 
        v3 <-(t*exp(anaMb$coef[3])/(1+exp(anaMb$coef[3])) -exp(anaMb$coef[3])/(1+exp(anaMb$coef[3])-exp(anaMb$coef[2]))) 
        v <- NMb*c(v1,v2,v3)
        erreurtypeMb <-sqrt((t(v)%*%varcovMb%*%v)-NMb) # calcul de l erreur type 
        Mb <- c(NMb,erreurtypeMb,anaMb$dev,anaMb$df.residual,anaMb$aic)
        # Autres paramètres
        parMb<-matrix(c(NMb,1-exp(anaMb$coef[2])/(1+exp(anaMb$coef[3])),exp(anaMb$coef[3])/(1+exp(anaMb$coef[3]))),nrow=1)
        colnames(parMb) <- c("N","p","c")


        # modele Mbh
        if (t>2)
        {
            mXMbh. <- cbind(inv_c1,nbcap_av,nbcap_ap)
            colnames(mXMbh.) <- c("eta","beta1","beta2")
            anaMbh <- glm(Y~mXMbh.,family=poisson)
            NMbh <- exp(anaMbh$coef[1])*((1+exp(anaMbh$coef[4]))^(t-1))*(1+exp(anaMbh$coef[2]+anaMbh$coef[3])/(1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3])))  # calcul de la taille de la population N
            varcovMbh <- summary(anaMbh)$cov.unscaled
            v1 <- (1+exp(anaMbh$coef[2]+anaMbh$coef[3])/(1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3])))
            v2 <- exp(anaMbh$coef[2]+anaMbh$coef[3])/(1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3]))
            v3 <- exp(anaMbh$coef[2]+anaMbh$coef[3])/(1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3])) + exp(anaMbh$coef[2]+2*anaMbh$coef[3])/((1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3]))^2)
            v4 <- (t-1)*exp(anaMbh$coef[4])/(1+exp(anaMbh$coef[4])) - exp(anaMbh$coef[2]+anaMbh$coef[3]+anaMbh$coef[4])/((1+exp(anaMbh$coef[4])-exp(anaMbh$coef[3]))^2)
            v <- exp(anaMbh$coef[1])*((1+exp(anaMbh$coef[4]))^(t-1))*c(v1,v2,v3,v4)
            erreurtypeMbh <- sqrt((t(v)%*%varcovMbh%*%v) - NMbh)  # calcul de l erreur type
            Mbh <- c(NMbh,erreurtypeMbh,anaMbh$dev,anaMbh$df.residual,anaMbh$aic)
            # Autres paramètres
            pMbh <- sum(anaMbh$fitted.values[ifirstcap==1])/NMbh
            cMbh <- (anaMbh$fitted.values[1]/(NMbh*pMbh))^(1/(t-1))
            parMbh<-matrix(c(NMbh,pMbh,cMbh),nrow=1)
            colnames(parMbh) <- c("N","p","c")
        } else {
            anaMbh <- NULL
            Mbh <- rep(NA,5)
            parMbh <- rep(NA,2)
        }

    
        # Préparation des sorties
        tableau <- rbind(M0,Mt,MhC,MhP,MhD,MthC,MthP,MthD,Mb,Mbh)
        dimnames(tableau) <- list(c("M0","Mt","Mh Chao","Mh Poisson2","Mh Darroch","Mth Chao","Mth Poisson2","Mth Darroch","Mb","Mbh"),c("abundance","stderr","deviance","df","AIC"))
        ans <- list(n=sum(Y),t=t,results=tableau,glmM0=anaM0,glmMt=anaMt,glmMhC=anaMhC,glmMhP=anaMhP,glmMhD=anaMhD,glmMthC=anaMthC,glmMthP=anaMthP,glmMthD=anaMthD,glmMb=anaMb,glmMbh=anaMbh,
                                                 parM0=parM0,parMt=parMt,parMhC=parMhC,parMhP=parMhP,parMhD=parMhD,parMthC=parMthC,parMthP=parMthP,parMthD=parMthD,parMb=parMb,parMbh=parMbh,negMhC=posMhC,negMthC=posMthC)
        class(ans) <- "closedp"
        ans

}


print.closedp <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Abundance estimations and model fits:\n")
        tableau <- x$results
        tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
        tableau[,4] <- round(tableau[,4],0)
        tableau[,c(3,5)] <- round(tableau[,c(3,5)],3)       
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE)
        if (isTRUE(all.equal(length(x$negMhC),1))) cat("\nNote:",length(x$negMhC),"eta parameter has been set to zero in the Mh Chao model")
        if (isTRUE(length(x$negMhC)>1)) cat("\nNote:",length(x$negMhC),"eta parameters has been set to zero in the Mh Chao model")
        if (isTRUE(all.equal(length(x$negMthC),1))) cat("\nNote:",length(x$negMthC),"eta parameter has been set to zero in the Mth Chao model")
        if (isTRUE(length(x$negMthC)>1)) cat("\nNote:",length(x$negMthC),"eta parameters has been set to zero in the Mth Chao model")
        cat("\n\n")
        invisible(x)
}


boxplot.closedp <- function(x, ...){
        boxplot.default((x$glmM0$y-fitted(x$glmM0))/sqrt(fitted(x$glmM0)),(x$glmMt$y-fitted(x$glmMt))/sqrt(fitted(x$glmMt)),(x$glmMhC$y-fitted(x$glmMhC))/sqrt(fitted(x$glmMhC)),
                (x$glmMhP$y-fitted(x$glmMhP))/sqrt(fitted(x$glmMhP)),(x$glmMhD$y-fitted(x$glmMhD))/sqrt(fitted(x$glmMhD)),(x$glmMthC$y-fitted(x$glmMthC))/sqrt(fitted(x$glmMthC)),
                (x$glmMthP$y-fitted(x$glmMthP))/sqrt(fitted(x$glmMthP)),(x$glmMthD$y-fitted(x$glmMthD))/sqrt(fitted(x$glmMthD)),(x$glmMb$y-fitted(x$glmMb))/sqrt(fitted(x$glmMb)),
                (x$glmMbh$y-fitted(x$glmMbh))/sqrt(fitted(x$glmMbh)),names=c("M0","Mt","MhC","MhP","MhD","MthC","MthP","MthD","Mb","Mbh"),main="Boxplots of Pearson Residuals")
				abline(h=0,lty=3)
}
