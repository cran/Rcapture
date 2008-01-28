profileCI <- function(X, dfreq=FALSE, m="M0", h="Chao", a=2, mX=NULL, mname="Customized model", neg=TRUE, alpha=0.05)
{
        #####################################################################################################################################
        # Validation des arguments fournis en entrée
        
        # Argument dfreq
        if(!is.logical(dfreq)||length(dfreq)!=1) stop("'dfreq' must be a logical object of length 1")
    
            X <- as.matrix(X)
            t <- if(dfreq) dim(X)[2]-1 else dim(X)[2]
    
        # Argument X
        if (dfreq)
        {
            if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("Every columns of 'X' but the last one must contain only zeros and ones")
            if (any((X[,t+1]%%1)!=0)) stop("The last column of 'X' must contain capture histories frequencies, therefore integers")
        } else {
            if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
        }
       
        # Argument m
        if(length(m)!=1) stop("'m' must be of length 1")
        if(m!="M0"&&m!="Mt"&&m!="Mh"&&m!="Mth")
        stop("'m' can only take the value 'M0', 'Mt', 'Mh' or 'Mth'")
    
        # Argument h
        if(length(h)!=1) stop("'h' must be of length 1")
        if(!is.function(h)&&h!="Chao"&&h!="Darroch"&&h!="Poisson")
        stop("'h' must be a function or a charater strings taking the value 'Chao', 'Darroch' or 'Poisson'")
        
        # Argument a
        if(length(a)!=1) stop("'a' must be of length 1")
        if (!is.numeric(a)) stop("'a' must be a numeric value")
    
        # Argument mX
        if(!is.null(mX))
        {
            mX<-as.matrix(mX)
            if (dim(mX)[1]!=2^t-1) stop("'mX' must have 2^t-1 rows")
        }

        # Argument neg
        if(!is.logical(neg)||length(neg)!=1) stop("'neg' must be a logical object of length 1")

        #####################################################################################################################################

        
        Y <- histfreq.t(X,dfreq=dfreq)
        mXsans <- if (is.null(mX)) Xclosedp(t,m,h,a)$mat else mX
        rd.call <- match.call()
                      
######### Obtention de l'estimateur de N (sans mu_0). Note : c'est approximatif, on n'ajuste donc pas pour les eta négatifs modèles Chao
        anasans <- glm(Y~mXsans,family=poisson)
        Np <- sum(na.rm=TRUE,Y)+exp(anasans$coef[1])
        varcov <- summary(anasans)$cov.unscaled
        SEp <- sqrt(exp(anasans$coef[1])+(exp(2*anasans$coef[1]))*varcov[1,1])
                               

#-------------------------------------------------------------------------------------------------------------------------------------------
# Idée d'algorithme:
#     - D'abord fixer une borne maximale pour la borne supérieure d'un intervalle de confiance pour mu_0.
#       J'utilise l'intervalle de confiance pour mu_0 chapeau de la régression poisson, construit en postulant la normalité asymptotique.
#     - Ensuite, calculer la log vraisemblance multinommiale profile en 21 points également répartis entre 0 et la borne max inclusivement
#     - Positionner la valeur maximale parmis ces 21 valeurs et calculer la log vraisemblance multinommiale en 9 points répartis également
#       entre mu_0_max et mu_0_max+ et en 9 points répartis également entre mu_0_max- et mu_0_max.
#     - Continuer ce "zoom" jusqu'à ce que mes points soient des entiers consécutifs.
#-------------------------------------------------------------------------------------------------------------------------------------------


######### Détermination des 21 points de départ
        saut <- round((Np-sum(na.rm=TRUE,Y)+2*qnorm(1-alpha/2)*SEp)/20)
        points <- (0:20)*saut
        
######### Boucle de calcul de la vraisemblance multinomiale profile
        continue <- TRUE
        iter <- 1
        while(continue)
        {
            # Calule de la vraisemblance poisson en tous les points
            loglik <- -Inf
            for (n0 in points)
            {
                    Yavec <- c(Y,n0)
                    mXavec <- rbind(mXsans,rep(0,dim(mXsans)[2]))
                    anaavec <- glm(Yavec~mXavec,family=poisson)
                    
                    #############################################################################################################
                    # Vérification eta négatifs pour modèles Chao
                    if(neg)
                    {
                        if ((m=="Mh"||m=="Mth")&&h=="Chao")
                        {
                            # Réajustement du modèle en enlevant les eta négatifs
                            ppositions <- 0
                            param <- anaavec$coef
                            indic <- as.vector(c(rep(0,length(param)-t+2),ifelse(param[(length(param)-t+3):length(param)]<0,1,0)))
                            while(sum(na.rm=TRUE,indic)>0) # Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
                            {
                                # Détermination de la position du premier eta négatif
                                pos <- 1
                                while(indic[pos]==0) pos <- pos + 1
                                ppositions <- c(ppositions,pos)
                                # Retrait de la bonne colonne de mX et réajustement du modèle
                                mXavec <- matrix(mXavec[,-(pos-sum(na.rm=TRUE,ppositions<pos))],nrow=dim(mXavec)[1])
                                anaavec <- glm(Yavec~mXavec,family=poisson)
                                # Ajout de zéros dans le vecteur des paramètres loglinéaires
                                positions <- sort(ppositions[-1])                
                                param <- c(anaavec$coef[1:(positions[1]-1)],0)
                                if(length(positions)>1)
                                {
                                    for ( i in 2:length(positions))
                                    {
                                        if(positions[i]==positions[i-1]+1) {
                                            param <- c(param,0)
                                        } else {
                                            param <- c(param,anaavec$coef[(positions[i-1]-i+2):(positions[i]-i)],0)
                                        }
                                    }
                                }
                                if((positions[length(positions)]-length(positions)+1)<=length(anaavec$coef)) param <- c(param,anaavec$coef[(positions[length(positions)]-length(positions)+1):length(anaavec$coef)])                   
                                indic <- as.vector(c(rep(0,length(param)-t+2),ifelse(param[(length(param)-t+3):length(param)]<0,1,0)))
                            }
                        }
                    }
                    #############################################################################################################
                        
                    # calcule du terme correctif
                    Nn0 <- sum(na.rm=TRUE,Yavec)
                    if(Nn0>100)
                    {
                    ct <- if (n0==0||n0==1) -Nn0+0.5*log(2*pi*Nn0) else n0-Nn0-0.5*log(n0/Nn0)
                    } else { ct <- log((n0^n0)*factorial(Nn0)/((Nn0^Nn0)*factorial(n0))) } 
                    
                    # log vraisemblance multinomiale profile
                    loglik <- c(loglik,(anaavec$deviance - 2*ct)/(-2))                
            }
            loglik <- loglik[-1]
                              
            # Pour former le vecteur complet de log vraisemblance profile
            if(iter==1)
            {
                loglik_c <- loglik
                points_c <- points
            } else {
                loglik_c <- c(loglik_c,loglik)
                points_c <- c(points_c,points)
                loglik_c <- loglik_c[order(points_c)]
                points_c <- points_c[order(points_c)]
            }

            # Objets relatifs à l'itération
            if(saut==1) continue <- FALSE
            iter <- iter + 1
            
            # Détermination du maximum et des points suivants
            lmax <- max(loglik_c)
            pos_lmax <- which.max(loglik_c)
            saut <- ceiling(saut/10)
            points <- c(points_c[pos_lmax]-(9:1)*saut,points_c[pos_lmax]+(1:9)*saut)
            points <- points[points>=0]
            
        }


######### Déterminiation du N qui maximise la log vraisemblance multinomiale profile 
        x <- sum(na.rm=TRUE,Y)+points_c
        N <- x[pos_lmax]
        
######### Interpolation linéaire pour trouver la valeur approximative des bornes de l'intervalle
        sup<-loglik_c>(lmax-qchisq(1-alpha,1)/2)
        posInfMax<-1
        while(!sup[posInfMax]) { posInfMax<-posInfMax+1 }
        posInfMin<-posInfMax-1
        posSupMax<-posInfMax
        while(sup[posSupMax]) { posSupMax<-posSupMax+1 }
        InfCL <- if(posInfMax==1) sum(na.rm=TRUE,Y) else x[posInfMax-1]+((lmax-qchisq(1-alpha,1)/2-loglik_c[posInfMax-1])/(loglik_c[posInfMax]-loglik_c[posInfMax-1]))*(x[posInfMax]-x[posInfMax-1])
        SupCL <- x[posSupMax-1]+(loglik_c[posSupMax-1]-(lmax-qchisq(1-alpha,1)/2))/(loglik_c[posSupMax-1]-loglik_c[posSupMax])*(x[posSupMax]-x[posSupMax-1])
        
######### Production du graphique
        plot(x[posInfMin:posSupMax],loglik_c[posInfMin:posSupMax],type="l",ylab="multinomial profile loglikelihood",xlab="N",main="Profile Likelihood Confidence Interval")
        lInf <- if (InfCL==sum(na.rm=TRUE,Y)) loglik_c[1] else lmax-qchisq(1-alpha,1)/2
        segments(x0=InfCL,y0=min(loglik_c[posInfMin:posSupMax]),x1=InfCL,y1=lInf)
        text(InfCL,min(min(loglik_c[posInfMin:posSupMax])),round(InfCL,2),pos=1,offset=0.2)
        segments(x0=SupCL,y0=min(loglik_c[posInfMin:posSupMax]),x1=SupCL,y1=lmax-qchisq(1-alpha,1)/2)
        text(SupCL,min(loglik_c[posInfMin:posSupMax]),round(SupCL,2),pos=1,offset=0.2)
        segments(x0=N,y0=min(loglik_c[posInfMin:posSupMax]),x1=N,y1=lmax,lty=2)
        text(N,min(loglik_c[posInfMin:posSupMax]),round(N,2),pos=1,offset=0.2)

######### Production sortie non graphique
        results<-matrix(c(N,InfCL,SupCL),nrow=1)
        mname<- if(!is.null(mX)) mname else if (m=="M0"||m=="Mt") m else if (is.function(h)) paste(m,rd.call$h) else if(h=="Poisson") paste(m,paste(h,a,sep="")) else paste(m,h)
        dimnames(results) <- list(mname,c("abundance","InfCL","SupCL"))

        ans <- list(n=sum(na.rm=TRUE,Y),results=results,alpha=alpha)
        class(ans) <- "profileCI"
        ans

}


print.profileCI <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat(paste((1-x$alpha)*100,"%",sep=""),"Profile likelihood confidence interval:\n")
        print.default(x$results, print.gap = 2, quote = FALSE, right=TRUE)
        cat("\n")
        invisible(x)
}
