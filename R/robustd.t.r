"robustd.t" <- function(X, dfreq=FALSE, vt, vm="M0", vh=list("Chao"), va=2, neg=TRUE)
{

        I <- length(vt) # nombre de periodes primaires

############################################################################################################################################
# Validation des arguments fournis en entr�e et changement de leur forme si n�cessaire
############################################################################################################################################

        # Argument dfreq
        if(!is.logical(dfreq)||!isTRUE(all.equal(length(dfreq),1))) stop("'dfreq' must be a logical object of length 1")
    
        # Argument X
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
    
        # Argument vm
        if(!(isTRUE(all.equal(length(vt),length(vm)))||isTRUE(all.equal(length(vm),1)))) stop("'vm' must be of length 1 or have the same length than 'vt'")
        for (i in 1:length(vm))
        {
            if(!identical(vm[i],"none")&&!identical(vm[i],"M0")&&!identical(vm[i],"Mt")&&!identical(vm[i],"Mh")&&!identical(vm[i],"Mth"))
            stop("'vm' components can only take the value 'none', 'M0', 'Mt', 'Mh' or 'Mth'")
            if((identical(vm[i],"Mh")||identical(vm[i],"Mth"))&&isTRUE(vt[i]<3))
            stop("Mh and Mth models require at least 3 capture occasions")
        }
        if(identical(vm[1],"none")||identical(vm[length(vm)],"none")) stop("The 'no model' cannot be chosen for the first or the last period")
        # Forme
        if (isTRUE(all.equal(length(vm),1))) vm <- rep(vm,length(vt))       
    
        # Argument vh
        nh <- length(vm[vm=="Mh"|vm=="Mth"])
        if (isTRUE(nh>0))
        {
            vh <- as.list(vh)
            if(!(isTRUE(all.equal(length(vh),nh))||isTRUE(all.equal(length(vh),1)))) stop("'vh' must be of length 1 or of length equal to the number of heterogenous models specified in 'vm'")
            for (i in 1:length(vh))
            {
                if(!is.function(vh[[i]])&&!identical(vh[[i]],"Chao")&&!identical(vh[[i]],"Darroch")&&!identical(vh[[i]],"Poisson"))
                stop("The 'vh' elements must be functions or charater strings taking the value 'Chao', 'Darroch' or 'Poisson'")
            }
            # Forme
            if (isTRUE(all.equal(length(vh),1))) vh <- rep(vh,nh)
            vht<-rep(NA,I)
            j<-1
            for (i in 1:I)
            { 
                if (identical(vm[i],"Mh")||identical(vm[i],"Mth"))
                {
                    vht[[i]] <- vh[[j]]
                    j <- j+1
                }
            }
            vh <- vht
        }
        
        # Argument va
        nP <- length(vh[vh=="Poisson"])
        if (isTRUE(nP>0))
        {
            if(!(isTRUE(all.equal(length(va),nP))||isTRUE(all.equal(length(va),1)))) stop("'va' must be of length 1 or of length to equals the number of Poisson heterogenous models specified in 'vh'")
            if (!is.numeric(va)) stop("The 'va' components must be numeric values")
            # Forme
            if (isTRUE(all.equal(length(va),1))) va <- rep(va,nP)
            vat<-rep(1,I)
            j<-1
            for (i in 1:I)
            { 
                if (identical(vh[i],"Poisson"))
                {
                    vat[i] <- va[j]
                    j <- j+1
                }
            }
            va <- vat
        }
    
        # Argument neg
        if(!is.logical(neg)||!isTRUE(all.equal(length(neg),1))) stop("'neg' must be a logical object of length 1")



############################################################################################################################################
# AJUSTEMENT DU MOD�LE
############################################################################################################################################


#-------------------------------------#
# �laboration et ajustement du mod�le #
#-------------------------------------#

        Y <- histfreq.t(X,dfreq=dfreq)
        rd.call <- match.call()
        Xw <- Xomega.t(vt,vm,vh,va,rd.call)     # deuxieme composante (celle de Beta) dans le modele loglineaire du robust design
        Xp <- histpos.t(sum(vt))
        Xdelta <- matrix(0,dim(Xp)[1],I)
        for (i in 1:I)
        {
                if (isTRUE(all.equal(i,1))) { Xs <- Xp[,c(1:vt[i])] } else
                { Xs <- Xp[,c((sum(vt[1:(i-1)])+1):sum(vt[1:i]))] }
                Xdelta[,i] <- apply(Xs,1,max)
        }
        Zw <- Zdelta(Xdelta)     # premiere composante (celle de Alpha) dans ce meme modele
        gammanames <- vector("character",2*I-2)
        for (i in 1:(2*I-2)) { gammanames[i]<-paste("gamma",i,sep="") }
        colnames(Zw) <- gammanames
        mX <- cbind(Zw,Xw$mat)      # on fusionne ces 2 composantes explicatives
               
        # Ajustement du mod�le
        anaMrd <- glm(Y~mX,family=poisson)

    #V�rification du bon ajustement du mod�le loglin�aire
    if(!anaMrd$converged) stop("The algorithm did not converged")
    if(any(is.na(anaMrd$coef))) warning("Some loglinear parameter estimations cannot be evaluated")


#-----------------------------------------------------------------------------------#
# R�ajustement du mod�le en enlevant les gamma et les eta (si mod�le Chao) n�gatifs #
#-----------------------------------------------------------------------------------#

        param<-anaMrd$coef
        ppositions <- 0
        if (neg)
        {
            # Vecteur d'indicatrices pour les param�tres d'int�r�t n�gatifs
            indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0)))
            for (i in 1:I)
            {
                if ((identical(vm[i],"Mh")||identical(vm[i],"Mth"))&&identical(vh[[i]],"Chao")) {
                    indic <- as.vector(c(indic,rep(0,Xw$nbparam[i]-vt[i]+2),ifelse(param[(2*I+sum(Xw$nbparam[1:i])-vt[i]+2):(2*I+sum(Xw$nbparam[1:i])-1)]<0,1,0)))
                } else {
                    indic <- as.vector(c(indic,rep(0,Xw$nbparam[i])))
                }
            }
            while(isTRUE(sum(indic)>0)) # R�p�ter la boucle jusqu'� ce qu'aucun gamma appropri� ne soit n�gatif
            {
                # D�termination de la position du premier gamma appropri� n�gatif
                pos <- 1
                while(isTRUE(all.equal(indic[pos],0))) pos <- pos + 1
                ppositions <- c(ppositions,pos)
                # Retrait de la bonne colonne de mX et r�ajustement du mod�le
                mX <- mX[,-(pos-sum(ppositions<pos))]
                anaMrd <- glm(Y~mX,family=poisson)        
                # Ajout de z�ros dans le vecteur des param�tres loglin�aires
                positions <- sort(ppositions[-1])                
                param <- c(anaMrd$coef[1:(positions[1]-1)],0)
                if(isTRUE(length(positions)>1))
                {
                    for ( i in 2:length(positions))
                    {
                        if(isTRUE(all.equal(positions[i],positions[i-1]+1))) {
                            param <- c(param,0)
                        } else {
                            param <- c(param,anaMrd$coef[(positions[i-1]-i+2):(positions[i]-i)],0)
                        }
                    }
                }
                param <- c(param,anaMrd$coef[(positions[length(positions)]-length(positions)+1):length(anaMrd$coef)])
                # Vecteur d'indicatrices pour les param�tres d'int�r�t n�gatifs
                indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0)))
                for (i in 1:I)
                {
                    if ((identical(vm[[i]],"Mh")||identical(vm[[i]],"Mth"))&&identical(vh[[i]],"Chao")) {
                        indic <- as.vector(c(indic,rep(0,Xw$nbparam[i]-vt[[i]]+2),ifelse(param[(2*I+sum(Xw$nbparam[1:i])-vt[[i]]+2):(2*I+sum(Xw$nbparam[1:i])-1)]<0,1,0)))
                    } else {
                        indic <- as.vector(c(indic,rep(0,Xw$nbparam[i])))
                    }
                }
            }
        }
        positions <- sort(ppositions[-1]) 


#------------------------------------------------------------------------#
# Ajustement d'un mod�le pour tester la pr�sence d'�migration temporaire #
#------------------------------------------------------------------------#

        Inono <- length(vm[vm!="none"])
        idpemig <- 0
        for (i in 2:(I-1)) if(!identical(vm[i],"none")) idpemig <- c(idpemig,i)
        idpemig <- idpemig[-1]
        
        if(isTRUE(all.equal(Inono,3)))
        {
            anaMrd2 <- NULL
            parap2 <- NULL

            mX3<-cbind(mX,Xdelta[,idpemig])
            anaMrd3 <- glm(Y~mX3,family=poisson)
            parap3<-summary(anaMrd3)$coef[length(anaMrd3$coef),1:2]
        } else if(isTRUE(Inono>3))
            {
                mX2<-cbind(mX,apply(Xdelta[,idpemig],1,sum))
                anaMrd2 <- glm(Y~mX2,family=poisson)
                parap2<-summary(anaMrd2)$coef[length(anaMrd2$coef),1:2]

                mX3<-cbind(mX,Xdelta[,idpemig])
                anaMrd3 <- glm(Y~mX3,family=poisson)
                nn<-dim(summary(anaMrd3)$coef)[1]
                parap3<-summary(anaMrd3)$coef[(nn-Inono+3):nn,1:2]
            } else {
                anaMrd2 <- NULL
                parap2 <- NULL
    
                anaMrd3 <- NULL
                parap3 <- NULL
            }
     

#---------------------------------------#
# Formation des vecteurs des param�tres #
#---------------------------------------#

        # valeur de l intercept
        interc <- param[1]


        # creation des vecteurs de parametres alpha et beta
        Alpha <-rep(0,2*I-2)
        for (i in (1:(2*I-2)))
        {
                Alpha[i] <- param[i+1]
        }

        Beta <- rep(0,length(Xw$mat[1,]))
        for (i in (1:length(Xw$mat[1,])))
        {
                Beta[i] <- param[2*I-1+i]
        }


        # V�rification de la pr�sence de param�tres gamma n�gatifs si l'option "neg"=FALSE
        if(!neg)
        {
            if (any(Alpha<0)) warning(paste("One or more gamma parameters are negative.","\n",
                        "You can set them to zero with the 'neg' option.","\n",sep=""))
        }


#--------------------------------------------------------------#
# Matrice de variances-covariances des param�tres loglin�aires #
#--------------------------------------------------------------#

        if(isTRUE(length(positions)>0))
        {
            # Insertion de colonnes de z�ros
            varcovc <- cbind(summary(anaMrd)$cov.unscaled[,1:(positions[1]-1)],rep(0,dim(summary(anaMrd)$cov.unscaled)[1]))
            if(isTRUE(length(positions)>1))
            {
                for ( i in 2:length(positions))
                {
                    if(isTRUE(all.equal(positions[i],positions[i-1]+1))) {
                        varcovc <- cbind(varcovc,rep(0,dim(summary(anaMrd)$cov.unscaled)[1]))
                    } else {
                        varcovc <- cbind(varcovc,summary(anaMrd)$cov.unscaled[,(positions[i-1]-i+2):(positions[i]-i)],rep(0,dim(summary(anaMrd)$cov.unscaled)[1]))
                    }
                }
            }
            varcovc <- cbind(varcovc,summary(anaMrd)$cov.unscaled[,(positions[length(positions)]-length(positions)+1):dim(summary(anaMrd)$cov.unscaled)[2]])
            # Insertion de lignes de z�ros
            varcov <- rbind(varcovc[1:(positions[1]-1),],rep(0,dim(varcovc)[2]))
            if(isTRUE(length(positions)>1))
            {
                for ( i in 2:length(positions))
                {
                    if(isTRUE(all.equal(positions[i],positions[i-1]+1))) {
                        varcov <- rbind(varcov,rep(0,dim(varcovc)[2]))
                    } else {
                        varcov <- rbind(varcov,varcovc[(positions[i-1]-i+2):(positions[i]-i),],rep(0,dim(varcovc)[2]))
                    }
                }
            }
            varcov <- rbind(varcov,varcovc[(positions[length(positions)]-length(positions)+1):dim(varcovc)[1],]) 
        } else { varcov <- summary(anaMrd)$cov.unscaled }



############################################################################################################################################
# Estimation des param�tres d�mographiques
############################################################################################################################################

#--------------------------------------------#
# calcul des probabilites de capture (pstar) #
#--------------------------------------------#

        pstar <- rep(0,I)
        dpstar<-matrix(rep(0,I*length(param)),ncol=I)
        beta <- Beta
        for (i in (1:I))
        {
                if (identical(vm[i],"none")) # cas du no model
                {
                        pstar[i]<- exp(beta[1]+log(2^vt[i]-1))/(1+exp(beta[1]+log(2^vt[i]-1)))
                        dpstar[(2*I+ifelse(i>1,sum(Xw$nbparam[1:(i-1)]),0)):(2*I-1+sum(Xw$nbparam[1:i])),i] <- exp(beta[1]+log(2^vt[i]-1))/(1+exp(beta[1]+log(2^vt[i]-1)))^2
                        beta <- beta[-1]
                } else if (identical(vm[i],"Mt")) # cas du Modele Mt
                {
                        pstar[i] <- 1-prod((1+exp(beta[c(1:vt[i])]))^-1)
                        dpstar[(2*I+ifelse(i>1,sum(Xw$nbparam[1:(i-1)]),0)):(2*I-1+sum(Xw$nbparam[1:i])),i] <- prod((1+exp(beta[c(1:vt[i])]))^-1)*exp(beta[c(1:vt[i])])/(1+exp(beta[c(1:vt[i])]))
                        beta <- beta[-c(1:vt[i])]
                } else # Tous les autres mod�les
                {
                        Xpf <- Xclosedp(vt[i],vm[i],vh[[i]],va[i])
                        pstar[i] <- sum(exp(Xpf$mat%*%beta[c(1:Xpf$nbparam)]))/(1+ sum(exp(Xpf$mat%*%beta[c(1:Xpf$nbparam)])))
                        dpstar[(2*I+ifelse(i>1,sum(Xw$nbparam[1:(i-1)]),0)):(2*I-1+sum(Xw$nbparam[1:i])),i] <- t(Xpf$mat)%*%exp(Xpf$mat%*%beta[c(1:Xpf$nbparam)])/(1+sum(exp(Xpf$mat%*%beta[c(1:Xpf$nbparam)])))^2
                        beta <- beta[-c(1:Xpf$nbparam)]
                }
        }
        varcovpstar <- t(dpstar)%*%varcov%*%dpstar
        pstarStderr <- sqrt(diag(varcovpstar))    


#---------------#
# calcul des Ui #
#---------------#

        uv <- rep(1,(I-1))
        duv<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        for (i in (1:(I-1)))
        {
                eAlpha <- exp(Alpha[I:(I+i-1)])
                unmpstar <- (1-pstar[I:(I-i+1)])
                uv[i] <- prod(eAlpha*unmpstar)*(1-exp(-Alpha[I+i-1]))
                duv[,i] <- prod(eAlpha*unmpstar)*c(rep(0,I+i-1),exp(-Alpha[I+i-1]),rep(0,(length(param)-I-i)))
                for (j in 1:i)
                {
                    duv[,i] <- duv[,i] + (1-exp(-Alpha[I+i-1]))*prod(eAlpha[-j]*unmpstar[-j])*((1-pstar[I-j+1])*c(rep(0,I+j-1),exp(Alpha[I+j-1]),rep(0,(length(param)-I-j)))-exp(Alpha[I+j-1])*dpstar[,I-j+1])
                }
        }
        uv <- c(1,uv)
        duv <- cbind(rep(0,length(param)),duv)
        varcovuv <- t(duv)%*%varcov%*%duv
        uvStderr <- sqrt(diag(varcovuv))    


#---------------#
# calcul des Vi #
#---------------#

        vv <- rep(1,(I-1))
        dvv<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        for (i in (1:(I-1)))
        {
                eAlpha <- exp(Alpha[1:i])
                unmpstar <- (1-pstar[1:i])
                vv[i] <- prod(eAlpha*unmpstar)*(1-exp(-Alpha[i]))
                dvv[,i] <- prod(eAlpha*unmpstar)*c(rep(0,i),exp(-Alpha[i]),rep(0,(length(param)-i-1)))
                for (j in 1:i)
                {
                    dvv[,i] <- dvv[,i] + (1-exp(-Alpha[i]))*prod(eAlpha[-j]*unmpstar[-j])*((1-pstar[j])*c(rep(0,j),exp(Alpha[j]),rep(0,(length(param)-j-1)))-exp(Alpha[j])*dpstar[,j])
                }
        }
        vv <- c(1,vv)
        dvv <- cbind(rep(0,length(param)),dvv)
        varcovvv <- t(dvv)%*%varcov%*%dvv
        vvStderr <- sqrt(diag(varcovvv))    


#--------------------------------------------------------------#
# calcul des probabilites de survie entre chaque periode (phi) #
#--------------------------------------------------------------#

        phi <- rep(0,(I-1))
        dphi<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        phi[(I-1):1] <-1/(1+uv[2:I]/cumsum(uv[1:(I-1)]))   
        for ( i in 1:(I-1))
        {
            if(isTRUE(all.equal(i,I-1)))
            {
            dphi[,i] <- -(phi[i]^2)*duv[,I-i+1]
            } else {
            dphi[,i] <- -(phi[i]^2)*(duv[,I-i+1]*sum(uv[1:(I-i)])-uv[I-i+1]*rowSums(duv[,1:(I-i)]))/sum(uv[1:(I-i)])^2
            }
        }
        varcovphi <- t(dphi)%*%varcov%*%dphi
        phistderr <- sqrt(diag(varcovphi))    
  

#---------------------------------------#
# calcul des taille de populations (Ni) #
#---------------------------------------#

        Npop <- rep(0,I)
        dNpop<-matrix(rep(0,I*length(param)),ncol=I)
        
        Npop[1] <- exp(interc)/(prod(1-pstar)*prod(phi))   
        dprodpstar <-rep(0,length(param))
        for (i in 1:I)
        {
            dprodpstar<-dprodpstar-prod(1-pstar[-i])*dpstar[,i]
        }
        dprodphi <-rep(0,length(param))
        for (i in 1:(I-1))
        {
            dprodphi<-dprodphi+prod(phi[-i])*dphi[,i]
        }
        dNpop[,1] <- (prod(1-pstar)*prod(phi)*c(exp(interc),rep(0,(length(param)-1))) - exp(interc)*(prod(phi)*dprodpstar+prod(1-pstar)*dprodphi))/ (prod(1-pstar)*prod(phi))^2
        
        Npop[2:I]<-Npop[1]*cumprod((vv[2:I]/cumsum(vv[1:(I-1)])+1)*phi)
        for (i in 2:I)
        {
            if(isTRUE(all.equal(i,2))) {
                dNpop[,i] <- phi[i-1]*(vv[i]+1)*dNpop[,i-1] + Npop[i-1]*(vv[i]+1)*dphi[,i-1] + Npop[i-1]*phi[i-1]*dvv[,i]          
            } else {
                dNpop[,i] <- phi[i-1]*(vv[i]/sum(vv[1:(i-1)])+1)*dNpop[,i-1] + Npop[i-1]*(vv[i]/sum(vv[1:(i-1)])+1)*dphi[,i-1] + Npop[i-1]*phi[i-1]*(sum(vv[1:(i-1)])*dvv[,i]-vv[i]*rowSums(dvv[,1:(i-1)]))/sum(vv[1:(i-1)])^2
            }
        }
        varcovtpop <- t(dNpop)%*%varcov%*%dNpop
        NpopStderr <- sqrt(diag(varcovtpop))    
        NpopStderr <- sqrt(pmax(NpopStderr^2-Npop,0))          


#----------------------------#
# calcul des naissances (Bi) #
#----------------------------#

        B<-Npop[2:I]-Npop[1:(I-1)]*phi
        dB <- dNpop[,2:I] - t(phi*t(dNpop[,1:(I-1)])) - t(Npop[1:(I-1)]*t(dphi))
        varcovB <- t(dB)%*%varcov%*%dB
        BStderr <- sqrt(diag(varcovB))
        
        
#--------------------------------------------------------------------#
# Calcul du nombre total d'individus qui ont pass� sur le territoire #
#--------------------------------------------------------------------#

        Ntot <- Npop[1]+sum(B)
        dNtot <- dNpop[,1] + rowSums(dB)    
        NtotStderr <- sqrt(max(t(dNtot)%*%varcov%*%dNtot-Ntot,0))


############################################################################################################################################
# Pr�sentation des r�sultats
############################################################################################################################################

        modelfit <- matrix(c(anaMrd$deviance,anaMrd$df.residual,anaMrd$aic),nrow=1)
        dimnames(modelfit) <- list("fitted model",c("deviance","    df","      AIC"))
        
        emigfit <- cbind(c(anaMrd2$deviance,anaMrd3$deviance),c(anaMrd2$df.residual,anaMrd3$df.residual),c(anaMrd2$aic,anaMrd3$aic))
        if (isTRUE(all.equal(Inono,3))) { dimnames(emigfit) <- list(c("model with temporary emigration"),c("deviance","    df","      AIC"))
        } else if (isTRUE(Inono>3)) { dimnames(emigfit) <- list(c("model with homogeneous temporary emigration","model with temporary emigration"),c("deviance","    df","      AIC"))}    
        
        titre.periode.emig<-rep(0,length(idpemig))
        for (i in 1:length(idpemig)){ titre.periode.emig[i]<-paste("period",idpemig[i])}
        titre.periode<-rep(0,I)
        for (i in 1:I){titre.periode[i]<-paste("period",i)}
        titre.inter.periode<-rep(0,I-1)
        for (i in 1:(I-1)){titre.inter.periode[i]<-paste("period",i,"->",i+1)}
        titre.i<-rep(0,I)
        for (i in 1:I){titre.i[i]<-paste("i =",i-1)}
        
        parap <- rbind(parap3,parap2)     
        pstar <- cbind(pstar,pstarStderr)
        phi <- cbind(phi,phistderr)
        Npop <- cbind(Npop,NpopStderr)
        B <- cbind(round(B,digit=6),round(BStderr,digit=6))
        Ntot <- cbind(Ntot,NtotStderr)
        loglinearpara <- cbind(param,sqrt(diag(varcov)))
        uv <- cbind(uv,uvStderr)
        vv <- cbind(vv,vvStderr)

        if (isTRUE(all.equal(Inono,3))) { dimnames(parap) <- list(titre.periode.emig,c("estimate","stderr")) 
        } else if (isTRUE(Inono>3)) { dimnames(parap) <- list(c(titre.periode.emig,"homogenous"),c("estimate","stderr")) }
        dimnames(pstar) <- list(titre.periode,c("estimate","stderr"))
        dimnames(phi) <- list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Npop) <- list(titre.periode,c("estimate","stderr"))
        dimnames(B) <- list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Ntot) <- list("all periods",c("estimate","stderr"))
        dimnames(loglinearpara) <- list(c("intercept",gammanames,Xw$paramnames),c("estimate","stderr"))
        dimnames(uv) <- list(titre.i,c("estimate","stderr"))
        dimnames(vv) <- list(titre.i,c("estimate","stderr"))

        models<-matrix(Xw$models,nrow=1)
        dimnames(models) <- list("model",titre.periode)

        # Matrice de variances-covariances des param�ters pstar, phi, Npop, B et Ntot
        dP <- cbind(dpstar,dphi,dNpop,dB,dNtot)
        covP <- t(dP)%*%varcov%*%dP
        titre.P <- c(rep(0,4*I-2),"Ntot")
        for (i in 1:I){
            titre.P[i]<-paste("p*",i)
            titre.P[2*I-1+i] <- paste("Npop",i)
        }
        for (i in 1:(I-1)){
            titre.P[I+i]<-paste("phi",i)
            titre.P[3*I-1+i] <- paste("B",i)
        }
        dimnames(covP)<-list(titre.P,titre.P)
      
      
        ans<-list(n=sum(Y),models=models,model.fit=modelfit,emig.fit=emigfit,emig.param=parap,
                  capture.prob=pstar,survivals=phi,N=Npop,birth=B,Ntot=Ntot,
                  loglin.param=loglinearpara,u.vector=uv,v.vector=vv,cov=covP, neg=positions)
        class(ans) <- "robustd"
        ans                
}


print.robustd <- function(x, ...){
        test <- rep(x$models[1],length(x$models))
        if(isTRUE(all.equal(test,as.character(x$models))))
        {
            cat("\nClosed population model for every period:",x$models[1],"\n")
        } else {
            cat("\nClosed population models for each period:\n")
            print.default(x$models, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")        
        }
        cat("\nModel fit:\n")
        x$model.fit[c(1,3)] <- round(x$model.fit[c(1,3)],3)
        x$model.fit[2] <- round(x$model.fit[2],0)
        print.default(x$model.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nTest for temporary emigration:\n")
        x$emig.fit[,c(1,3)] <- round(x$emig.fit[,c(1,3)],3)
        x$emig.fit[,2] <- round(x$emig.fit[,2],0)
        print.default(x$emig.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nCapture probabilities:\n")
        print.default(round(x$capture.prob,4), print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nSurvival probabilities:\n")
        print.default(round(x$survivals,4), print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nAbundances:\n")
        x$N <- round(x$N,1)
        print.default(x$N, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nNumber of new arrivals:\n")
        x$birth <- round(x$birth,1)
        print.default(x$birth, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nTotal number of units who ever inhabited the survey area:\n")
        x$Ntot <- round(x$Ntot,1)
        print.default(x$Ntot, print.gap = 2, quote = FALSE, right=TRUE, na.print="--")
        cat("\nTotal number of captured units:",x$n,"\n")
        if (isTRUE(all.equal(length(x$neg[x$neg<2*dim(x$N)[1]]),1))) cat("\nNote:",length(x$neg[x$neg<2*dim(x$N)[1]]),"gamma parameter has been set to zero")
        if (isTRUE(length(x$neg[x$neg<2*dim(x$N)[1]])>1)) cat("\nNote:",length(x$neg[x$neg<2*dim(x$N)[1]]),"gamma parameters has been set to zero")
        if (isTRUE(all.equal(length(x$neg[x$neg>=2*dim(x$N)[1]]),1))) cat("\nNote:",length(x$neg[x$neg>=2*dim(x$N)[1]]),"eta parameter has been set to zero")
        if (isTRUE(length(x$neg[x$neg>=2*dim(x$N)[1]])>1)) cat("\nNote:",length(x$neg[x$neg>=2*dim(x$N)[1]]),"eta parameters has been set to zero")
        cat("\n\n")
        invisible(x)
}