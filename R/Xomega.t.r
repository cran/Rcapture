"Xomega.t" <-
function(vt,vm,vh,va,rd.call)
{

        histpos <- histpos.t(sum(vt))  # matrice des historiques de captures possibles pour le nombre d'occasions de capture total
        I <- length(vt) # nombre de periodes primaires
        nrows <- dim(histpos)[1]
        M <- rep(0,nrows)        # aide a la fabrication de la matrice
        nbparam <- rep(0,I)
        models <- rep(0,I)


        # on cree la matrice periode par periode en respectant les models demandes en entree
        for (i in (1:I))
        {
            # Création de vecteurs de noms qui permettent de rendre plus clair les objects glm fournis en sortie
            betanames <- vector("character",vt[i])
            for (j in 1:vt[i]) { betanames[j]<-paste("beta",i,".",j,sep="") }
            etanames <- vector("character",vt[i]-2)
            for (j in 3:vt[i]) { etanames[j-2]<-paste("eta",i,".",j,sep="") }

            # selection des colonnes correspondantes a la periode etudiee dans cette boucle
            if (isTRUE(all.equal(i,1))) { histposp <- histpos[,c(1:vt[i])]   
            } else { histposp <- histpos[,c((sum(vt[1:(i-1)])+1):sum(vt[1:i]))] }

            # creation de la matrice Xomega selon le modele choisi
            if (identical(vm[i],"none")) # no model
            {
                mXp <- as.matrix(apply(histposp,1,max))
                colnames(mXp) <- betanames[1]
                models[i] <- "none"
            } else if (identical(vm[i],"M0")) # modele M0
                    {
                        mXp <- as.matrix(apply(histposp, 1, sum))
                        colnames(mXp) <- betanames[1]
                        models[i] <- "M0"
                    } else if (identical(vm[i],"Mt")) # modele Mt
                            {
                                mXp <- histposp     # la matrice reste la meme
                                colnames(mXp) <- betanames
                                models[i] <- "Mt"
                            } else if (identical(vm[i],"Mh")) # modele Mh
                                    {
                                        mXp1 <- apply(histposp,1,sum) # nombre de capture par historique possible
                                        if (identical(vh[[i]],"Chao"))
                                        {
                                            mXp2 <- matrix(0,nrows,vt[i]-2)
                                            for (j in (3:vt[i])) { mXp2[,j-2]<-pmax(mXp1-j+1,0) }
                                        } else
                                        if (identical(vh[[i]],"Poisson")) mXp2 <- va[i]^mXp1 - 1 else
                                        if (identical(vh[[i]],"Darroch")) mXp2 <- (mXp1^2)/2 else
                                        if(is.function(vh[[i]])) mXp2 <- vh[[i]](mXp1)
                                        mXp <- cbind(mXp1,mXp2)
                                        colnames(mXp) <- if (identical(vh[[i]],"Chao")) c(betanames[1],etanames) else c(betanames[1],paste("tau",i,".",1,sep=""))
                                        models[i] <- if (is.function(vh[[i]])) paste("Mh",rd.call$vh[i+1]) else if(identical(vh[[i]],"Poisson")) paste("Mh",paste(vh[[i]],va[i],sep="")) else paste("Mh",vh[[i]])
                                    } else if (identical(vm[i],"Mth")) # modele Mth
                                            {
                                                mXp1 <- apply(histposp,1,sum)
                                                if (identical(vh[[i]],"Chao"))
                                                {
                                                    mXp2 <- matrix(0,nrows,vt[i]-2)
                                                    for (j in (3:vt[i])) { mXp2[,j-2]<-pmax(mXp1-j+1,0) }
                                                } else
                                                if (identical(vh[[i]],"Poisson")) mXp2 <- va[i]^mXp1 - 1 else
                                                if (identical(vh[[i]],"Darroch")) mXp2 <- (mXp1^2)/2 else
                                                if(is.function(vh[[i]])) mXp2 <- vh[[i]](mXp1)
                                                mXp <- cbind(histposp,mXp2)
                                                colnames(mXp) <- if (identical(vh[[i]],"Chao")) c(betanames,etanames) else c(betanames,paste("tau",i,".",1,sep=""))
                                                models[i] <- if (is.function(vh[[i]])) paste("Mth",rd.call$vh[i+1]) else if(identical(vh[[i]],"Poisson")) paste("Mth",paste(vh[[i]],va[i],sep="")) else  paste("Mth",vh[[i]])
                                            }
            M <- cbind(M,mXp)
            nbparam[i] <- dim(mXp)[2]
        }

        M <- M[,-1]     # on supprime la premiere colonne de zero inutile du debut

        list(mat=M,nbparam=nbparam,models=models,paramnames=colnames(M))
}
