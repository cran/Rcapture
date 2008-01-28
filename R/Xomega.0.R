"Xomega.0" <-
function(vt,vm,vh,va,rd.call)
{

        histpos <- histpos.0(vt)  # matrice des historiques de captures possibles pour le nombre d'occasions de capture total
        I <- length(vt)   # nombre de periodes primaires
        nrows <- dim(histpos)[1]
        M <- rep(0,nrows)        # aide a la fabrication de la matrice
        nbparam <- rep(0,I)
        models <- rep(0,I)


        # on cree la matrice periode par periode en respectant les models demandes en entree
        for (i in (1:I))
        {
            # Création de vecteurs de noms qui permettent de rendre plus clair les objects glm fournis en sortie
            beta <- paste("beta",i,".",1,sep="")
            etanames <- vector("character",vt[i]-2)
            for (j in 3:vt[i]) { etanames[j-2]<-paste("eta",i,".",j,sep="") }

            # selection des colonnes correspondantes a la periode etudiee dans cette boucle
            histposp<-histpos[,i]

            # creation de matrice selon le modele choisi
            if (vm[i]=="none") # no modele
            {
                    mXp <- as.matrix(ifelse(histposp>0,1,0))
                    colnames(mXp) <- beta
                    models[i] <- "none"
            } else if (vm[i]=="M0") # modele M0
                    {
                        mXp <- as.matrix(histposp)     # la matrice reste la meme
                        colnames(mXp) <- beta
                        models[i] <- "M0"
                    } else  if (vm[i]=="Mh") # modele Mh
                            {
                                if (vh[[i]]=="Chao")
                                {
                                    mXp2 <- matrix(0,nrows,vt[i]-2)
                                    for (j in (3:vt[i])) { mXp2[,j-2]<-pmax(histposp-j+1,0) }
                                } else
                                if (vh[[i]]=="Poisson") mXp2 <- va[i]^histposp - 1 else
                                if (vh[[i]]=="Darroch") mXp2 <- (histposp^2)/2 else
                                if(is.function(vh[[i]])) mXp2 <- vh[[i]](histposp)
                                mXp <- cbind(histposp,mXp2)
                                colnames(mXp) <- if (vh[[i]]=="Chao") c(beta,etanames) else c(beta,paste("tau",i,".",1,sep=""))
                                models[i] <- if (is.function(vh[[i]])) paste("Mh",rd.call$vh[i+1]) else if(vh[[i]]=="Poisson") paste("Mh",paste(vh[[i]],va[i],sep="")) else paste("Mh",vh[[i]])
                            }
            M <- cbind(M,mXp)
            nbparam[i] <- dim(mXp)[2]
        }

        M <- M[,-1]     # on supprime la premiere de colonne de zero inutile du debut

        list(mat=M,nbparam=nbparam,models=models,paramnames=colnames(M))
}
