"descriptive" <- function(X,dfreq=FALSE,dtype=c("hist","nbcap"),t)
{
    ############################################
    # Validation des arguments fournis en entrée
    valid.one(dfreq,"logical")
    dtype <- dtype[1]
    if(!dtype%in%c("hist","nbcap")) stop("'dtype' must be given the value 'hist' or 'nbcap'")
    if (missing(t)) t <- NULL else valid.one(t,"numeric")
    Xvalid<-valid.X(X,dfreq,dtype,t,warn=FALSE)
        X <- Xvalid$X
        t <- Xvalid$t    
    ############################################

# Statistiques descriptives de base générées automatiquement    
        
    Y <- histfreq.0(X=X,dfreq=dfreq,dtype=dtype,vt=t) # Nombre d'individus capturés i fois
    nbrecapt <- rev(Y)
    if (dtype=="hist") {        
          premcapt <- getui(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés pour la première fois à l'occasion i        
          derncapt <- getvi(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés pour la dernière fois à l'occasion i
          captoccas <- getni(X=X,dfreq=dfreq,t=t) # Nombre d'individus capturés à l'occasion i
    }       
        
    titre.i<-paste("i =",1:t)
    if (dtype=="hist") {        
        tableau<-cbind(nbrecapt,premcapt,derncapt,captoccas)
        dimnames(tableau) <- list(titre.i,c("fi","ui","vi","ni"))
    } else {
        tableau<-matrix(nbrecapt,ncol=1)
        dimnames(tableau) <- list(titre.i,"fi")
    }       
        
    nbre <- sum(na.rm=TRUE,tableau[,1]) # le nombre d unites etudiees dans cette matrice


# Matrice de recapture générée mais pas dans le print

    if (dtype=="hist") {    
        recap<-matrix(rep(NA,t*t),ncol=t)
        for(i in 1:(t-1)){
                Xsous<-matrix(X[X[,i]==1,],ncol=dim(X)[2])
                for(j in (i+1):t){
                        recap[i,j-1] <- if(dfreq) sum(na.rm=TRUE,Xsous[Xsous[,j]==1,t+1]) else sum(na.rm=TRUE,Xsous[,j])
                        Xsous<-matrix(Xsous[Xsous[,j]!=1,],ncol=dim(X)[2])
                }
                recap[i,t] <- if(dfreq) sum(na.rm=TRUE,Xsous[,t+1]) else dim(Xsous)[1]
        }
        recap[t,t]<-captoccas[t]

        table.recap<-cbind(captoccas,recap)
        dimnames(table.recap) <- list(titre.i,c("ni",paste("c",2:t,sep=""),"not recapt"))
    }

# Préparation des sorties

    ans <- if (dtype=="hist") list(n=nbre,base.freq=tableau,m.array=table.recap) else list(n=nbre,base.freq=tableau)
    class(ans) <- "descriptive"
    ans
}


print.descriptive <- function(x, ...)
{
        cat("\nNumber of captured units:",x$n,"\n","\nFrequency statistics:\n")
        print.default(format(x$base.freq), print.gap = 2, quote = FALSE, ...)
           cat("fi: number of units captured i times\n")
           if (dim(x$base.freq)[2]==4) {
               cat("ui: number of units captured for the first time on occasion i\n")
               cat("vi: number of units captured for the last time on occasion i\n")
               cat("ni: number of units captured on occasion i\n\n")
           }    
        invisible(x)
}


plot.descriptive <- function(x,main="Exploratory Heterogeneity Graph", ...){
        graph1<-cbind(1:dim(x$base.freq)[1],x$base.freq[,1],log(x$base.freq[,1]/choose(dim(x$base.freq)[1],1:dim(x$base.freq)[1])))
        graph1<-graph1[graph1[,2]!=0,,drop=FALSE]
        if (dim(x$base.freq)[2]==4) {
             graph2<-cbind(1:dim(x$base.freq)[1],x$base.freq[,2],log(x$base.freq[,2]))
             graph2<-graph2[graph2[,2]!=0,,drop=FALSE]
        }
        ngraph <- if (dim(x$base.freq)[2]==1||dim(graph1)[1]==1) 1 else 2 
        mar <- if (ngraph==2) c(3, 5.5, 5.5, 2) else c(5, 5.5, 5.5, 2)
        op<-par(mfrow=c(ngraph,1),mar=mar)
        on.exit(par(op))
        if (dim(graph1)[1]!=1||dim(x$base.freq)[2]==1) { # Si trop de fi nuls, on ne les illustre pas
             plot(graph1[,1],graph1[,3],type="b",ann=0,...)
             mtext("fi: number of units captured i times",side=3,line=0.5,adj=0,font=2)
             mtext(expression("log"*bgroup("(",frac("fi",bgroup("(", atop(n, i), ")")),")")),side=2,line=2.5,las=1)
             mtext("i: number of captures",side=1,line=2.5)
        }
        if (ngraph==2) mtext(main,side=3,line=2.7,cex=1.8) # S'il y a deux graphiqes, le texte doit être sjouté après le premier graphique
        if (dim(x$base.freq)[2]==4) {
             if (ngraph==2) par(mar=c(5, 5.5, 3.5, 2))
             plot(graph2[,1],graph2[,3],type="b",ann=0,...)
             mtext("ui: number of units captured for the first time on occasion i",side=3,line=0.5,adj=0,font=2)
             mtext("log(ui)",side=2,line=2.5,las=1)
             mtext("i: capture occasion identification number",side=1,line=2.5)
        }
        if (ngraph==1) mtext(main,side=3,line=2.7,cex=1.8)
}

#        op<-par(mfrow=c(ngraph,1),mar=c(3ou5, 5.5, 5.5, 2))
#        on.exit(par(op))
#        if (dim(graph1)[1]==1) { # Si trop de fréquences nulles
#             plot(graph2[,1],graph2[,3],type="b",ann=0,...)
#             mtext("ui: number of units captured for the first time on occasion i",side=3,line=0.5,adj=0,font=2)
#             mtext("log(ui)",side=2,line=2.5,las=1)
#             mtext("i: capture occasion identification number",side=1,line=2.5)
#        } else {          
#             op<-par(mfrow=c(2,1),mar=c(3, 5.5, 5.5, 2))
#             plot(graph1[,1],graph1[,3],type="b",ann=0,...)
#             mtext("fi: number of units captured i times",side=3,line=0.5,adj=0,font=2)
#             mtext(expression("log"*bgroup("(",frac("fi",bgroup("(", atop(n, i), ")")),")")),side=2,line=2.5,las=1)
#             mtext("i: number of captures",side=1,line=2.5)
#             par(mar=c(5, 5.5, 3.5, 2))
#             plot(graph2[,1],graph2[,3],type="b",ann=0,...)
#             mtext("ui: number of units captured for the first time on occasion i",side=3,line=0.5,adj=0,font=2)
#             mtext("log(ui)",side=2,line=2.5,las=1)
#             mtext("i: capture occasion identification number",side=1,line=2.5)
#        }
#        mtext(main,side=3,line=2.7,cex=1.8)
#


##################################################################################################
## Sous-fonctions pour le calcul des certaines stat descriptives

getfi <- function(X,dfreq,t)
{ 
        # Nombre d'individus capturés i fois
        nbrecapt <- rep(0,t) # On veut avoir les fréquences pour tous les nbcap, même les fréquences nulles
        for (i in (1: dim(X[,1:t])[1])) {
                v <- sum(na.rm=TRUE,X[i,1:t])
                nbrecapt[v] <- nbrecapt[v] + if(dfreq) X[i,t+1] else 1
        }
        return(nbrecapt)
}

getui <- function(X,dfreq,t)
{ 
        # Nombre d'individus capturés pour la première fois à l'occasion i
        premcapt <- rep(0,t)
        for(i in (1:dim(X[,1:t])[1])) {
                k <- 0
                j <- 1
                while(k==0&&j<=t) {
                        if (X[i,j]==1) {
                                k<-1
                                premcapt[j]<- premcapt[j] + if(dfreq) X[i,t+1] else 1
                        } else {
                                j<-j+1
                        }
                }
        }
        return(premcapt)
}

getvi <- function(X,dfreq,t)
{ 
        # Nombre d'individus capturés pour la dernière fois à l'occasion i
        derncapt <- rep(0,t)
        for(i in (1:dim(X[,1:t])[1])) {
                k <- 0
                j <- t
                while(k==0&&j>=1) {
                        if (X[i,j]==1) {
                                k<-1
                                derncapt[j]<- derncapt[j] + if(dfreq) X[i,t+1] else 1
                        } else {
                                j<-j-1
                        }
                }
        }
        return(derncapt)
}

getni <- function(X,dfreq,t)
{ 
        # Nombre d'individus capturés à l'occasion i
        if (dfreq) {
            captoccas <- rep(0,t)
            for (i in 1:t) { captoccas[i] <- sum(na.rm=TRUE,X[X[,i]==1,t+1]) }
        } else captoccas <- colSums(X)
        return(captoccas)
}
