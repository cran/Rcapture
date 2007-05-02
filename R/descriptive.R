"descriptive" <- function(X,dfreq=FALSE)
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


# Statistiques descriptives de base générées automatiquement

        # Nombre d'individus capturés i fois
        nbrecapt <- rep(0,t)
        for (i in (1: dim(X[,1:t])[1]))
        {
                v <- sum(X[i,1:t])
                nbrecapt[v] <- nbrecapt[v] + ifelse(dfreq,X[i,t+1],1)
        }

        # Nombre d'individus capturés pour la première fois à l'occasion i
        premcapt <- rep(0,t)
        for(i in (1:dim(X[,1:t])[1]))
        {
                k <- 0
                j <- 1
                while(isTRUE(all.equal(k,0))&&isTRUE(j<=t))
                {
                        if (isTRUE(all.equal(X[i,j],1)))
                        {
                                k<-1
                                premcapt[j]<- premcapt[j]+ifelse(dfreq,X[i,t+1],1)
                        } else
                        {
                                j<-j+1
                        }
                }
        }


        # Nombre d'individus capturés pour la dernière fois à l'occasion i
        derncapt <- rep(0,t)
        for(i in (1:dim(X[,1:t])[1]))
        {
                k <- 0
                j <- t
                while(isTRUE(all.equal(k,0))&&isTRUE(j>=1))
                {
                        if (isTRUE(all.equal(X[i,j],1)))
                        {
                                k<-1
                                derncapt[j]<- derncapt[j]+ifelse(dfreq,X[i,t+1],1)
                        } else
                        {
                                j<-j-1
                        }
                }
        }


        # Nombre d'individus capturés à l'occasion i
        if (dfreq)
        {
            captoccas <- rep(0,t)
            for (i in 1:t) { captoccas[i] <- sum(X[X[,i]==1,t+1]) }
        } else captoccas <- apply(X,2,sum)
        
        tableau<-cbind(nbrecapt,premcapt,derncapt,captoccas)
        titre.i<-rep(0,t)
        for (i in 1:t){titre.i[i]<-paste("i =",i)}
        dimnames(tableau) <- list(titre.i,c("fi","ui","vi","ni"))


        # le nombre d unites etudiees dans cette matrice
        nbre <- sum(tableau[,1])


# Matrice de recapture générée mais pas dans le print

        recap<-matrix(rep(NA,t*t),ncol=t)
        for(i in 1:(t-1))
        {
                Xsous<-matrix(X[X[,i]==1,],ncol=dim(X)[2])
                for(j in (i+1):t)
                {
                        recap[i,j-1]<-ifelse(dfreq,sum(Xsous[Xsous[,j]==1,t+1]),sum(Xsous[,j]))
                        Xsous<-matrix(Xsous[Xsous[,j]!=1,],ncol=dim(X)[2])
                }
                recap[i,t]<-ifelse(dfreq,sum(Xsous[,t+1]),dim(Xsous)[1])
        }
        recap[t,t]<-captoccas[t]

        table.recap<-cbind(captoccas,recap)
        c.i<-rep(0,t-1)
        for (i in 2:t){c.i[i-1]<-paste("c",i,sep="")}
        dimnames(table.recap) <- list(titre.i,c("ni",c.i,"not recapt"))


# Préparation des sorties

        ans<-list(n=nbre,base.freq=tableau,m.array=table.recap)
        class(ans) <- "descriptive"
        ans
}


print.descriptive <- function(x, ...)
{
        cat("\nNumber of captured units:",x$n,"\n","\nFrequency statistics:\n")
        print.default(format(x$base.freq), print.gap = 2, quote = FALSE)
                cat("fi: number of units captured i times\n")
                cat("ui: number of units captured for the first time on occasion i\n")
                cat("vi: number of units captured for the last time on occasion i\n")
                cat("ni: number of units captured on occasion i\n\n")
        invisible(x)
}


plot.descriptive <- function(x, ...){
        graph1<-cbind(1:dim(x$base.freq)[1],x$base.freq[,1],log(x$base.freq[,1]/choose(dim(x$base.freq)[1],1:dim(x$base.freq)[1])))
        graph1<-graph1[graph1[,2]!=0,]
        graph2<-cbind(1:dim(x$base.freq)[1],x$base.freq[,2],log(x$base.freq[,2]))
        graph2<-graph2[graph2[,2]!=0,]
        par(mfrow=c(2,1),mar=c(3, 5.5, 5.5, 2))
        plot(graph1[,1],graph1[,3],type="b",ann=0,xaxp=c(1,x$n,x$n-1))
        mtext("Exploratory Heterogeneity Graphs",side=3,line=2.7,cex=1.8)      
        mtext("fi: number of units captured i times",side=3,line=0.5,adj=0,font=2)
        mtext(expression("log"*bgroup("(",frac("fi",bgroup("(", atop(n, i), ")")),")")),side=2,line=2.5,las=1)
        mtext("i: number of captures",side=1,line=2.5)
        par(mar=c(5, 5.5, 3.5, 2))
        plot(graph2[,1],graph2[,3],type="b",ann=0,xaxp=c(1,x$n,x$n-1))
        mtext("ui: number of units captured for the first time on occasion i",side=3,line=0.5,adj=0,font=2)
        mtext("log(ui)",side=2,line=2.5,las=1)
        mtext("i: capture occasion identification number",side=1,line=2.5)
        par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
}
