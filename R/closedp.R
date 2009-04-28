"closedp" <- "closedp.t" <- function(X,dfreq=FALSE,neg=TRUE,trace=FALSE)
{    
     eval(valid.t)
     eval(valid.closedp)
     eval(model.t)
     # Identification des modèles à ajuster
     if (t==2) {
          lmn <- smn <- c("M0","Mt","Mb")
     } else {
          lmn <- c("M0","Mt","Mh Chao","Mh Poisson2","Mh Darroch","Mh Gamma3.5","Mth Chao","Mth Poisson2","Mth Darroch","Mth Gamma3.5","Mb","Mbh")
          smn <- c("M0","Mt","MhC","MhP","MhD","MhG","MthC","MthP","MthD","MthG","Mb","Mbh")
     }
     eval(closedp.internal)
     class(ans) <- c("closedp","closedp.t")
     return(ans)
}

"closedp.0" <- function(X,dfreq=FALSE,dtype=c("hist","nbcap"),t,t0=t,neg=TRUE,trace=FALSE)
{    
     eval(valid.0)
     eval(valid.closedp)
     eval(model.0)
     # Identification des modèles à ajuster
     if (t0==2) {
          lmn <- smn <- c("M0")
     } else {
          lmn <- c("M0","Mh Chao","Mh Poisson2","Mh Darroch","Mh Gamma3.5")
          smn <- c("M0","MhC","MhP","MhD","MhG")
     }
     eval(closedp.internal)
     ansf <- c(ans,list(t0=t0))
     class(ansf) <- "closedp"
     return(ansf)
}

# Bouts de code communs aux deux fonctions : expressions R prêtes à être évaluées

valid.closedp <- quote({ # Validation des arguments des fonctions closedp
     valid.one(dfreq,"logical")
     Xvalid<-valid.X(X,dfreq,dtype,t,warn=typet)
         X <- Xvalid$X
         t <- Xvalid$t
     if (t<2) stop("the number of capture occasions 't' must be at least 2")    
     if (!typet) {    
         valid.one(t0,"numeric")
         if ((t0%%1)!=0||t0>t||t0<2) 
               stop("'t0' must be an integer between 't', the number of capture occasion, and 2 inclusively")
     }
     valid.one(neg,"logical")
     valid.one(trace,"logical")
})

closedp.internal <- quote({ # calculs principaux des fonctions closed
    n <- sum(na.rm=TRUE,Y)
    nbcap <- rowSums(histpos)      
    nm <- length(lmn) 

    # Initialisation d'objets de la sortie
    tableau<-matrix(nrow=nm,ncol=5)
    dimnames(tableau) <- list(lmn,c("abundance","stderr","deviance","df","AIC"))
    converge<-vector(mode="logical",length=nm)
    names(converge) <- lmn
    glm.<-vector(mode="list",length=nm)
    param<-vector(mode="list",length=nm)
    names(glm.) <- names(param) <- smn
    neg.eta<-vector(mode="list")
    nmC<-1
    
    # Boucle qui ajuste tous les modèles
    for (j in 1:nm)
    {          
        # Construction de la matrice X
        if (smn[j]=="M0")   Xclosedp.out <- Xclosedp(t=t,m="M0",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="Mt")   Xclosedp.out <- Xclosedp(t=t,m="Mt",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MhC")  Xclosedp.out <- Xclosedp(t=t,m="Mh",h="Chao",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MhP")  Xclosedp.out <- Xclosedp(t=t,m="Mh",h="Poisson",theta=2,histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MhD")  Xclosedp.out <- Xclosedp(t=t,m="Mh",h="Darroch",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MhG")  Xclosedp.out <- Xclosedp(t=t,m="Mh",h="Gamma",theta=3.5,histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MthC") Xclosedp.out <- Xclosedp(t=t,m="Mth",h="Chao",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MthP") Xclosedp.out <- Xclosedp(t=t,m="Mth",h="Poisson",theta=2,histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MthD") Xclosedp.out <- Xclosedp(t=t,m="Mth",h="Darroch",histpos=histpos,nbcap=nbcap)  
        if (smn[j]=="MthG") Xclosedp.out <- Xclosedp(t=t,m="Mth",h="Gamma",theta=3.5,histpos=histpos,nbcap=nbcap)  
        if (smn[j]=="Mb")   Xclosedp.out <- Xclosedp(t=t,m="Mb",histpos=histpos,nbcap=nbcap)  
        if (smn[j]=="Mbh")  Xclosedp.out <- Xclosedp(t=t,m="Mbh",histpos=histpos,nbcap=nbcap)
        if (smn[j]=="MhC") {
               mX. <- Xclosedp.out$mat
               mX.names <- Xclosedp.out$paramnames               
        } else {
               mX. <- cbind(Xclosedp.out$mat,mXomit)
               mX.names <- c(Xclosedp.out$paramnames,colnames(mXomit))
        }
        colnames(mX.) <- mX.names
        
        # Ajustement du modèle
        if (trace) {
            iwarn<-getOption("warn")
            options(warn=1)  # Afin que les warnings soient affichés dès qu'ils surviennent
            message(paste("The",lmn[j],"model is fitted."))
        }
        glmo <- glm.call(Y,mX.,cst)
        ##### Réajustement du modèle en enlevant les eta négatifs pour les modèles de Chao
        if (smn[j]%in%c("MhC","MthC")&&neg) {
            nca <- if (smn[j]=="MhC") 2 else t+1
            glmo <- Chao.neg(glmo,nca)
            # Pour identifier les eta fixés à zéro
            pnames <- if(is.null(glmo$offset)) colnames(glmo$model)[-1] else colnames(glmo$model)[-(1:2)]  
            idnegeta<-setdiff(mX.names,pnames)
            if (length(idnegeta)) {
                neg.eta[[nmC]] <- idnegeta 
                names(neg.eta)[nmC] <- smn[j]
                nmC<-nmC+1
            } 
        }
        if (trace) options(warn=iwarn)  # Afin de remettre l'option warn à sa valeur initiale
        glm.[[j]] <- glmo
        converge[j] <- glmo$converge
        
        
        # Calcul de N et de son erreur type
        varcov <- summary(glmo)$cov.unscaled 
        if (smn[j]=="Mb") {
            N <-(exp(glmo$coef[1])*(1+exp(glmo$coef[3]))^t)/(1+exp(glmo$coef[3])-exp(glmo$coef[2]))
            v1 <- 1
            v2 <- exp(glmo$coef[2])/(1+exp(glmo$coef[3])-exp(glmo$coef[2]))
            v3 <-(t*exp(glmo$coef[3])/(1+exp(glmo$coef[3])) -exp(glmo$coef[3])/(1+exp(glmo$coef[3])-exp(glmo$coef[2]))) 
            v <- N*c(v1,v2,v3)
            erreurtype <-sqrt((t(v)%*%varcov%*%v)-N) # calcul de l'erreur type 
        } else
        if (smn[j]=="Mbh") {
            N <- exp(glmo$coef[1])*((1+exp(glmo$coef[4]))^(t-1))*(1+exp(glmo$coef[2]+glmo$coef[3])/(1+exp(glmo$coef[4])-exp(glmo$coef[3])))  # calcul de la taille de la population N
            v1 <- (1+exp(glmo$coef[2]+glmo$coef[3])/(1+exp(glmo$coef[4])-exp(glmo$coef[3])))
            v2 <- exp(glmo$coef[2]+glmo$coef[3])/(1+exp(glmo$coef[4])-exp(glmo$coef[3]))
            v3 <- exp(glmo$coef[2]+glmo$coef[3])/(1+exp(glmo$coef[4])-exp(glmo$coef[3])) + exp(glmo$coef[2]+2*glmo$coef[3])/((1+exp(glmo$coef[4])-exp(glmo$coef[3]))^2)
            v4 <- (t-1)*exp(glmo$coef[4])/(1+exp(glmo$coef[4])) - exp(glmo$coef[2]+glmo$coef[3]+glmo$coef[4])/((1+exp(glmo$coef[4])-exp(glmo$coef[3]))^2)
            v <- exp(glmo$coef[1])*((1+exp(glmo$coef[4]))^(t-1))*c(v1,v2,v3,v4)
            erreurtype <- sqrt((t(v)%*%varcov%*%v) - N)  # calcul de l erreur type        
        } else { 
            N <- n + exp(glmo$coef[1])
            erreurtype <- sqrt(exp(glmo$coef[1])+(exp(2*glmo$coef[1]))*varcov[1,1]) 
        }       
        tableau[j,] <- c(N,erreurtype,glmo$dev,glmo$df.residual,glmo$aic)

        
        # Calcul des probabilités de capture
        fi<-tapply(Y,list(nbcap),sum)
        if (smn[j]=="M0") {
            param[[j]]<-c(N,exp(glmo$coef[2])/(1+exp(glmo$coef[2])))
            nomparam <- c("N","p")
        }
        if (smn[j]=="Mt") {
            param[[j]]<-c(N,exp(glmo$coef[2:(t+1)])/(1+exp(glmo$coef[2:(t+1)])))
            nomparam <- c("N",paste("p",1:t,sep=""))
        }
        if (smn[j]%in%c("MhC","MhP","MhD","MhG")) {
#            ifirstcap <- rep(1:t,2^(t-1:t))
#            param[[j]]<-c(N,sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==1])/N)
            param[[j]]<-c(N,sum(fi*1:t)/(t*N))
            nomparam <- c("N","p")
        }
        if (smn[j]%in%c("MthC","MthP","MthD","MthG")) {
            ifirstcap <- rep(1:t,2^(t-1:t))
            upred <- rep(0,t)
            for ( i in 1:t ) { upred[i] <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==i]) }
            deno <- N
            for ( i in 2:t ) { deno <- c(deno,N-sum(na.rm=TRUE,upred[1:(i-1)])) }          
            param[[j]]<-c(N,upred/deno)
            nomparam <- c("N",paste("p",1:t,sep=""))
        }
        if (smn[j]=="Mb") {
            param[[j]]<-c(N,1-exp(glmo$coef[2])/(1+exp(glmo$coef[3])),exp(glmo$coef[3])/(1+exp(glmo$coef[3])))
            nomparam <- c("N","p","c")
        }
        if (smn[j]=="Mbh") {
            ifirstcap <- rep(1:t,2^(t-1:t))
            p <- sum(na.rm=TRUE,glmo$fitted.values[ifirstcap==1])/N
            cpar <- (glmo$fitted.values[1]/(N*p))^(1/(t-1))
            param[[j]]<-c(N,p,cpar)
            nomparam <- c("N","p","c")
        }
        param[[j]] <- matrix(param[[j]],nrow=1)  # Meilleur format de sortie, pas de notation scientifique
        dimnames(param[[j]]) <- list("estimate:",nomparam)  

    }    

    # Préparation des sorties
    rownames(tableau)<-paste(lmn,ifelse(converge,"","**"),sep="")
    ans <- list(n=n,t=t,results=tableau,converge=converge,glm=glm.,parameters=param,neg.eta=neg.eta,X=X,dfreq=dfreq)
})


#####################################################################################################
# Méthodes pour objets de type closedp et closedp.0

"print.closedp" <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        
        if (!is.null(x$results)) {
          cat("Abundance estimations and model fits:\n")
          tableau <- x$results
          tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
          if (tableau[1,3]>0.001) tableau[,3] <- round(tableau[,3],3)       
          tableau[,4] <- round(tableau[,4],0)
          tableau[,5] <- round(tableau[,5],3)       
          print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
        
          for ( i in 1:length(x$converge))
          {
              if (!x$converge[i]) cat(paste("\n ** : The",names(x$converge[i]),"model did not converge"))
          }
          if (length(x$neg.eta$MhC)==1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameter has been set to zero in the Mh Chao model")
          if (length(x$neg.eta$MhC)>1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameters has been set to zero in the Mh Chao model")
          if (length(x$neg.eta$MthC)==1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameter has been set to zero in the Mth Chao model")
          if (length(x$neg.eta$MthC)>1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameters has been set to zero in the Mth Chao model")
        }
        
        if (dim(tableau)[1]==3) cat("\nNote: When there is 2 capture occasions, only models M0, Mt and Mb are fitted")
        if (dim(tableau)[1]==1) cat("\nNote: When there is 2 capture occasions, only model M0 is fitted")

        cat("\n\n")
        invisible(x)
}

"boxplot.closedp" <- function(x,main="Boxplots of Pearson Residuals", ...){
        if (!is.null(x$results)) {
            model<-which(x$converge)   
            nmodel <- length(model)
            liste <- vector("list",length=nmodel)
            names(liste) <- names(x$glm)[model]
            pres<-function(x){(x$y-fitted(x))/sqrt(fitted(x))}
            for (i in 1:nmodel) liste[[i]] <- pres(x$glm[[model[i]]])
            cex.axis <- if (nmodel>10) 0.75 else 1
            boxplot.default(liste, main=main, cex.axis=cex.axis, ...)
            abline(h=0,lty=3)
        } else cat("Note: There is no residuals to plot\n")
}

"plot.closedp" <- function(x,main="Residual plots for some heterogeneity models", ...){
     typet <- if(any(class(x)=="closedp.t")) TRUE else FALSE
     t <- if(typet) x$t else x$t0
     converge <- x$converge[c("Mh Poisson2","Mh Darroch","Mh Gamma3.5")]
     if (sum(converge)==0) stop("models did not converged, there is no data to plot")  
     pres<-function(x,typet) {
          if (typet) {
               nbcap <- rowSums(histpos.t(t))
               fi <- tapply(x$y,nbcap,sum,na.rm=TRUE)
               fipred <- tapply(fitted(x),nbcap,sum,na.rm=TRUE)
          } else {
               fi <- rev(x$y)
               fipred <- rev(fitted(x))
          }
          (fi-fipred)/sqrt(fipred)
     }
     plotres<-function(res,main) {
          plot(1:t,res[1:t],type="b",ann=FALSE)
          mtext(main,side=3,line=0.5,adj=0,font=2)
          abline(h=0,lty=2)
     }
     op <- par(mfrow=c(sum(converge),1),mar=c(2.1, 3.1, 4.1, 2.1),oma=c(3,2,3,0))
     on.exit(par(op))
     if(converge[1]) plotres(pres(x$glm$MhP,typet),main="Mh Poisson2")
     if(converge[2]) plotres(pres(x$glm$MhD,typet),main="Mh Darroch")
     if(converge[3]) plotres(pres(x$glm$MhG,typet),main="Mh Gamma3.5")
     mtext(main,side=3,cex=1.8,outer=TRUE)
     ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
     mtext(ylab,side=2,outer=TRUE)
     mtext("number of captures",side=1,line=1.5,outer=TRUE)
}

