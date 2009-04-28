"closedpCI.t" <-
function(X, dfreq=FALSE, m=c("M0","Mt","Mh","Mth"), h=c("Chao","Poisson","Darroch","Gamma"), theta=2, mX=NULL, mname, neg=TRUE, alpha=0.05)
{
     eval(valid.t)
     eval(valid.closedpCI)
     eval(model.t)
     eval(closedpCI.internal)
     class(ans) <- c("closedpCI","closedpCI.t")
     return(ans)
}

"closedpCI.0" <-
function(X, dfreq=FALSE, dtype=c("hist","nbcap"), t, t0=t, m=c("M0","Mh"), h=c("Chao","Poisson","Darroch","Gamma"), theta=2, mX=NULL, mname, neg=TRUE, alpha=0.05)
{
     eval(valid.0)
     eval(valid.closedpCI)
     eval(model.0)
     eval(closedpCI.internal)
     ansf <- c(ans,list(t0=t0))
     class(ansf) <- "closedpCI"
     return(ansf)
}

# Bouts de code communs aux deux fonctions : expressions R prêtes à être évaluées

valid.closedpCI <- quote({ # Validation des arguments des fonctions closedpCI
     valid.one(dfreq,"logical")
     Xvalid<-valid.X(X,dfreq,dtype,t,warn=typet)
        X <- Xvalid$X
        t <- Xvalid$t
     if (t<2) stop("the number of capture occasions 't' must be at least 2")    
     if (!typet) {    
        valid.one(t0,"numeric")
        if ((t0%%1)!=0||t0>t||t0<2) 
              stop("'t0' must be an integer between t, the number of capture occasion, and 2 inclusively")
     }
     mpos <- if(typet) c("M0","Mt","Mh","Mth") else c("M0","Mh")
     m<-valid.vm(m,mpos,t)
     hname <- deparse(substitute(h))
     h<-valid.vh(h,c("Chao","Poisson","Darroch","Gamma"),m)[[1]]
     hname <- if(is.function(h)) hname else h
     valid.one(theta,"numeric")
     valid.one(neg,"logical")
     valid.one(alpha,"numeric")
     # Argument mname
     if (missing(mname)) {
       mname <- if(!is.null(mX)) deparse(substitute(mX)) else if (m%in%c("M0","Mt")) m else
                if(hname%in%c("Poisson","Gamma")) paste(m,paste(hname,theta,sep="")) else paste(m,hname)
     } else valid.one(mname,"character")
     # Argument mX
     if(!is.null(mX)){
       mX<-as.matrix(mX)
       if(typet) {
            if (dim(mX)[1]!=2^t-1) stop("'mX' must have 2^t-1 rows")
       } else {
            if (dim(mX)[1]!=t) stop("'mX' must have t rows")
       } 
     }
})

closedpCI.internal <- quote({ # calculs principaux des fonctions closedCI
        n <- sum(na.rm=TRUE,Y)
        if (is.null(mX)) {
            Xclosedp.out <- Xclosedp(t=t,m=m,h=h,theta=theta,histpos=histpos)
            mX. <- Xclosedp.out$mat
            mX.names <- Xclosedp.out$paramnames
        } else {
            mX. <- mX
            mX.names <- colnames(mX.)
            if(is.null(mX.names)||any(mX.names=="")) mX.names <- paste("col",1:dim(mX.)[2],sep="")
        }
        if (!(m=="Mh"&&hname=="Chao")) {
            mX. <- cbind(mX.,mXomit)
            mX.names <- c(mX.names,colnames(mXomit))
        }
        colnames(mX.) <- mX.names
        mXavec <- rbind(mX.,rep(0,dim(mX.)[2]),deparse.level=0)
        cstavec <- c(cst,0)
        neg.eta <- NULL
        Nval <- NULL
        loglikval <- NULL
                      
######### Ajustement du modèle avec un modèle loglinéaire poisson -> Obtention de l'estimateur de N (sans mu_0).
        glmo <- glm.call(Y,mX.,cst)
        ##### Réajustement du modèle en enlevant les eta négatifs pour les modèles de Chao
        if (m%in%c("Mh","Mth")&&hname=="Chao"&&neg) {
            nca <- if (m=="Mh") 2 else t+1
            glmo <- Chao.neg(glmo,nca)
            neg.eta<-setdiff(Xclosedp.out$paramnames,colnames(glmo$model[,-1]))
        }
        N <- n+exp(glmo$coef[1])
        varcov <- summary(glmo)$cov.unscaled
        erreurtype <- sqrt(exp(glmo$coef[1])+(exp(2*glmo$coef[1]))*varcov[1,1])
        results <- matrix(c(N,erreurtype,glmo$dev,glmo$df.residual,glmo$aic),nrow=1)

######### Fonction de calcul de la log vraisemblance multinomiale profile à optimiser.
          loglikemult <- function(N,lobj=0)
          {
               n0 <- N-n
               Yavec <- c(Y,n0)
               glmoavec <- suppressWarnings(glm.call(Yavec,mXavec,cstavec))
                    # On omet les warnings pour ne pas voir plus d'une fois l'avertissement pour algo non convergent
               if (m%in%c("Mh","Mth") && hname=="Chao" && neg) {
                   nca <- if (m=="Mh") 2 else t+1
                   glmoavec <- suppressWarnings(Chao.neg(glmoavec,nca))
               }
                   
               # Calcul du terme correctif (Cormack 1992)
               Nn0 <- sum(na.rm=TRUE,Yavec)
               if(Nn0>100){
                   ct <- if (n0==0||n0==1) -Nn0+0.5*log(2*pi*Nn0) else n0-Nn0-0.5*log(n0/Nn0)
               } else ct <- log((n0^n0)*factorial(Nn0)/((Nn0^Nn0)*factorial(n0))) 
               
               # log vraisemblance multinomiale profile
               loglik <- (glmoavec$deviance - 2*ct)/(-2) - lobj
               loglikval <<- c(loglikval,loglik+lobj)
               Nval <<- c(Nval,N)
               return(loglik)                
          }

######### Détermination du maximum
          opmax <- optimize(loglikemult, c(n, 3*N), tol = 0.0001, maximum=TRUE)
          Nmax <- opmax$maximum
          lmax <- opmax$objective
          lminCI <- lmax-qchisq(1-alpha,1)/2

######### Détermination de la borne inférieure
          infroot <- try(uniroot(loglikemult, c(n, Nmax), lobj=lminCI, tol = 0.0001),silent=TRUE)
          InfCL <- if(!inherits(infroot, "try-error")) infroot$root else n

######### Détermination de la borne supérieure
          suproot <- try(uniroot(loglikemult, c(Nmax, 3*N), lobj=lminCI, tol = 0.0001),silent=TRUE)
          SupCL <- if(!inherits(suproot, "try-error")) suproot$root else paste(">",round(3*N,1),sep="")

######### Préparation de la sortie
        loglikval <- loglikval[order(Nval)]
        Nval <- Nval[order(Nval)]
        if(!inherits(infroot, "try-error")) { # Si InfCI > n, je vais tronquer mes vecteurs pour les graphiques
               posInf <- which(Nval < InfCL - Nmax*0.02)
               gInf <- if(length(posInf)==0) 1 else max(posInf)
               loglikval <- loglikval[gInf:length(loglikval)]
               Nval <- Nval[gInf:length(Nval)]
        }                     
        if(!inherits(suproot, "try-error")) { # Si SupCI < 3*N, je vais tronquer mes vecteurs pour les graphiques
               posSup <- which(Nval > SupCL + Nmax*0.02)
               gSup <- if(length(posSup)==0) length(Nval) else  min(posSup)
               loglikval <- loglikval[1:gSup]
               Nval <- Nval[1:gSup]
        }                     
        CI <- if(!inherits(suproot, "try-error")) matrix(c(Nmax,InfCL,SupCL),nrow=1) else data.frame(Nmax,InfCL,SupCL)
        mname <- paste(mname,ifelse(glmo$converge,"","**"),sep="")
        dimnames(CI) <- list(mname,c("abundance","InfCL","SupCL"))
        dimnames(results) <- list(mname,c("abundance","stderr","deviance","df","AIC"))

        ans <- list(n=n,t=t,results=results,glm=glmo,neg.eta=neg.eta,CI=CI,alpha=alpha,N.CI=Nval,loglik.CI=loglikval)
})


#####################################################################################################
# Méthodes pour objets de type closedpCI et closedpCI.0

"print.closedpCI" <- function(x, ...) {
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Poisson estimation and model fit:\n")
        tableau <- x$results
        tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
        tableau[,4] <- round(tableau[,4],0)
        tableau[,c(3,5)] <- round(tableau[,c(3,5)],3)       
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
        if (length(x$neg.eta)==1) cat("\nNote:",length(x$neg.eta),"eta parameter has been set to zero\n")
        if (length(x$neg.eta)>1) cat("\nNote:",length(x$neg.eta),"eta parameters has been set to zero\n")

        cat("\nMultinomial estimation,",paste((1-x$alpha)*100,"%",sep=""),"profile likelihood confidence interval:\n")
        if (is.data.frame(x$CI)) {
               tableau2 <- x$CI
               tableau2[,c(1,2)] <- round(tableau2[,c(1,2)],1)
               print.data.frame(tableau2, ...)
        } else {
               tableau2 <- round(x$CI,1)
               print.default(tableau2, print.gap = 2, quote = FALSE, right=TRUE, ...)
        }  
        if (!x$glm$converge) cat(paste("\n ** : The model did not converge"))
        cat("\n")
        invisible(x)
}

"plotCI" <- function (x, ...) UseMethod("plotCI")

"plotCI.closedpCI" <- function(x,main="Profile Likelihood Confidence Interval", ...) {
        plot.default(x$N.CI,x$loglik.CI,type="l",ylab="multinomial profile loglikelihood",xlab="N",main=main, ...)
        # Ajout de lignes verticales pour identifier les borne et l'estimation ponctuelle
        lmax <- max(x$loglik.CI); lmin <- min(x$loglik.CI); 
        N <- x$CI[1,1]; InfCL <- x$CI[1,2]; SupCL <- x$CI[1,3]  
        lInf <- if (InfCL==x$n) x$loglik.CI[1] else lmax-qchisq(1-x$alpha,1)/2
        segments(x0=InfCL,y0=lmin,x1=InfCL,y1=lInf)
        text(InfCL,lmin,round(InfCL,2),pos=1,offset=0.2,xpd=NA)
        if (!is.factor(SupCL)) {
               segments(x0=SupCL,y0=lmin,x1=SupCL,y1=lmax-qchisq(1-x$alpha,1)/2)
               text(SupCL,lmin,round(SupCL,2),pos=1,offset=0.2,xpd=NA)
        }
        segments(x0=N,y0=lmin,x1=N,y1=lmax,lty=2)
        text(N,lmin,round(N,2),pos=1,offset=0.2,xpd=NA)
}

"boxplot.closedpCI" <- function(x,main="Boxplots of Pearson Residuals", ...) {
        boxplot.default((x$glm$y-fitted(x$glm))/sqrt(fitted(x$glm)),main=main, ...)     
}

"plot.closedpCI" <- function(x,main="Scatterplot of Pearson Residuals", ...){
     typet <- if(any(class(x)=="closedpCI.t")) TRUE else FALSE
     t <- if(typet||rownames(x$results)=="Mh Chao") x$t else x$t0 
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
     res<-pres(x$glm,typet) 
     ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
     plot(1:t,res[1:t],type="b",main=main,xlab="number of captures",ylab=ylab)
     abline(h=0,lty=2)
}
