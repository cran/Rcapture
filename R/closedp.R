"closedp" <- "closedp.t" <- function(X,dfreq=FALSE,neg=TRUE)
{    
  call <- match.call()
  closedp.internal(X=X, dfreq=dfreq, neg=neg, call=call)  
}

"closedp.0" <- function(X,dfreq=FALSE,dtype=c("hist","nbcap"),t=NULL,t0=NULL,neg=TRUE)
{    
  call <- match.call()
  closedp.internal(X=X, dfreq=dfreq, dtype=dtype[1], t=t, t0=t0, neg=neg, call=call)  
}

closedp.internal <- function(X, dfreq=FALSE, dtype="hist", t=NULL, t0=NULL, neg=TRUE, call) 
{
  # Initialisation de variables
  typet <- substr(paste(call[1]), nchar(paste(call[1])), nchar(paste(call[1]))) %in% c("t", "p")
           # différent des autres fonctions car closedp existe (=closedp.t) 
  tinf <- if(is.null(t)) FALSE else is.infinite(t)
  
  ######### Validation des arguments en entrée #########
  valid.one(dfreq,"logical")
  valid.dtype(dtype)
  valid.t(t=t, pInf=!typet)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t, warn=typet)
    X <- Xvalid$X
    t <- Xvalid$t  ## t est modifié s'il prennait la valeur NULL ou Inf
  t0 <- valid.t0(t0=t0, typet=typet, t=t) # doit être soumis après valid.X qui modifie t
  valid.one(neg,"logical")
  ########### Fin de la validation des arguments ###########
  
  
  #### Préparation pour l'ajustement du modèle

  # Création du vecteur de variable réponse Y
  getY.out <- getY(typet=typet, X=X, dfreq=dfreq, dtype=dtype, t=t, t0=t0) 
  Y <- getY.out$Y
  n <- getY.out$n 
  # Préliminaire à la création de la matrice X 
  histpos <- gethistpos(typet=typet, t=t, t0=t0)
  nbcap <- getnbcap(histpos)
  # Création de la variable offset
  cst <- getcst(typet=typet, tinf=tinf, t=t, t0=t0, nbcap=nbcap)
 
  # Identification des modèles à ajuster
  if (typet) {
    if (t==2) {
      lmn <- smn <- c("M0","Mt","Mb")
    } else {
      lmn <- c("M0","Mt","Mh Chao (LB)","Mh Poisson2","Mh Darroch","Mh Gamma3.5","Mth Chao (LB)","Mth Poisson2","Mth Darroch","Mth Gamma3.5","Mb","Mbh")
      smn <- c("M0","Mt","MhC","MhP","MhD","MhG","MthC","MthP","MthD","MthG","Mb","Mbh")
    }
  } else {
    if (t0==2) {
      lmn <- smn <- c("M0")
    } else {
      lmn <- c("M0","Mh Chao (LB)","Mh Poisson2","Mh Darroch","Mh Gamma3.5")
      smn <- c("M0","MhC","MhP","MhD","MhG")
    }
    
  }
  nm <- length(lmn) 
  
  # Initialisation d'objets de la sortie
  tableau <- matrix(nrow=nm,ncol=5)
  dimnames(tableau) <- list(lmn,c("abundance","stderr","deviance","df","AIC"))
  converge <- vector(mode="logical",length=nm)
  names(converge) <- lmn
  glm. <- vector(mode="list",length=nm)
  glm.warn <- vector(mode="list",length=nm)
  param <- vector(mode="list",length=nm)
  names(glm.) <- names(glm.warn) <- names(param) <- smn
  neg.eta <- vector(mode="list")
  nmC <- 1
  
  # Boucle qui ajuste tous les modèles
  t. <- if (typet) t else t0
  for (j in 1:nm)
  {          
    # Construction de la matrice X
    if (smn[j] %in% c("M0", "Mt", "Mb", "Mbh")) { 
      m <- smn[j]
      h <- NULL; theta <- NULL
    } else {
      m <- substr(smn[j],1,nchar(smn[j])-1)
      if (smn[j] %in% c("MhC", "MthC")) {h <- "Chao";    theta <- NULL}
      if (smn[j] %in% c("MhP", "MthP")) {h <- "Poisson"; theta <- 2}
      if (smn[j] %in% c("MhD", "MthD")) {h <- "Darroch"; theta <- NULL}
      if (smn[j] %in% c("MhG", "MthG")) {h <- "Gamma";   theta <- 3.5}
    }  
    Xclosedp.out <- Xclosedp(t=t., m=m, h=h, theta=theta, histpos=histpos, nbcap=nbcap)
    mX. <- Xclosedp.out$mat
    mX.names <- Xclosedp.out$paramnames
    colnames(mX.) <- mX.names
    
    # Ajustement du modèle
    glm.out <- glm.call(Y=Y, mX.=mX., cst=cst, mname=lmn[j])
    glmo <- glm.out$glmo
    glmo.warn <- glm.out$warnings
    
    ##### Réajustement du modèle en enlevant les eta négatifs pour les modèles de Chao
    if (smn[j] %in% c("MhC","MthC") && neg) {
      nca <- if (smn[j] == "MhC") 2 else t+1
      neg.out <- Chao.neg(glmo=glmo, nca=nca, mname=lmn[j])
      glmo <- neg.out$glmo
      glmo.warn <- neg.out$warnings
      # Pour identifier les eta fixés à zéro
      pnames <- if(is.null(glmo$offset)) colnames(glmo$model)[-1] else colnames(glmo$model)[-(1:2)]  
      idnegeta<-setdiff(mX.names,pnames)
      if (length(idnegeta)) {
        neg.eta[[nmC]] <- idnegeta 
        names(neg.eta)[nmC] <- smn[j]
        nmC<-nmC+1
      } 
    }
    glm.[[j]] <- glmo
    if (!is.null(glmo.warn)) glm.warn[[j]] <- glmo.warn  # car affecter NULL à un élément d'une liste efface cet élément, je ne veux pas ça
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
      fi<-tapply(Y,list(nbcap),sum)
      param[[j]]<-c(N,sum(fi*1:t.)/(t.*N))
      ##### À vérifier avec Louis-Paul
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
  ### 22 mai 2012 : On a décidé de ne pas afficher cet avertissement.
  #if (!all(sapply(glm.warn,is.null))) warning("the 'glm.fit' function generated one or more warnings (see the output value 'glm.warn')")
  ###
  rownames(tableau)<-paste(lmn,ifelse(converge,"","**"),sep="")
  ans <- list(n=n, t=t, t0=t0, results=tableau, converge=converge, glm=glm., glm.warn=glm.warn, parameters=param,
              neg.eta=neg.eta, X=X, dfreq=dfreq)
  if (typet) ans$t0 <- NULL
  class(ans) <- if (typet) c("closedp", "closedp.t") else "closedp"
  return(ans)  
}


#####################################################################################################
# Méthodes pour objets de type closedp

"print.closedp" <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")
  
  if (!is.null(x$results)) {
    cat("Abundance estimations and model fits:\n")
    tableau <- x$results
    tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
    if (!is.na(tableau[1,3])) if (tableau[1,3]>0.001) tableau[,3] <- round(tableau[,3],3)       
    tableau[,4] <- round(tableau[,4],0)
    tableau[,5] <- round(tableau[,5],3)       
    print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
    
    for ( i in 1:length(x$converge))
    {
      if (!x$converge[i]) cat(paste("\n ** model",names(x$converge[i]),"did not converge\n"))
    }
    ###################################################
    ### 22 mai 2012 : On a décidé de ne plus imprimer ces notes car l'utilisateur ne comprend pas quel
    ### impact des parametres eta fixés à zéro ont sur ses résultats. Ça l'embête plus qu'autre chose.
    #if (length(x$neg.eta$MhC)==1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameter has been set to zero in the Mh Chao model")
    #if (length(x$neg.eta$MhC)>1) cat("\nNote:",length(x$neg.eta$MhC),"eta parameters have been set to zero in the Mh Chao model")
    #if (length(x$neg.eta$MthC)==1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameter has been set to zero in the Mth Chao model")
    #if (length(x$neg.eta$MthC)>1) cat("\nNote:",length(x$neg.eta$MthC),"eta parameters have been set to zero in the Mth Chao model")
    ###################################################
  }
  
  if (dim(tableau)[1]==3) cat("\nNote: When there is 2 capture occasions, only models M0, Mt and Mb are fitted.\n")
  if (dim(tableau)[1]==1) cat("\nNote: When there is 2 capture occasions, only model M0 is fitted.\n")
  
  cat("\n")
  invisible(x)
}

"boxplot.closedp" <- function(x,main="Boxplots of Pearson Residuals", ...){
  if (!is.null(x$results)) {
    model<-which(x$converge)   
    nmodel <- length(model)
    liste <- vector("list",length=nmodel)
    names(liste) <- names(x$glm)[model]
    pres2<-function(x){(x$y-fitted(x))/sqrt(fitted(x))}
    for (i in 1:nmodel) liste[[i]] <- pres2(x$glm[[model[i]]])
    cex.axis <- if (nmodel>10) 0.75 else 1
    boxplot.default(liste, main=main, cex.axis=cex.axis, ...)
    abline(h=0,lty=3)
  } else cat("Note: There is no residuals to plot.\n")
}

"plot.closedp" <- function(x,main="Residual plots for some heterogeneity models", ...){
  typet <- if(any(class(x)=="closedp.t")) TRUE else FALSE
  t <- if(typet) x$t else x$t0
  converge <- x$converge[c("Mh Poisson2","Mh Darroch","Mh Gamma3.5")]
  if (sum(converge)==0) stop("models did not converge, there is no data to plot")  
  plotres<-function(res,main) {
    plot(1:t,res[1:t],type="b",ann=FALSE,...)
    mtext(main,side=3,line=0.5,adj=0,font=2)
    abline(h=0,lty=2)
  }
  op <- par(mfrow=c(sum(converge),1),mar=c(2.1, 3.1, 4.1, 2.1),oma=c(3,2,3,0))
  on.exit(par(op))
  if(converge[1]) plotres(pres(x$glm$MhP,typet,t),main="Mh Poisson2")
  if(converge[2]) plotres(pres(x$glm$MhD,typet,t),main="Mh Darroch")
  if(converge[3]) plotres(pres(x$glm$MhG,typet,t),main="Mh Gamma3.5")
  mtext(main,side=3,cex=1.8,outer=TRUE)
  ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
  mtext(ylab,side=2,outer=TRUE)
  mtext("number of captures",side=1,line=1.5,outer=TRUE)
}

