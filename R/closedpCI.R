"closedpCI.t" <- function(X, dfreq=FALSE, m=c("M0","Mt","Mh","Mth"), mX=NULL, 
    h=NULL, h.control=list(), mname=NULL, alpha=0.05, fmaxSupCL=3, ...)
{
  call <- match.call()
  closedpCI.internal(X=X, dfreq=dfreq, m=m, mX=mX, h=h, h.control=h.control, mname=mname, 
      alpha=alpha, fmaxSupCL=fmaxSupCL, call=call, ...)  
}

"closedpCI.0" <- function(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=NULL, m=c("M0","Mh"), 
    mX=NULL, h=NULL, h.control=list(), mname=NULL, alpha=0.05, fmaxSupCL=3, ...)
{
  call <- match.call()
  closedpCI.internal(X=X, dfreq=dfreq, dtype=dtype[1], t=t, t0=t0, m=m, mX=mX, h=h, 
      h.control=h.control, mname=mname, alpha=alpha, fmaxSupCL=fmaxSupCL, call=call, ...)  
}

closedpCI.internal <- function(X, dfreq=FALSE, dtype="hist", t=NULL, t0=NULL, m="M0", 
    mX=NULL, h=NULL, h.control=list(), mname=NULL, alpha=0.05, fmaxSupCL=3, call, ...)
{  
  ### Initialisation de variables
  typet <- substr(paste(call[1]), nchar(paste(call[1])), nchar(paste(call[1]))) == "t"
  tinf <- if(is.null(t)) FALSE else is.infinite(t)
  if(!is.list(h.control)) stop("'h.control' must be a list")
  theta <- h.control$theta
  neg <- h.control$neg
  initcoef <- h.control$initcoef
  initsig <- h.control$initsig
  method <- h.control$method
  neg.eta <- NULL

  ### Erreur si un argument est inutilisé
  ### (Cette erreur n'est pas générée automatiquement à cause de l'argument ...)
  nargs <- names(call)  
  nargs.accept <- c("", "X", "dfreq", "m", "mX", "h", "h.control", "mname", "alpha", 
                    "fmaxSupCL", "method", "lower", "upper", "control", "hessian")
  if (!typet) nargs.accept <- nargs.accept <- c(nargs.accept, "dtype", "t", "t0")
  # Il y a trois arguments acceptés de plus avec closedpCI.0
  if (any(!(nargs %in% nargs.accept))) {
    # Pour préparer la liste des arguments inutilisés du message d'erreur. 
    unaccept <- nargs[!(nargs %in% nargs.accept)]
    unaccept.value <- vector(length=length(unaccept))
    for (i in seq_along(unaccept)) {
      unaccept.value[i] <- deparse(call[[unaccept[i]]])
    }
    unaccept.msg.vct <- paste(unaccept, "=", unaccept.value)
    unaccept.msg <- paste(unaccept.msg.vct, collapse=", ")
    stop("unused argument(s) (", unaccept.msg, ")")
  }
  
  ######### Validation des arguments en entrée #########
  valid.one(dfreq,"logical")
  valid.dtype(dtype)
  valid.t(t=t, pInf=!typet)
  Xvalid <- valid.X(X=X, dfreq=dfreq, dtype=dtype, t=t, warn=typet)
  X <- Xvalid$X
  t <- Xvalid$t  ## t est modifié s'il prennait la valeur NULL ou Inf
  t0 <- valid.t0(t0=t0, typet=typet, t=t, tinf=tinf) # doit être soumis après valid.X qui modifie t
  mX <- valid.mX(mX=mX, typet=typet, t=t, t0=t0)  
  if(is.null(mX)) {
    m <- valid.vm(vm=m, values=c("M0","Mt","Mh","Mth"), vt=t, typet=typet)
  } else {
    m <- NULL
    if (!is.null(call[["m"]])) warning("the argument (m = ", call[["m"]], ") was not used since an argument 'mX' was given")
  }
  valid.h.out <- valid.h(h=h, values=c("Chao","LB","Poisson","Darroch","Gamma","Normal"), m=m, call=call)
  h <- valid.h.out$h
  hname <- valid.h.out$hname
  
  if (is.null(theta)) { theta <- if (hname=="Poisson") 2 else if (hname=="Gamma") 3.5 } else valid.theta(theta) 
  if (is.null(neg)) { neg <- if (hname %in% c("Chao", "LB")) TRUE } else valid.one(neg,"logical")
  if (is.null(initsig) || hname=="Normal") initsig <- 0.2 else valid.initsig(initsig)
  if (is.null(method) || hname=="Normal") method <- "BFGS"  # pas de validation car optim va le faire 

  mname <- valid.mname(mname=mname, typet=typet, m=m, hname=hname, theta=theta, call=call)
  valid.alpha(alpha)
  valid.fmaxSupCL(fmaxSupCL)
  
  ########### Fin de la validation des arguments ###########
  
  
  #### Préparation pour l'ajustement du modèle
  
  ### Création du vecteur de variable réponse Y
  getY.out <- getY(typet=typet, X=X, dfreq=dfreq, dtype=dtype, t=t, t0=t0) 
  Y <- getY.out$Y
  n <- getY.out$n 
  
  ### Création de la matrice X
  getmX.out <- getmX(typet=typet, t=t, t0=t0, m=m, h=h, theta=theta, mX=mX)
  mX. <- getmX.out$mX. 
  nbcap <- getmX.out$nbcap
  
  ### Création de la variable offset
  cst <- getcst(typet=typet, tinf=tinf, t=t, t0=t0, nbcap=nbcap)

  ### Indicateur pour le modèle normal
  # (Je n'utilise pas hname=="Normal" pour créer cet indicateur car
  #  il se pourrait que h soit une fonction nommée Normal.)
  indicNormal <- !is.null(h) && !is.function(h) && h == "Normal"
  
  ### Arrêt de la fonction si le modèle a plus de paramètres que le nombre
  ### d'observations pour ajuster le modèle
  nrows <- if(typet) 2^t-1 else t0
  nparams <- if (indicNormal) ncol(mX.) + 2 else ncol(mX.) + 1 
  if (nparams > nrows){ 
    stop("some model's parameters are not estimable since the number of parameters (", 
         nparams, ")\n  is higher than the number of observations to fit the model (", 
         nrows, ")", 
         ## Si c'est à cause du t0 donné j'ajoute : 
         ngettext(!typet && (t0 < t) && (nparams <= t), 
         ",\n  the given t0 argument is too small for the requested model", "") )
  }
  
  if (indicNormal) {  ## Si un modèle hétérogène normal a été demandé :
    
    # Ajustement du modèle loglinéaire avec matrice X : mX. + colonne Darroch
    # Utile pour déterminer les valeurs initiales des paramètres loglinéaire si 'initcoef' n'est pas fourni
    # Aussi utile pour déterminer le rang de la matrice X dans le calcul du AIC  
    mX.D <- cbind(mX.,tau=(nbcap^2)/2)
    glmD <- glm.call(Y=Y, mX.=mX.D, cst=cst, 
                     mname=paste(mname, "+ h=Darroch to get initial values for the loglinear parameters"))$glmo
    
    
    #### Validation des arguments en entrée partie 2 ####
    # Argument initcoef
    if(!is.null(initcoef)) {
      if(length(initcoef)!=ncol(mX.)+1) stop("'initcoef' must be of length 'ncol(mX)+1'")
      if(!is.numeric(initcoef)) stop("'initcoef' must be numeric")
    }
    #######################################################
    
    
    #### Détermination des valeurs initiales des paramètres dans optim
    if(is.null(initcoef)) { # si non fournies : mX. + colonne Darroch
      if (!glmD$converge) stop("the Darroch's model, fitted to find initial values for the loglinear coefficients, did not converge : please input argument 'initcoef'")
      initparamll <- glmD$coef[1:(ncol(mX.)+1)]
    } else initparamll <- initcoef
    names(initparamll) <- c("(Intercept)", colnames(mX.))
    initparam <- c(initparamll, sigma=initsig)
    
    ###  Ajustement du modèle normal avec optim
    Max <- optim(initparam, devPois, Y=Y, mX=mX., nbcap=nbcap, cst=cst, method=method, hessian=TRUE, ...)
    if(Max$par[ncol(mX.)+2]<=0) warning("The sigma estimate is non-positive. It seems there is not heterogeneous catchability in the data.")
    
    
    ######### Calculs pour préparation de la sortie
    N <- n + exp(Max$par[1])
    varcov <- solve(Max$hessian)
    erreurtype <- sqrt(exp(Max$par[1])+(exp(2*Max$par[1]))*varcov[1,1])
    # Chao 1987 formule 12
    qZ <- qnorm(1-alpha/2)
    C <- exp(qZ*sqrt(log(1+erreurtype^2/(N-n)^2))) 
    InfCL <- n + (N-n)/C
    SupCL <- n + (N-n)*C
    dev <- Max$value
    dff <- nrow(mX.) - ncol(mX.) - 2
    param <- rbind(estimate=Max$par,stderr=diag(varcov))
    mu <- getmu(param=Max$par,Y,mX.,nbcap,cst)
    aic <- -2*sum(dpois(Y, mu, log=TRUE)) + 2*glmD$rank
    # le premier élément de la somme est dérivé de family.R ligne 170
    # le deuxième élément de la somme est dérivé de glm.R ligne 369
    # on utilise la même valeur du rang pour les modèles Darroch et normal
    results <- matrix(c(N,erreurtype,InfCL,SupCL,dev,dff,aic),nrow=1)
    mname <- paste(mname,ifelse(Max$converge==0,"","**"),sep="")
    dimnames(results) <- list(mname,c("abundance","stderr","InfCL","SupCL","deviance","df","AIC"))
    
    fit <- list(parameters=param,varcov=varcov,y=Y,fitted.values=mu,initparam=initparam,optim.out=Max)
    
  } else { ## Si on n'a pas demandé un modèle hétérogène normal : 
    
    # Indicateur hétéro Chao (ou LB)
    hChao <- if(is.function(h) || is.null(h)) { FALSE } else { if(hname %in% c("Chao", "LB")) TRUE else FALSE}  
    
    ######### Ajustement du modèle avec un modèle loglinéaire poisson -> Obtention de l'estimateur de N (sans mu_0).
    glm.out <- glm.call(Y=Y, mX.=mX., cst=cst, mname=mname)
    glmo <- glm.out$glmo
    glmo.warn <- glm.out$warnings
    ##### Réajustement du modèle en enlevant les eta négatifs pour les modèles de Chao
    if (hChao && neg) {
      nca <- getmX.out$nca
      neg.out <- Chao.neg(glmo=glmo, nca=nca, mname=mname)
      glmo <- neg.out$glmo
      glmo.warn <- neg.out$warnings
      neg.eta <- setdiff(colnames(mX.),colnames(glmo$model[,-1]))
    }
    N <- n+exp(glmo$coef[1])
    varcov <- summary(glmo)$cov.unscaled
    erreurtype <- sqrt(exp(glmo$coef[1])+(exp(2*glmo$coef[1]))*varcov[1,1])
    results <- matrix(c(N,erreurtype,glmo$dev,glmo$df.residual,glmo$aic),nrow=1)
    
    ######### Fonction de calcul de la log vraisemblance multinomiale profile à optimiser.
    mXavec <- rbind(mX.,rep(0,ncol(mX.)),deparse.level=0)
    cstavec <- c(cst,0)
    Nval <- NULL
    loglikval <- NULL
    loglikemult <- function(N,lobj=0)
    {
      n0 <- N-n
      Yavec <- c(Y,n0)
      glmoavec <- glm.call(Y=Yavec, mX.=mXavec, cst=cstavec, mname=paste(mname, "for the profile CI computation"))$glmo
      if (is.null(mX) && m %in% c("Mh", "Mth") && hChao && neg) {
        nca <- if (m=="Mh") 2 else t+1
        glmoavec <- Chao.neg(glmo=glmoavec, nca=nca, mname=paste(mname, "for the profile CI computation"))$glmo
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
    # En théorie N (N Poisson) > Nmax (N multinomial), on pourrait donc faire la recherche 
    # sur l'intervalle (n, N). Mais si N est très proche de Nmax ça pourrait peut-être causer de problèmes.
    # On prend donc l'intervalle de recherche (n, N) qui contient assurémant le maximum.
    Nmax <- opmax$maximum
    lmax <- opmax$objective
    lminCI <- lmax-qchisq(1-alpha,1)/2
    
    ######### Détermination de la borne inférieure
    infroot <- try(uniroot(loglikemult, c(n, Nmax), lobj=lminCI, tol = 0.0001),silent=TRUE)
    InfCL <- if(!inherits(infroot, "try-error")) infroot$root else n
    
    ######### Détermination de la borne supérieure
    suproot <- try(uniroot(loglikemult, c(Nmax, fmaxSupCL*N), lobj=lminCI, tol = 0.0001),silent=TRUE)
    SupCL <- if(!inherits(suproot, "try-error")) suproot$root else paste(">",round(fmaxSupCL*N,1),sep="")
    
    ######### Préparation de la sortie
    loglikval <- loglikval[order(Nval)]
    Nval <- Nval[order(Nval)]
    if(!inherits(infroot, "try-error")) { # Si InfCI > n, je vais tronquer mes vecteurs pour les graphiques
      posInf <- which(Nval < InfCL - Nmax*0.02)
      gInf <- if(length(posInf)==0) 1 else max(posInf)
      loglikval <- loglikval[gInf:length(loglikval)]
      Nval <- Nval[gInf:length(Nval)]
    }                     
    if(!inherits(suproot, "try-error")) { # Si SupCL < fmaxSupCL*N, je vais tronquer mes vecteurs pour les graphiques
      posSup <- which(Nval > SupCL + Nmax*0.02)
      gSup <- if(length(posSup)==0) length(Nval) else  min(posSup)
      loglikval <- loglikval[1:gSup]
      Nval <- Nval[1:gSup]
    }                     
    CI <- if(!inherits(suproot, "try-error")) matrix(c(Nmax,InfCL,SupCL),nrow=1) else data.frame(Nmax,InfCL,SupCL)
    if (!glmo$converge) mname <- paste(mname,"**",sep="") 
    dimnames(CI) <- list(mname,c("abundance","InfCL","SupCL"))
    dimnames(results) <- list(mname,c("abundance","stderr","deviance","df","AIC"))
    
    fit <- glmo
    
  }
  
  ans <- list(n=n,t=t,t0=t0,results=results,fit=fit)
  if (!indicNormal) {
    ### 22 mai 2012 : On a décidé de ne pas afficher cet avertissement.
    #if (!is.null(glmo.warn)) warning("the 'glm.fit' function generated one or more warnings (see the output value 'glm.warn')")
    ###
    ans <- c(ans, list(glm.warn=glmo.warn))
    if (hChao) ans <- c(ans, list(neg.eta=neg.eta))
    ans <- c(ans, list(CI=CI,alpha=alpha,N.CI=Nval,loglik.CI=loglikval)) 
  } else ans <- c(ans, list(alpha=alpha))
  if (typet) ans$t0 <- NULL
  class(ans) <- if (typet) c("closedpCI", "closedpCI.t") else "closedpCI"
  return(ans)  
}


getmu <- function(param,Y,mX,nbcap,cst){
  #############################################################################################################
# Fonction qui calcule les valeurs prédites par le "modèle normal"
# On doit lui fournir :
# param: la valeur de tous les paramètres du modèle (l'ordonnée à l'origine + un paramètre par colonne de la 
#          matrice mX + sigma, en respectant cet ordre)
# Y: les fréquences observées (respectant le même ordre que mX)
# mX: la matrice X de la partie log-linéaire du modèle (celle donnée en entrée à closedpNormal.t, 
#       sans colonne de 1 pour l'ordonnée à l'origine)
#   nbcap: le nombre de captures par historique de capture (respectant le même ordre que mX)
# cst: l'offset à ajouter au modèle (0 pour fct .t, mais utile pour fct .0)
# Elle retourne un vecteur de longueur 2^t-1: les valeurs prédites par le "modèle normal" (mu)
  #############################################################################################################
  
  mXi <- cbind(rep(1,nrow(mX)),mX)
  mXi_nbcap <- cbind(mXi,nbcap)
  t <- max(nbcap)
  intercept <- param[1]
  paramll <- param[-(ncol(mX)+2)]
  sig <- param[ncol(mX)+2]
  
# Les points et les poids pour l'intégration numérique par la méthode de quadrature de Gauss
# tiré de Aramowitz et Stegun (1972), page 924 (Hermite integration)
# On peut aussi retrouver ces valeurs en ligne http://www.efunda.com/math/num_integration/findgausshermite.cfm
  xnor<-c(-5.38748089001 ,-4.60368244955 ,-3.94476404012 ,-3.34785456738 ,-2.78880605843 ,-2.25497400209 ,
      -1.73853771212 ,-1.2340762154  ,-0.737473728545,-0.245340708301,0.245340708301 ,
      0.737473728545 ,1.2340762154   ,  1.73853771212  ,2.25497400209  ,2.78880605843  ,
      3.34785456738  ,3.94476404012  ,4.60368244955  ,5.38748089001)
  wnor<-c(2.22939364554E-013,4.39934099226E-010,1.08606937077E-007,7.8025564785E-006 ,0.000228338636017 ,
      0.00324377334224 ,0.0248105208875 , 0.10901720602 ,0.286675505363 ,0.462243669601 ,  
      0.462243669601 ,   0.286675505363 ,  0.10901720602 ,  0.0248105208875 , 
      0.00324377334224  ,0.000228338636017 ,7.8025564785E-006 ,1.08606937077E-007,
      4.39934099226E-010,2.22939364554E-013)
  
# Calcul des valeurs de la constante C pour les 20 points de l'intégration numérique (C(sigma*z_i) dans la formule de Rivest (2011), page 7)
# C(alpha)=p(0|alpha) conditionnal probability of not appearing in any list
  C<-vector(length=20)
  for ( i in (1:20)){ C[i]<-exp(intercept)/(exp(intercept)+sum(exp(cst+mXi_nbcap%*%c(paramll,xnor[i]*sqrt(2)*sig)))) }
  
# Calcul de la fonction varphi (formule de Rivest (2011), page 7)
  phi<-vector(length=t)
  for (i in (1:t)){ phi[i]<-log(sum(wnor*exp(i*sqrt(2)*sig*xnor)*C)/sum(wnor*C)) }
  
# Calcul des valeurs prédites
  mu<-exp(cst + mXi%*%paramll + phi[nbcap]) # eq3, Rivest (2011), page 3
  mu <- as.vector(mu)
  names(mu) <- 1:length(mu)
  return(mu)
}


devPois <- function(param,Y,mX,nbcap,cst) {
  #############################################################################################################
# Fonction qui sera donnée en entrée à optim.
# Elle calcule la valeur de la statistique à minimiser (la déviance Poisson) pour le "modèle normal".
# On doit lui fournir :
# param: la valeur de tous les paramètres du modèle (l'ordonnée à l'origine + un paramètre par colonne de la 
#          matrice mX + sigma, en respectant cet ordre)
# Y: les fréquences observées (respectant le même ordre que mX)
# mX: la matrice X de la partie log-linéaire du modèle (celle donnée en entrée à closedpNormal.t, 
#       sans colonne de 1 pour l'ordonnée à l'origine)
#   nbcap: le nombre de captures par historique de capture (respectant le même ordre que mX)
# cst: l'offset à ajouter au modèle (0 pour fct .t, mais utile pour fct .0)
# Elle retourne un seul chiffre: la valeur de la déviance poisson
  #############################################################################################################
  
  mu <- getmu(param,Y,mX,nbcap,cst)
  dev.resids <- function(y, mu) { 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu)) }   #  dérivé de family.R ligne 169
  dev <- sum(dev.resids(y=Y, mu=mu)) # dérivé de glm.R ligne 183
  return(dev)
}




#####################################################################################################
# Méthodes pour objets de type closedpCI et closedpCI.0

"print.closedpCI" <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")  
  if (is.null(x$N.CI)) {
    cat("Abundance estimation,",paste((1-x$alpha)*100,"%",sep=""),"confidence interval and model fit:\n")
    tableau <- x$results
    tableau[,1:4] <- round(tableau[,1:4],1)
    tableau[,6] <- round(tableau[,6],0)
    tableau[,c(5,7)] <- round(tableau[,c(5,7)],3)       
    print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
    if (x$fit$optim.out$converge!=0) cat(paste("\n ** the model did not converge (see the output value 'fit$optim.out' for more details)\n"))  
  } else {    
    cat("Poisson estimation and model fit:\n")
    tableau <- x$results
    tableau[,c(1,2)] <- round(tableau[,c(1,2)],1)
    tableau[,4] <- round(tableau[,4],0)
    tableau[,c(3,5)] <- round(tableau[,c(3,5)],3)       
    print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
    ###################################################
    ### 22 mai 2012 : On a décidé de ne plus imprimer ces notes car l'utilisateur ne comprend pas quel
    ### impact des parametres eta fixés à zéro ont sur ses résultats. Ça l'embête plus qu'autre chose.
    #if (length(x$neg.eta)==1) cat("\nNote:",length(x$neg.eta),"eta parameter has been set to zero\n")
    #if (length(x$neg.eta)>1) cat("\nNote:",length(x$neg.eta),"eta parameters has been set to zero\n")
    ###################################################
    
    cat("\nMultinomial estimation,",paste((1-x$alpha)*100,"%",sep=""),"profile likelihood confidence interval:\n")
    if (is.data.frame(x$CI)) {
      tableau2 <- x$CI
      tableau2[,c(1,2)] <- round(tableau2[,c(1,2)],1)
      print.data.frame(tableau2, ...)
    } else {
      tableau2 <- round(x$CI,1)
      print.default(tableau2, print.gap = 2, quote = FALSE, right=TRUE, ...)
    }  
    if (!x$fit$converge) cat(paste("\n ** the model did not converge\n"))
  }  
  cat("\n")
  invisible(x)
}

"plotCI" <- function(x.closedpCI, main = "Profile Likelihood Confidence Interval", ...) {
  ############################################################################################################################
  # Validation de l'argument fourni en entrée
  if(!any(class(x.closedpCI)=="closedpCI")) stop("'x.closedpCI' must be an object produced with 'closedpCI.t' or 'closedpCI.0")
  ############################################################################################################################
  
  if (is.null(x.closedpCI$N.CI)) {
    message("For normal heterogeneous models, a log-transformed confidence interval is constructed instead of a profile likelihood one.",
        "\nTherefore, 'plotCI' cannot produce a plot of the multinomial profile likelihood for the given 'x.closedpCI'." )
  } else {   
    plot.default(x.closedpCI$N.CI,x.closedpCI$loglik.CI,type="l",ylab="multinomial profile loglikelihood",xlab="N",main=main, ...)
    # Ajout de lignes verticales pour identifier les borne et l'estimation ponctuelle
    lmax <- max(x.closedpCI$loglik.CI); lmin <- min(x.closedpCI$loglik.CI); 
    N <- x.closedpCI$CI[1,1]; InfCL <- x.closedpCI$CI[1,2]; SupCL <- x.closedpCI$CI[1,3]  
    lInf <- if (InfCL==x.closedpCI$n) x.closedpCI$loglik.CI[1] else lmax-qchisq(1-x.closedpCI$alpha,1)/2
    segments(x0=InfCL,y0=lmin,x1=InfCL,y1=lInf)
    text(InfCL,lmin,round(InfCL,2),pos=1,offset=0.2,xpd=NA)
    if (!is.factor(SupCL)) {
      segments(x0=SupCL,y0=lmin,x1=SupCL,y1=lmax-qchisq(1-x.closedpCI$alpha,1)/2)
      text(SupCL,lmin,round(SupCL,2),pos=1,offset=0.2,xpd=NA)
    }
    segments(x0=N,y0=lmin,x1=N,y1=lmax,lty=2)
    text(N,lmin,round(N,2),pos=1,offset=0.2,xpd=NA)
  }
}

"boxplot.closedpCI" <- function(x,main="Boxplots of Pearson Residuals", ...) {
  boxplot.default((x$fit$y-x$fit$fitted.values)/sqrt(x$fit$fitted.values), main=main, ...)     
}

"plot.closedpCI" <- function(x,main="Scatterplot of Pearson Residuals", ...){
  typet <- if(any(class(x)=="closedpCI.t")) TRUE else FALSE
  t <- if (typet) x$t else x$t0
  res <- pres(x=x$fit, typet=typet, t=t) 
  ylab <- if(typet) "Pearson residuals in terms of fi (number of units captured i times)" else "Pearson residuals"
  plot(1:t,res,type="b",main=main,xlab="number of captures",ylab=ylab, ...)
  abline(h=0,lty=2)
}