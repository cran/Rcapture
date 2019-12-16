#--------------------------------------------------------------------------------------------------#
#### Fonctions pour des bouts de code repetes ####
# se nomme .DRY, qui represente la philosophie de programmation "Don't Repeat Yourself"

getY <- function(typet, X, dfreq, dtype, t, t0) {
#?#  Creation du vecteur de variable reponse Y
  Y <- if (typet) {
        histfreq.t(X=X,dfreq=dfreq)
      } else {
        histfreq.0(X=X,dfreq=dfreq,dtype=dtype[1],vt=t)
      }
  # Valeur derivee de ce calcul : n
  n <- sum(na.rm=TRUE,Y)
  # Modification de Y pour fonction .0 (typet=FALSE) au besoin (t0 different de t).
  if (!typet) Y <- modifYt0(Y=Y, t=t, t0=t0)
  return(list(Y=Y, n=n))
}


modifYt0 <- function(Y, t, t0){
#?#  Modification de Y si t0 est different de t
  if (t0 < t) {
    Y <- Y[(1 + t - t0):t] ## On enleve les premiere frequences du vecteur, soit celle pour nbcap>t0.    
  } else if (t0 > t) { ## cas ou l'argument t=Inf et t0 est superieur a tmax (maintenant = t)
    Y <- c(rep(0, t0 - t), Y) ## on doit ajouter des zeros lorsque plus de t captures
  }
  return(Y)
}


gethistpos <- function(typet, t, t0=t) {
#?#  Obtention de la matrice des historiques de capture possibles
  histpos <- if (typet) {
    histpos.t(t=t)
  } else {
    histpos.0(vt=t0)
  }
  return(histpos)
}


getnbcap <- function(histpos) {
#?#  Calcul du nombre de captures associe a chaque historique de capture possible
  rowSums(histpos)
}


getmX <- function(typet, t, t0, m, h=NULL, theta=NULL, mX=NULL) {
#?# Description de la fonction : Creation de la matrice X
  
  ### Creation preliminaire d'objets necessaires aux calculs   
  histpos <- gethistpos(typet=typet, t=t, t0=t0)
  nbcap <- getnbcap(histpos)
  t. <- if (typet) t else t0
  
  ### Si un mX sous forme de formule a ete donne, on doit d'abord
  ### transformer la formule en la matrice equivalente.
  # (On sait a cause de la validation de mX que si mX est une formule
  #  alors necessairement typet=TRUE, donc t0=NULL) 
  if (inherits(mX, "formula")){   
    dfhistpos <- as.data.frame(histpos)
    colnames(dfhistpos) <- paste("c",1:t,sep="")    
    Terms <- terms(mX, data=dfhistpos)
    mX <- model.matrix(mX, data=dfhistpos)
    # Pour traiter l'ordonnee a l'origine
    if (attr(Terms, "intercept")==0){  ## si l'intercept a ete enleve : erreur
      stop("the intercept must not be removed in the formula given as 'mX' argument", call. = FALSE)
    } else {  ## si l'intercept est la c'est bien, mais on l'enleve il sera remis plus tard
      mX <- mX[,-1, drop = FALSE]  
    }
    # Pour changer les noms des colonnes dans cette matrice : 
    if (ncol(mX) > 0) {
      rsplit <- strsplit(colnames(mX), split="c")
      vrsplit <- sapply(rsplit,paste,collapse="")
      id <- nchar(vrsplit)
      # si un seul chiffre, on paste "beta" devant
      # si plus d'un chiffre, on paste "lambda" devant
      # (Je suis ici la notation de Rivest (2011), mais sans les parantheses
      #  car R pense que lambda est une fonction et la virgule est remplacee par _.)
      pnames <- paste(ifelse(id==1,"beta","lambda"), vrsplit, sep="")
      # puis on remplace le ":" par "_".
      colnames(mX) <- gsub(":", "_", pnames)
    }
  }
  
  ### Ensuite on peut obtenir la matrice de design mX.
  # (mX et mX. ne sont pas necessairement equivalentes puisque des colonnes
  #  peuvent etre ajoutees a mX pour l'heterogeneite.)  
  if (!is.null(mX) && ncol(mX) == 0 && is.null(h)) {
    # Si le modele comprend seulement une ordonnee a l'origine
    mX. <- mX    
  } else {
    Xclosedp.out <- Xclosedp(t=t., m=m, h=h, theta=theta, histpos=histpos, nbcap=nbcap, mX=mX)
    mX. <- Xclosedp.out$mat
    colnames(mX.) <- Xclosedp.out$coeffnames
  }
  
  ### Sortie des resultats
  nca <- if (is.null(m)) ncol(mX) + 1 else if (m %in% c("M0", "Mh")) 2 else t+1   
    ## nca : le nombre de colonnes dans la matrice de design (incluant une colonne
    ## pour l'ordonnee a l'origine) ne modelisant pas l'heterogeneite
    ## utile pour la correction des eta negatifs pour le modele Chao ou LB 
    
  return(list(mX.=mX., nbcap=nbcap, nca = nca))
}


getcst <- function(typet, tinf, t, t0, nbcap){
#?# Creation de la variable offset
  cst <- if(typet) {
        rep(0,length(nbcap))
      } else { 
        if (tinf) -log(factorial(nbcap)) else log(choose(t, t0:1)) 
      }
  return(cst)
}


tryCatch.W.E <- function(expr) { 
#?# Fonction pour stocker les erreurs et les warnings
#?# tiree de demo(error.catching), legerement modifiee
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W, w$message) ## ma modif ici pour stocker tous les warnings plutot que seulement le dernier
    invokeRestart("muffleWarning")
  }
  e.handler <- function(e){ # error handler
    class(e) <- "erreur"
    return(e)
  }
  list(value = withCallingHandlers(tryCatch(expr, error = e.handler), warning = w.handler),
       warnings = W)
}


glm.call <- function(Y, mX., cst=NULL, ...){
#?# Fonction qui appelle la fonction glmRcapture et conserve une trace
#?# des erreurs et avertissements obtenus, s'il y en a
#?# + ajoute un warning en cas de matrice de design pas de plein rang
  
  # (Je ne donne pas une matrice a droite de la formule car glm va alors nommer
  # les coefficients "nom_matrice.nom_colonne" et je veux avoir le nom d'une matrice
  # en prefixe uniquement si un argument mX a ete donne.)
  
  # (Aussi, je fournis manuellement l'intercept afin que glm() croit que le modele ne
  #  contient pas d'intercept et qu'il n'appelle pas une deuxieme fois glm.fit() afin
  #  de calculer correctement la deviance du modele NULL. Ainsi, il n'a pas d'ambuguite
  #  sur la provenance des avertissements ou erreurs generes par glm())
  mX.glm <- cbind(Intercept = rep(1, nrow(mX.)), mX.)
  data.glm <- as.data.frame(cbind(Y , mX.glm))
  formula.glm <- formula(paste("Y ~ 0 +", paste(colnames(mX.glm), collapse = " + ")))
  if (all(cst==0)) cst <- NULL
  
  ### On soumet la commande mais en enregistrant les warnings
  glm.out <- tryCatch.W.E(glm(formula = formula.glm, data = data.glm, family = poisson, offset = cst, ...))

  # Trace des warnings s'il y a lieu
  warn <- NULL
  if (!is.null(glm.out$warnings))
    warn <- c(warn, paste("glm:", glm.out$warnings))
  
  ### Trace de l'erreur s'il y a lieu
  if (inherits(glm.out$value, "erreur")) {
    # Si la commande a genere une erreur
    glmo <- NA
    error <- paste("glm:", glm.out$value$message)
  } else {
    # Si la commande a roule sans erreur
    glmo <- glm.out$value
    error <- NULL
    if (any(is.na(coef(glmo)))) {
      warn <- c(warn, "glm: design matrix not of full rank") 
      if (nrow(mX.) < ncol(mX.) + 1) {
        warn[length(warn)] <- paste0(warn[length(warn)], " (fewer rows than columns)")
      }
    }
  }
  
  # Sortie des resultats
  return(list(glmo = glmo, error = error, warnings = warn))   
}


Chao.neg <- function(glmo, mX., nca, ...)
{
#?# Fonction qui reajuste un modele de Chao s'il contient des parametres eta negatifs
  
  # On enleve les eta negatifs, mais certains eta peuvent prendre la valeur NA
  indicneg <- function(eta) { !is.na(eta) & eta < 0 }
  # cette fonction retourne un vecteur de logique de la meme longueur que le 
  # vecteur des eta donne en entree : TRUE si le parametre eta est negatif, FALSE autrement
  
  # Initialisation
  testneg <- any(indicneg(eta = coef(glmo)[-(1:nca)]))
  warn <- NULL
  nbrWarnPrev <- 0
  
  # Boucle
  while(testneg) {# Repeter la boucle jusqu'a ce qu'aucun eta ne soit negatif
  
    # Determination de la colonne a enlever dans mX.
    pos <- nca
    while(!indicneg(eta = coef(glmo)[pos+1])) pos <- pos + 1
    
    # Retrait de la bonne colonne de mX. et reajustement du modele
    mX. <- mX.[, -pos, drop = FALSE]  # pas pos + 1 car mX. n'a pas de colonne pour l'intercept
    glm.out <- glm.call(Y=glmo$y, mX.=mX., cst=glmo$offset, ...)
    
    # Si aucune erreur n'a ete generee, on poursuit normalement les calculs
    if (is.null(glm.out$error)) {
      etaRemoved <- names(coef(glmo))[pos+1]
      glmo <- glm.out$glmo
      testneg <- if(length(coef(glmo))>nca) any(indicneg(eta = coef(glmo)[-(1:nca)])) else FALSE

      # Si un warning a ete genere, on en garde une trace
      if (!is.null(glm.out$warnings)) {
        intro1 <- if (testneg) "warning not from the last fit" else "warning from the last fit"
        intro2 <- paste0("(after removing ", etaRemoved, "):")
        warn <- c(warn, paste(intro1, intro2, glm.out$warnings))
      }
      
    # Si une erreur a ete generee, on ne modifie pas glmo et on arrete la boucle
    } else {
      warn <- c(warn, paste("the negative eta parameter", names(coef(glmo))[pos+1], 
          "was not removed because of this error from glm when fitting the model without it:",
          glm.out$error))
      testneg <- FALSE
    }

  }
  
  return(list(glmo = glmo, warnings = warn))
# Note : possiblement a modifier pour les fonctions openp, robustd.t et robustd.0.
# pour enlever les gamma en plus des eta negatifs 
}

# Si je revise mes fonctions openp, robustd.t et robustd.0, j'envisage d'ajouter des
# fonctions, comme pour le calcul des parametres.


optim.call <- function(par, fn, method, ...) {
#?# Fonction qui appelle la fonction optim et conserve une trace
#?# des erreurs et avertissements obtenus, s'il y en a, 
#?# + ajoute un warning en cas de non-convergence
  
  optim.out <- tryCatch.W.E(optim(par = par, fn = fn, method = method, hessian=TRUE, ...))
  
  # Trace des warnings s'il y a lieu
  warn <- NULL
  if (!is.null(optim.out$warnings))
    warn <- c(warn, paste("optim:", optim.out$warnings))
  # warning non-convergence seulement si aucune erreur n'est generee

  # Trace de l'erreur s'il y a lieu
  if (inherits(optim.out$value, "erreur")) {
    # Si la commande a genere une erreur
    optimo <- NA
    error <- paste("optim:", optim.out$value$message)
  } else {
    # Si la commande a roule sans erreur
    optimo <- optim.out$value
    error <- NULL   
    if(optimo$converge!=0)
      warn <- c(warn, "optim: algorithm did not converge")
  }
  
  # Sortie des resultats
  return(list(optimo = optimo, error = error, warnings = warn))
}


getmu <- function(param,mX,nbcap,cst){
########################################################################################################
#?# Fonction qui calcule les valeurs predites par le "modele normal"
#?# On doit lui fournir :
#?# param: la valeur de tous les parametres du modele (l'ordonnee a l'origine 
#?#        + un parametre par colonne de la matrice mX + sigma, en respectant cet ordre)
#?# mX: la matrice X de la partie log-lineaire du modele (celle donnee en entree a closedpCI, 
#?#     sans colonne de 1 pour l'ordonnee a l'origine)
#?# nbcap: le nombre de captures par historique de capture (respectant le meme ordre que mX)
#?# cst: l'offset a ajouter au modele (0 pour fct .t, mais utile pour fct .0)
#?# Elle retourne un vecteur de longueur 2^t-1: les valeurs predites par le "modele normal" (mu)
########################################################################################################
  
  mXi <- cbind(rep(1,nrow(mX)),mX)
  mXi_nbcap <- cbind(mXi,nbcap)
  t <- max(nbcap)
  intercept <- param[1]
  coeffll <- param[-(ncol(mX)+2)]
  sig <- param[ncol(mX)+2]
  
  # Les points et les poids pour l'integration numerique par la methode de quadrature de Gauss
  # tire de Aramowitz et Stegun (1972), page 924 (Hermite integration)
  # On peut aussi retrouver ces valeurs en ligne :
  # http://www.efunda.com/math/num_integration/findgausshermite.cfm
  xnor <- c(-5.38748089001, -4.60368244955, -3.94476404012, -3.34785456738, -2.78880605843,
            -2.25497400209, -1.73853771212, -1.2340762154, -0.737473728545, -0.245340708301,
            0.245340708301, 0.737473728545, 1.2340762154, 1.73853771212, 2.25497400209,
            2.78880605843, 3.34785456738, 3.94476404012, 4.60368244955, 5.38748089001)
  wnor <- c(2.22939364554E-013, 4.39934099226E-010, 1.08606937077E-007, 7.8025564785E-006,
            0.000228338636017, 0.00324377334224, 0.0248105208875, 0.10901720602, 0.286675505363,
            0.462243669601, 0.462243669601, 0.286675505363, 0.10901720602, 0.0248105208875,
            0.00324377334224, 0.000228338636017, 7.8025564785E-006, 1.08606937077E-007,
            4.39934099226E-010, 2.22939364554E-013)
  
  # Calcul des valeurs de la constante C pour les 20 points de l'integration numerique 
  # (C(sigma*z_i) dans la formule de Rivest (2011), page 7)
  # C(alpha)=p(0|alpha) conditionnal probability of not appearing in any list
  C<-vector(length=20)
  for ( i in (1:20)){ C[i]<-exp(intercept)/(exp(intercept)+sum(exp(cst+mXi_nbcap%*%c(coeffll,xnor[i]*sqrt(2)*sig)))) }
  
  # Calcul de la fonction varphi (formule de Rivest (2011), page 7)
  phi<-vector(length=t)
  for (i in (1:t)){ phi[i]<-log(sum(wnor*exp(i*sqrt(2)*sig*xnor)*C)/sum(wnor*C)) }
  
  # Calcul des valeurs predites
  mu<-exp(cst + mXi%*%coeffll + phi[nbcap]) # eq3, Rivest (2011), page 3
  mu <- as.vector(mu)
  names(mu) <- 1:length(mu)
  return(mu)
}


devPois <- function(param,Y,mX,nbcap,cst) {
#######################################################################################################
#?# Fonction qui sera donnee en entree a optim.
#?# Elle calcule la valeur de la statistique a minimiser (la deviance Poisson) pour le "modele normal".
#?# On doit lui fournir :
#?# param: la valeur de tous les parametres du modele (l'ordonnee a l'origine 
#?#        + un parametre par colonne de la matrice mX + sigma, en respectant cet ordre)
#?# Y: les frequences observees (respectant le meme ordre que mX)
#?# mX: la matrice X de la partie log-lineaire du modele (celle donnee en entree a closedpCI, 
#?#     sans colonne de 1 pour l'ordonnee a l'origine)
#?# nbcap: le nombre de captures par historique de capture (respectant le meme ordre que mX)
#?# cst: l'offset a ajouter au modele (0 pour fct .t, mais utile pour fct .0)
#?# Elle retourne un seul chiffre: la valeur de la deviance poisson
#######################################################################################################
  
  mu <- getmu(param,mX,nbcap,cst)
  dev.resids <- function(y, mu) { 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu)) }   
  # derive de family.R ligne 169
  dev <- sum(dev.resids(y=Y, mu=mu)) # derive de glm.R ligne 183
  return(dev)
}


closedp.fitone <- function(n, Y, mX., nbcap, nca, cst, htype, neg, 
                           initcoef = NULL, initsig = NULL, method = NULL, ...) {    
#?# Fonction interne importante, appelee par les fonctions closedp, closedpCI, closedp.bc
#?# Elle ajuste un modele. a partir des estimations obtenues pour les parametres
#?# de ce modele, on estime N.
#?# L'ajustement se fait presque toujours avec la fonction glm (regression poisson),
#?# mais elle se fait avec optim pour des modeles heterogenes normaux.
#?# Elle conserve une trace des erreurs et les warnings obtenus lors de l'ajustement.
  
  fit.err <- fit.warn <- NULL
  
  if (htype == "Normal") {  ## Si un modele heterogene normal a ete demande :
        
    # Avant meme d'ajuster le modele, on verifie si le nombre de lignes 
    # d'observations pour l'ajustement du modele est bien inferieur ou egal 
    # au nombre de parametres a estimer (pour un modele ajuste avec glm,
    # cette situation genere seulement un avertissement car glm laisse de 
    # lui-meme tomber certains parametres a estimer).
    nrows <- nrow(mX.)
    nparams <- ncol(mX.) + 2 
    if (nrows < nparams) {
      fit.err <- paste0("the model was not fitted because the number of parameters (", 
                        nparams, ") is larger than the number of observations to fit the model (", 
                        nrows, "), therefore some model's coefficients are not estimable")
    }
      
    #########  Determination des valeurs initiales des parametres pour optim
    if(is.null(initcoef)) {       
      # Si non fournies : parametres obtenus du modele mX. + colonne Darroch
      mX.D <- cbind(mX.,tau=(nbcap^2)/2)
      glmD <- glm.call(Y=Y, mX.=mX.D, cst=cst)  
        # pas de ... ici car on ne peut pas a la fois passer des arguments a glm et a optim
      if (!is.null(glmD$error)) {
        fit.err <- c(fit.err, paste("error from a call to glm to determine 'initcoef':", 
                                    glmD$error))
      } else {
        initcoeffll <- coef(glmD$glmo)[1:(ncol(mX.)+1)]
        if (!is.null(glmD$warnings))
          fit.warn <- c(fit.warn, paste("warning from a call to glm to determine 'initcoef':",
                                        glmD$warnings))
        initcoef.converge <- glmD$glmo$converged
      }
    } else {
      # Si valeurs initiales fournies en entree, aucun calcul a faire
      initcoeffll <- initcoef
      initcoef.converge <- TRUE
    }

    # On fait les etapes suivantes uniquement si on a en main des valeurs initiales
    if (is.null(fit.err)) {

      names(initcoeffll) <- c("Intercept", colnames(mX.))
      # Ajout de la valeur initiale de sigma
      initparam <- c(initcoeffll, sigma=initsig)
      
      ######### Ajustement du modele normal avec optim
      
      ### On soumet la commande mais en enregistrant les warnings
      optim.out <- optim.call(par = initparam, fn = devPois, method=method,
                              Y=Y, mX=mX., nbcap=nbcap, cst=cst, ...)
      fit.err <- c(fit.err, optim.out$error)
      
      # Si la commande a roule sans erreur
      if (is.null(optim.out$error)) {
        
        optimo <- optim.out$optimo
        fit.warn <- c(fit.warn, optim.out$warnings)
        
        # Un warning de plus est possible
        if(optimo$par[length(optimo$par)] <= 0) 
          fit.warn <- c(fit.warn, "the sigma estimate is non-positive, it seems there is no heterogeneous catchability in the data")  
      
        ######### Calculs pour le tableau des resultats
        dev <- optimo$value
        dff <- nrow(mX.) - length(optimo$par)
        mu <- getmu(param=optimo$par,mX.,nbcap,cst)
        dobs <- dpois(Y, mu, log=TRUE)
        aic <- -2*sum(dobs) + 2*length(optimo$par)
        # le premier element de la somme est derive de family.R ligne 170
        # le deuxieme element de la somme est 2 fois le nombre de parametres dans le modele
        bic <- -2*sum(dobs) + log(n)*length(optimo$par)
        # On ne calcul pas de biais asymptotique pour le modele normal
        
        ######### Preparation de la sortie
        resultsFit <- c(NA, dev, dff, aic, bic)
        
        varcov <- solve(optimo$hessian)
        parameters <- rbind(estimate=optimo$par, stderr=diag(varcov))
        fit <- list(parameters=parameters,varcov=varcov,y=Y,fitted.values=mu,
                    initparam=initparam,optim=optimo)
      
      }
    } 
    
  } else { ## Si on n'a pas demande un modele heterogene normal : 

    ##### Ajustement du modele avec un modele loglineaire poisson 
    glm.out <- glm.call(Y=Y, mX.=mX., cst=cst, ...)
    fit.err <- c(fit.err, glm.out$error)
        
    # Si la commande a roule sans erreur
    if (is.null(glm.out$error)) {
      
      glmo <- glm.out$glmo
      fit.warn <- c(fit.warn, glm.out$warnings)
    
      ##### Reajustement du modele en enlevant les eta negatifs pour les modeles de Chao
      if (htype == "Chao" && neg) {
        neg.out <- Chao.neg(glmo=glmo, mX. = mX., nca = nca)
        neg.eta <- setdiff(names(coef(glmo)), names(coef(neg.out$glmo)))
        glmo <- neg.out$glmo
        # Trace des avertissements
        if(!is.null(neg.out$warnings)) {
          if(is.null(fit.warn) || (length(neg.out$warnings) == 1 && grepl("the negative eta", neg.out$warnings))) {
            intro <- NULL
          } else { 
            intro <- "warning not from the last fit (with every eta parameters):"
          }
          fit.warn <- c(paste(intro, fit.warn), neg.out$warnings)
        }
      } else { neg.eta <- NA }
      
      ##### Calcul du biais asymptotique (code de Louis-Paul (generalise))
      if (any(is.na(coef(glmo)))) {
        # Si la matrice de design est singuliere, on ne calcule pas le bias
        bias <- NA
      } else {
        mX.i <- model.matrix(glmo) 
        # cette matrice contient bien l'ordonnee a l'origine
        # on n'utilise pas directement mX. car elle peut avoir ete modifie
        # en enlevant des colonnes pour les eta negatifs
        inv.fisher.info <- vcov(glmo)
        xi <- -0.5 * diag(mX.i %*% inv.fisher.info %*% t(mX.i))
        psi <- -2 * t(mX.i) %*% diag(fitted(glmo)) %*% xi
        psi[1] <- psi[1] - 1
        vec1 <- -inv.fisher.info %*% psi
        bias <- 0.5*exp(coef(glmo)[1])*vec1[1]
      }
      
      ##### Preparation de la sortie
      bic <- glmo$aic + (-2 + log(n))*glmo$rank
      resultsFit <- c(bias, glmo$dev,glmo$df.residual,glmo$aic, bic)
      fit <- glmo
      
    }
    
  }
  
  # Si on a obtenu une erreur qui a empeche l'ajustement du modele
  if (!is.null(fit.err)) {
    resultsFit <- c(NA, NA, NA, NA, NA)
    fit <- NA
    if (htype == "Chao") neg.eta <- NA
  }
  
  ##### Sortie
  names(resultsFit) <- c("bias","deviance","df","AIC","BIC")
  ans <- list(resultsFit = resultsFit, fit = fit, fit.err = fit.err, fit.warn = fit.warn)
  if (htype == "Chao") ans <- c(ans, list(neg.eta = neg.eta)) 
  ans 
}


getN <- function(n, intercept, stderr.intercept) {
#?# Calcul standard de N et de son erreur type (pas pour modeles Mb et Mbh)
  c("abundance" = as.vector(n + exp(intercept)), 
    "stderr" = as.vector(sqrt(exp(intercept) + exp(2*intercept)*stderr.intercept)))
}


pres <- function(x, typet, t) {
#?# Pour calculer les residus de Pearson
#?# x doit etre un objet de classe 'glm'
  if (typet) {
    nbcap <- rowSums(histpos.t(t))
    fi <- tapply(x$y,nbcap,sum,na.rm=TRUE)
    fipred <- tapply(x$fitted.values,nbcap,sum,na.rm=TRUE)
  } else {
    fi <- rev(x$y)
    fipred <- rev(x$fitted.values)
  }
  (fi-fipred)/sqrt(fipred)
}


getBiasWarn <- function(N, bias) {
  #?# Afin de tester si le bias est grand et d'ecrire l'avertissement associe  
  if(!is.na(bias) && abs(bias) > 0.1*N) "the asymptotic bias is large" else NULL
}


getInfo <- function(err, warn) {
#?# Afin de determiner la valeur de infoFit ou infoCI qui
#?# decrit le type des conditions (erreur ou avertissements)  
#?# generes par l'ajustement du modele ou la calcul du CI  
  
  if (!is.null(err)) {
    infoVect <- -1
  } else if (is.null(warn)) { 
    infoVect <- 0
  } else {
    iConverge <- grepl("converge", warn, fixed = TRUE)
    iSigma <- grepl("the sigma estimate is non-positive", warn, fixed = TRUE)
    iBias <- grepl("the asymptotic bias is large", warn, fixed = TRUE)
    iSingular <- grepl("design matrix not of full rank", warn, fixed = TRUE)
    iMultiN <- grepl("multinomial estimation", warn, fixed = TRUE)
    iInfCL <- grepl("lower bound", warn, fixed = TRUE)
    iSupCL <- grepl("upper bound", warn, fixed = TRUE)    
    infoVect <- 
      ifelse(iConverge | iSigma | iBias, 1,
        ifelse(iSingular, 2,
          ifelse(iMultiN, 4,
            ifelse(iInfCL, 5,
              ifelse(iSupCL, 6, 3)))))
  }
  
  # concatenation des chiffres distincts dans infoVect en un scalaire
  as.numeric(paste(sort(unique(infoVect)), collapse = ""))
}


tabprint <- function(tab, digits, warn, ...) {
#?# Fonction pour afficher comme je le veux un tableau
#?# tab <- matrice contenant le tableau
#?# digits <- vecteur de longueur egale au nombre de colonnes de tab,
#?#           contenant le nombre de digits pour l'arrondissement de chaque colonne,
#?#           la valeur speciale NA signifie que la colonne contient le code
#?#           d'information sur l'ajustement du modele.
  
  tab <- as.data.frame(tab)
  
  # Arrondissement des colonnes a arrondir
  for (i in (1:ncol(tab))[!is.na(digits)]) {
    if (is.numeric(tab[,i])) tab[,i] <- round(tab[,i], digits[i])
    # Le test numerique est la, car dans la sortie du profile multinomial 
    # likelihood confidence interval il est possible qu'une colonne soit
    # sous forme de caractere (ex.">200").
  } 
  
  # S'il y a une colonne avec un code
  if(any(is.na(digits))) {
    signe <- sign(tab[, is.na(digits)])
    code <- as.character(tab[, is.na(digits)])
    nbrWarn <- sapply(warn, length)
    tab[, is.na(digits)] <- 
      ifelse(is.na(signe), "no fit",
        ifelse(signe == 0, "OK",
          ifelse(signe == -1, "error",
            ifelse(nchar(code) == 1, paste0(ifelse(nbrWarn>1, "warnings", "warning"), " #", code),
              paste0("warnings", gsub("([1-9])", " #\\1", code))))))
  }
  
  # Affichage du tableau
  print(tab, ...)
}


getFormulaFromName <- function(name) {
#?# ecrit les formules a fournir en entree a closedpCI.t en argument mX a partir 
#?# de noms de modele
#?# Exemple : le nom [123,4] va donner la formule :
#?#           ~ 1 + c1 + c2 + c3 + c4 + c1:c2 + c1:c3 + c2:c3 + c1:c2:c3 
#?# @param name vecteur contenant les noms de modeles (doit respecter le format
#?#             decrit dans la documentation de la fonction getAllModels)
#?#             (***pas de validation faite car cette fonction demeurera interne,
#?#              IMPORTANT : ne doit pas contenir d'espaces) 
#?# @return liste contenant les formules pour chaque nom
  
  # Afin de produire une liste de vecteurs contenant les elements dans les noms
  high <- strsplit(substr(name, 2, nchar(name) - 1), split = ",")
  
  # Boucle sur tous les noms
  formules <- vector(mode = "list", length = length(high))
  for (i in 1:length(high)) {
    
    # terms = Liste des termes dans la formule
    # on y inclut d'abord les termes formant le nom du modele
    terms <- high[[i]]
    
    # Si l'ordre le plus eleve present dans les termes est superieur a 1,
    # boucle sur les ordres des termes, en partant de l'ordre le plus grand
    if(length(terms) != 0 && max(nchar(terms)) > 1) {
      for(j in max(nchar(terms)):2) {
        
        # Boucles sur tous les termes d'ordre j
        for(k in which(nchar(terms) == j)) {
          
          # enumeration des termes d'un ordre inferieur a inclure
          # dans la formule par hierarchie
          combin <- combn(unlist(strsplit(terms[k], split = "")), j-1)
          
          # Ajout de ces termes dans le vecteur des termes du modele
          terms <- c(terms, apply(combin, 2, paste, collapse=""))
        }
        
        # Il y a souvent des termes dupliques a la fin de cette
        # boucle, on doit retirer les doublons
        terms <- unique(terms)
      }
      
      # Afin d'ajouter les "c:" aux termes pour les ecrire de la facon
      # comprise par closedpCI.t
      terms <- gsub(pattern = "([1-9])", replacement = ":c\\1" , x = terms)
      terms <- substr(terms, 2, nchar(terms)) 
      
    } else {
      
      # Si le modele ne contient que des effets simples, 
      # pas besoin d'utiliser gsub pour ajouter les "c:".
      if (length(terms) != 0) {
        terms <- paste0("c", terms)
      }
      # Il n'y a rien a faire si le modele ne comprend qu'une ordonnee a l'origine
    }
    
    # Ordonner les termes d'abord en ordre croissant d'ordre des termes,
    # puis, pour un meme ordre, en ordre alphabetique des noms des termes
    permut <- order(nchar(terms), terms)
    # Formation de la formule, en incluant toujours une ordonnee a l'origine
    formules[[i]] <- formula(paste("~", paste(c("1", terms[permut]), collapse = " + ")))
  }
  
  # Sortie : on retourne la liste des formules
  names(formules) <- name
  formules
  
}

