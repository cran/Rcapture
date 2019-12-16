closedpMS.t <- function(X, dfreq = FALSE, h=NULL, h.control=list(), 
                        maxorder = t - 1, forced = 1:t, stopiflong = TRUE, ...) {
  
  call <- match.call()
  
  #########  Validation des arguments en entree et initialisation de variables  #########
  
  # Validation des arguments et obtention de la valeur de t
  valid.one(dfreq,"logical")
  Xvalid <- valid.X(X=X, dfreq=dfreq)
  X <- Xvalid$X
  t <- Xvalid$t
  
  valid.h.out <- valid.h(h=h, m=NULL, call=call,
                                    values=c("Chao","LB","Poisson","Darroch","Gamma","Normal"))
  h <- valid.h.out$h
  htype <- valid.h.out$htype
  
  if(!is.list(h.control)) stop("'h.control' must be a list")
  theta <- valid.theta(theta = h.control$theta, htype = htype)
  neg <- valid.neg(neg = h.control$neg, htype = htype)
  initsig <- valid.initsig(initsig = h.control$initsig, htype = htype)
  method <- valid.method(method = h.control$method, htype = htype)
  # initcoef ne peut pas etre fourni car le nombre de parametres change selon le modele

  # maxorder, forced et stopiflong sont valides lors de l'appel a la fonction getAllModels
  
#   if(is.null(call$maxorder)) {
#     maxorder <- t - 1
#   } else {
#     if (!(maxorder %in% 1:(t-1))) 
#       stop("'maxorder' must be an integer between 1 and 't' - 1 inclusively")
#   }
#   
#   valid.one(stopiflong, type = "logical")
#   if(stopiflong && (t==6 && maxorder >=3 || t>=7 && maxorder >=2))
#     stop("the number of models to fit is large, therefore this command should be long to run,",
#          "\nset 'stopiflong' to FALSE if you really want to run it")
  
  #########  Creation de variables communes a tous les modeles  #########   
  
  ### Creation du vecteur de variable reponse Y
  Y <- histfreq.t(X=X,dfreq=dfreq)
  n <- sum(na.rm=TRUE,Y)
  
  ### Creation de la variable offset
  cst <- 0  
  
  
  #########  Boucle sur tous les modeles  ######### 
  
  # Identification de tous les modeles a ajuster
  name <- getAllModels(t=t, maxorder=maxorder, forced=forced, stopiflong=stopiflong)
  modeles <- getFormulaFromName(name=name)
  
  # Objets pour stocker les resultats
  tableau <- matrix(NA_real_, nrow=length(modeles), ncol=7)
  dimnames(tableau) <- list(name,c("abundance","stderr","deviance","df","AIC","BIC","infoFit"))
  bias <- rep(NA_real_, length=length(modeles))
  fit.err <- vector(mode="list", length=length(modeles))
  fit.warn <- vector(mode="list", length=length(modeles))
  names(fit.err) <- names(fit.warn) <- name
  if (htype == "Chao") {
    neg.eta <- vector(mode="list", length=length(modeles))
    names(neg.eta) <- name
  }
  
  # Boucle qui ajuste tous les modeles
  for (j in 1:length(modeles)) {
    
    ### Creation de la matrice X (qui est propre au modele)
    getmX.out <- getmX(typet=TRUE, t=t, t0=t, m=NULL, h=h, theta=theta, 
                                  mX=modeles[[j]])
    mX. <- getmX.out$mX. 
    nbcap <- getmX.out$nbcap
    nca <- getmX.out$nca
    
    ### Ajustement du modele
    fit.out <- closedp.fitone(n = n, Y = Y, mX. = mX., nbcap = nbcap, nca = nca, 
                                         cst = cst, htype = htype, neg = neg, 
                                         initsig = initsig, method = method, ...)
    
    if(!is.null(fit.out$fit.err)) fit.err[[j]] <- fit.out$fit.err
    if(!is.null(fit.out$fit.warn)) fit.warn[[j]] <- fit.out$fit.warn
    if (htype == "Chao") neg.eta[[j]] <- fit.out$neg.eta
    
    # Calculs a faire seulement si glm a produit une sortie
    # (on laisse dans tableau et param les NA mis lors de l'initialisation
    #  en cas d'erreur lors de l'ajustement du modele)
    if (is.null(fit.out$fit.err)) {
      
      # On met les statistiques d'ajustement du modele dans tableau et bias
      tableau[j, 3:6] <- fit.out$resultsFit[-1]
      bias[j] <- fit.out$resultsFit[[1]]
      
      # On estime N
      intercept <- if(htype == "Normal") fit.out$fit$parameters[1, 1] else coef(fit.out$fit)[1]
      stderr.intercept <- if(htype == "Normal") fit.out$fit$varcov[1, 1] else vcov(fit.out$fit)[1, 1]
      tableau[j, 1:2] <- getN(n = n, intercept = intercept, stderr.intercept = stderr.intercept)
      
      # Avertissement pour grand biais au besoin
      biasWarn <- getBiasWarn(N = tableau[j, "abundance"], bias = bias[j])
      if(!is.null(biasWarn)) fit.warn[[j]] <- c(fit.warn[[j]], biasWarn) 
      
    }
    
    # code pour les conditions mis dans le tableau
    tableau[j, 7] <- getInfo(err = fit.err[[j]], warn = fit.warn[[j]])
    
  }
  
  # Sortie des resultats
  ans <- list(n = n, t = t, results = tableau, bias = bias, fit.err = fit.err, fit.warn = fit.warn)
  if (htype == "Chao") ans <- c(ans, neg.eta)
  class(ans) <- "closedpMS"
  ans
}


print.closedpMS <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")
  
  if (!is.null(x$results)) {
    cat("Abundance estimations and model fits for the models with the smallest BIC:\n")
    tableau <- x$results[order(x$results[, "BIC"]), ]
    tabprint(tab = tableau[1:min(10,nrow(tableau)), ], 
                        digits = c(1,1,3,0,3,3,NA), warn = x$fit.warn, ...)
  }
  
  cat("\n")
  invisible(x)
}


plot.closedpMS <- function(x, main="Models comparison based on BIC", omitOutliers = TRUE, ...){
  N <- x$results[, "abundance"]
  BIC <- x$results[, "BIC"]
  
  if (omitOutliers) {
    # On retire les valeurs extremes : < q1 - 1.5*IQR  ou  > q3 + 1.5*IQR
    Nq1 <- quantile(N,0.25)
    Nq3 <- quantile(N,0.75)
    keepN <- N > Nq1 - 1.5*(Nq3-Nq1) & N < Nq3 + 1.5*(Nq3-Nq1)
    
    BICq1 <- quantile(BIC,0.25)
    BICq3 <- quantile(BIC,0.75)
    keepBIC <- BIC > BICq1 - 1.5*(BICq3-BICq1) & BIC < BICq3 + 1.5*(BICq3-BICq1)
    
    N <- N[keepN & keepBIC]
    BIC <- BIC[keepN & keepBIC]
  }
  
  # Graphique
  plot(x = N, y = -BIC, xlab = "abundance", ylab = "-BIC", main = main, ...)
  abline(v = N[order(BIC)][1], col= "blue") 
}


# ---------------------------------------------------------------------------- #


getAllModels <- function(t, maxorder = t - 1, forced = 1:t, stopiflong = TRUE) {
  
  # Validations et initialisations
  if (!(t %in% 2:9)) stop("'t' must be an integer between 2 and 9 inclusively")
  
  if (!(maxorder %in% 1:(t-1))) 
    stop("'maxorder' must be an integer between 1 and 't' - 1 inclusively")
  
  forced <- unique(as.character(forced))
  if (!all(unlist(strsplit(forced, "")) %in% 1:t))
    stop("the only accepted characters in 'forced' are ", 
         paste(1:(t-1),collapse=", "), " and ", t)
  
  valid.one(stopiflong, type = "logical")
  if(stopiflong && (t==6 && maxorder >=3 || t>=7 && maxorder >=2))
    stop("the number of models to enumerate is large, ",
         "therefore this command should be long to run,",
         "\nset 'stopiflong' to FALSE if you really want to run it")
  
  # Informations sur les termes forces dans le modele
  forcedOrder <- nchar(forced)
  minorder <- if(length(forced) == 0) 1 else max(forcedOrder)
  if (minorder > maxorder) 
    stop("the order of one of the forced term is larger than 'maxorder'")
  
  # Liste qui contiendra tous les modeles possibles
  modeles <- list()
  # (un modele = un vecteur des termes de haut de hierarchie formant son nom)
  # (l'objet modeles change de taille au fil des iterations, car il n'y a pas de
  #  formule simple pour calculer le nombre total de modeles)
  
  #---------------------------------------------------------------------------#
  
  # Boucle sur les ordres
  # on insere dans les modeles les termes un ordre a la fois,
  # en debutant par les termes de l'ordre le plus eleve possible, soit t-1
  for (i in maxorder:1) {
    
    # Tous les termes possibles de cet ordre
    termesPossibles <- vapply(combn(1:t,i, simplify = FALSE), 
                              paste, collapse="", FUN.VALUE = "a")
    
    # Termes de cet ordre forces dans le modele
    iforced <- forced[forcedOrder == i]
    
    # etape A : integration de termes aux modeles precedents, un modele a la fois
    if (length(modeles) > 0) {
      for (j in 1:length(modeles)) {
        
        # Identification des termes d'ordre inferieur qui sont obligatoirement
        # presents dans le modele (par hierarchie), qui ne doivent donc pas
        # apparaitre dans le nom du modele
        termesInfForces <- NULL
        for (k in 1:length(modeles[[j]])) {
          combin <- combn(unlist(strsplit(modeles[[j]][k], split = "")), i, simplify = FALSE)
          termesInfForces <- c(vapply(combin, paste, collapse="", FUN.VALUE = "a"),
                               termesInfForces)
        }
        iforcedNonHierar <- iforced[!(iforced %in% termesInfForces)]
        
        # Identification des termes qui pourraient apparaitre dans le nom du modele
        termesInfLibres <- setdiff(termesPossibles, termesInfForces)
        termesInfLibresNotForced <- setdiff(termesInfLibres, iforced)
        
        # Boucle sur le nombre de termes d'ordre inferieur non forces 
        # integres au modele
        if (length(termesInfLibresNotForced) > 0) {          
          for (k in length(termesInfLibresNotForced):1) {
            
            # enumeration de toutes les combinaisons de termes d'ordre inferieur
            # libres et non forces
            termesInfChoisis <- combn(termesInfLibresNotForced, k, simplify = FALSE)
            
            # Ajout des termes forces s'il y en a et s'ils ne sont pas
            # presents par hierarchie
            termesInfChoisis <- lapply(termesInfChoisis, c, iforcedNonHierar)
            
            # Ajout de ces combinaisons aux termes deja dans le modele
            ajout <- lapply(termesInfChoisis, c, modeles[[j]])
            # Pour remettre les termes dans l'ordre voulu
            na <- k + length(iforcedNonHierar)
            ajout <- lapply(ajout, '[', c((na+1):(na+length(modeles[[j]])), 1:na))
            
            # On ajoute ces modeles a la liste de tous les modeles
            modeles <- c(modeles, ajout)
          }
        }
        
        # Ajout des termes forces de l'ordre i sans autres termes de l'ordre i
        if (length(iforcedNonHierar) > 0) {
          modeles <- c(modeles, list(c(modeles[[j]], iforcedNonHierar)))
        }
        
        # Il faut retirer le modele j de la liste si :
        # - au moins un terme force est dans termesInfLibres 
        #   (car le modele ne contient pas le terme en question)
        # - 
        if(length(termesInfLibres) > length(termesInfLibresNotForced)) { 
          modeles[[j]] <- character(0)
        }
        
      }  
    }
    
    # etape B : ajout de modeles sans termes d'ordre supperieur a i
    # On n'insere que les modeles contenant les termes forces de
    # l'ordre i traite a cette iteration, s'il y en a.   
    if(i >= minorder) {
      
      termesPossiblesNotForced <- setdiff(termesPossibles, iforced)
      
      # Boucle sur le nombre de termes possibles non forces
      if (length(termesPossiblesNotForced) > 0) {
        for (j in length(termesPossiblesNotForced):1) {
          modeles <- c(modeles, 
                       lapply(combn(termesPossiblesNotForced,
                                    j, simplify = FALSE), 
                              c, iforced))
        }
      }
      
      # Ajout du modele contenant uniquement les termes forces d'ordre i,
      # s'il y a des termes forces de cet ordre
      if (length(iforced) > 0) modeles <- c(modeles, list(iforced))
    }
    
    # etape C : Retirer les elements de la liste qui sont des 
    # vecteurs vides, s'il y a lieu
    modeles <- modeles[lapply(modeles, length) > 0]
    
  } 
  
  # Pour tous les modeles, on met les termes formant sont nom entre crochets,
  # separes par des virgules
  noms <- paste0("[", vapply(modeles, paste, collapse = ",", FUN.VALUE = "a"), "]")
  # On ajoute le modele contenant seulement une ordonnee a l'origine,
  # s'il n'y a pas de termes forces dans le modele
  if(length(forced) == 0) noms <- c(noms, "[]")
  # Sortie : on retourne les noms trouves
  noms
}


