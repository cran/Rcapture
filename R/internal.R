###################################################################################################
# Fonctions pour des bouts de code apparaissant plus d'une fois dans mes fonctions

getY <- function(typet, X, dfreq, dtype, t, t0) {
# Création du vecteur de variable réponse Y
  Y <- if (typet) {
        histfreq.t(X=X,dfreq=dfreq)
      } else {
        histfreq.0(X=X,dfreq=dfreq,dtype=dtype[1],vt=t)
      }
  # Valeur dérivée de ce calcul : n
  n <- sum(na.rm=TRUE,Y)
  # Modification de Y pour fonction .0 (typet=FALSE) au besoin (t0 différent de t).
  if (!typet) Y <- modifYt0(Y=Y, t=t, t0=t0)
  return(list(Y=Y, n=n))
}

modifYt0 <- function(Y, t, t0){
# Modification de Y si t0 est différent de t
  if (t0 < t) {
    Y <- Y[(1 + t - t0):t] ## On enlève les première fréquences du vecteur, soit celle pour nbcap>t0.    
  } else if (t0 > t) { ## cas ou l'argument t=Inf et t0 est supérieur à tmax (maintenant = t)
    Y <- c(rep(0, t0 - t), Y) ## on doit ajouter des zéros lorsque plus de t captures
  }
  return(Y)
}

gethistpos <- function(typet, t, t0=t) {
  histpos <- if (typet) {
        histpos.t(t=t)
      } else {
        histpos.0(vt=t0)
      }
  return(histpos)
}

getnbcap <- function(histpos) { rowSums(histpos) }

getmX <- function(typet, t, t0, m, h=NULL, theta=NULL, mX=NULL) {
##### Description de la fonction : Création de la matrice X
  
  ### Création préliminaire d'objets nécessaires aux calculs   
  histpos <- gethistpos(typet=typet, t=t, t0=t0)
  nbcap <- getnbcap(histpos)
  t. <- if (typet) t else t0
  
  ### Si un mX sous forme de formule a été donné, on doit d'abord
  ### transformer la formule en la matrice équivalente.
  # (On sait à cause de la validation de mX que si mX est une formule
  #  alors nécessairement typet=TRUE, donc t0=NULL) 
  if (class(mX)=="formula"){   
    dfhistpos <- as.data.frame(histpos)
    colnames(dfhistpos) <- paste("c",1:t,sep="")    
    Terms <- terms(mX, data=dfhistpos)
    mX <- model.matrix(mX, data=dfhistpos)
    # Pour traiter l'ordonnée à l'origine
    if (attr(Terms, "intercept")==0){  ## si l'intercept a été enlevé : erreur
      stop("the intercept must not be removed in the formula given as 'mX' argument", call. = FALSE)
    } else {  ## si l'intercept est là c'est bien, mais on l'enlève il sera remis plus tard
      mX <- mX[,-1]  
    }
    # Pour changer les noms des colonnes dans cette matrice : 
    rsplit <- strsplit(colnames(mX), split="c")
    vrsplit <- sapply(rsplit,paste,collapse="")
    id <- nchar(vrsplit)
    # si un seul chiffre, on paste "beta" devant
    # si plus d'un chiffre, on paste "lambda" devant
    # (Je suis ici la notation de Rivest (2011), mais sans les paranthèses
    #  car R pense que lambda est une fonction et la virgule est remplacée par _.)
    pnames <- paste(ifelse(id==1,"beta","lambda"), vrsplit, sep="")
    # puis on remplace le ":" par "_".
    colnames(mX) <- gsub(":", "_", pnames)
  }
  
  ### Ensuite on peut obtenir la matrice de design mX.
  # (mX et mX. ne sont pas nécessairement équivalentes puisque des colonnes
  #  peuvent être ajoutées à mX pour l'hétérogénéité.)  
  Xclosedp.out <- Xclosedp(t=t., m=m, h=h, theta=theta, histpos=histpos, nbcap=nbcap, mX=mX)
  mX. <- Xclosedp.out$mat
  colnames(mX.) <- Xclosedp.out$paramnames
  
  ### Sortie des résultats
  nca <- if (is.null(m)) ncol(mX) + 1 else if (m=="Mh") 2 else t+1   
    ## nca : le nombre de colonnes ne modélisant pas l'hétérogénéité
    ## utile pour la correction des eta négatifs pour le modèle Chao ou LB 
  return(list(mX.=mX., nbcap=nbcap, nca = nca))
}

getcst <- function(typet, tinf, t, t0, nbcap){
# Création de la variable offset
  cst <- if(typet) {
        rep(0,length(nbcap))
      } else { 
        if (tinf) -log(factorial(nbcap)) else log(choose(t, t0:1)) 
      }
  return(cst)
}

pres <- function(x, typet, t) {
# Pour calculer les résidus de Pearson
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

tryCatch.W.E <- function(expr)
{ # Fonction pour stocker les erreurs et les warnings
  # tirée de demo(error.catching), légèrement modifiée
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W, w$message) ## ma modif ici pour stocker tous les warnings plutôt que seulement le dernier
    invokeRestart("muffleWarning")
  }
  e.handler <- function(e){ # error handler
    class(e) <- "erreur"
    return(e)
  }
  list(value = withCallingHandlers(tryCatch(expr, error = e.handler),
          warning = w.handler),
      warnings = W)
}

"glm.call"<-function(Y, mX., cst=NULL, mname){
  ### Écriture de la commande glm
  # (Je ne donne pas une matrice à droite de la formule car glm va alors nommer
  # les paramètres "nom_matrice.nom_colonne" et je veux avoir le nom d'une matrice
  # en préfixe uniquement si un argument mX a été donné.)
  data.glm <- as.data.frame(cbind(Y,mX.))
  glmcall <- if (all(cst==0)||is.null(cst)) {
    paste("glmRcapture(Y ~ ",paste(colnames(mX.),collapse=" + "),
          ", data=data.glm, family=poisson, control=list(maxit = 200))",sep="")
  } else {
    data.glm <- cbind(data.glm, cst)
    paste("glmRcapture(Y ~ offset(cst) + ",paste(colnames(mX.),collapse=" + "),
          ", data=data.glm, family=poisson, control=list(maxit = 200))",sep="")
  }

  ### On soumet la commande mais en enregistrant les warnings
  glm.out <- tryCatch.W.E(eval(parse(text=glmcall)[[1]]))
  if (identical(class(glm.out$value), "erreur"))
    stop("when fitting the model ", mname, " : ", gettext(glm.out$value$message), call. = FALSE)
    # glm.out$value$message sera peut-être dans la langue de l'utilisateur. C'est bien.
  
  ### Préparation de la sortie 
  return(list(glmo=glm.out$value, warnings=glm.out$warnings))   
}

"Chao.neg" <- function(glmo, nca, mname)
{
    testneg<-any(glmo$coef[-(1:nca)]<0)
    nocorrect <- !testneg
    while(testneg) {# Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
        # Détermination de la colonne à enlever dans mX
        pos<-nca
        while(glmo$coef[pos+1]>0) pos <- pos + 1
        # Retrait de la bonne colonne de mX et réajustement du modèle
        mX. <- if(is.null(glmo$offset)) glmo$model[,-1] else glmo$model[,-(1:2)]
        mX. <- mX.[,-pos,drop=FALSE]
        glm.out <- glm.call(Y=glmo$y, mX.=mX., cst=glmo$offset, mname=mname)
        glmo <- glm.out$glmo
        testneg<- ifelse(length(glmo$coef)>nca,any(glmo$coef[-(1:nca)]<0),FALSE)
    }
    warn <- if (nocorrect) NULL else glm.out$warnings
    return(list(glmo=glmo, warnings=warn))
# Note : possiblement à modifier pour les fonctions openp, robustd.t et robustd.0.
# pour enlever les gamma en plus des eta négatifs 
}

# Si je révise mes fonctions openp, robustd.t et robustd.0, j'envisage d'ajouter des
# fonctions ici, comme pour le calcul des paramètres.

# J'ai pensé à ajouter une fonction pour le calcul de N et de son erreurtype pour 
# un modèle de population fermée standard. Par contre, cette fonction n'aurait eu
# que 2-3 lignes et elle n'aurait été appelée que dans 3 fonctions : closedp, 
# closedp.bc et closedp.CI. J'ai donc décidé de laisser faire. 





####################################################################################################
# Fonctions de validation

"valid.one" <- function(x,type)
{
    if(!eval(call(paste("is.",type,sep=""),x))||length(x)!=1) {
       stop("'",deparse(substitute(x)),"' must be a length-one object of type ",type, call. = FALSE)
    } 
}

"valid.dtype" <- function(dtype)
{
  if(!(dtype %in% c("hist","nbcap")) ) stop("'dtype' must be given the value \"hist\" or \"nbcap\"", call. = FALSE)
}

"valid.t" <- function(t, tmin=2, pInf=FALSE) 
{
# Validation de l'argument t
# t : argument t donné en entrée
# tmin : valeur minimale que peux prendre t (les cas particuliers sont traités ailleurs) 
# pInf : Est-ce que la valeur inf est acceptée? (TRUE pour closedp.0, closedpCI.0 et descriptive seulement)
  
### tmin n'est pas vraiment nécessaire pour l'instant. Il prend toujours sa valeur par défaut de 2.
### C'est dans valid.vm que je produis des erreurs si t est trop petit pour un modèle en particulier.
  
  if (!is.null(t)) {
    isInf <- is.infinite(t)
    if (isInf) {
      if (!pInf) stop("'t' can take the value 'Inf' only with the functions 'closedp.0', 'closedpCI.0' and 'descriptive'", call. = FALSE)
    } else {
      # un entier supérieur ou égal à tmin
      if (!(is.numeric(t) && (t %% 1) == 0 && t >= tmin && length(t) == 1))
        stop("if not Null, 't' must be an integer greater or equal to ", tmin, 
             ngettext(pInf, " or take the value 'Inf'", ""), call. = FALSE)
    }
  }
}

"valid.X" <- function(X, dfreq, dtype="hist", t=NULL, vt=t, warn=FALSE)
{
### Utilité : validation de l'argument X
# Liste des fonctions appellant valid.X
# openp et closedp.Mtb : seulement X et dfreq sont donnés
# robustd.t, robustd.0 et periodhist : seulement X, dfreq et vt sont donnés
#           (pour l'instant le format dtype="nbcap" n'est pas accepté par ces fonctions)  
# closedp.bc et descriptive : seuls vt et warn ne sont pas donnés
# closedp.internal et closedpCI.internal : seul vt n'est pas donné
  
  X <- as.matrix(X)
  
  # Validation de l'argument t et modification de sa valeur s'il prend la valeur NULL ou Inf
  if (dtype=="hist") {
    tverif <- if (dfreq) ncol(X) - 1 else ncol(X)
    if (is.null(t) || is.infinite(t)) {  
      t <- tverif
      if (t<2) stop("the data set 'X' must contain at least 2 capture occasions", call. = FALSE)
    } else {
      if(t != tverif)
        stop("'t' is not equal to the number of columns in 'X'", ngettext(dfreq," minus 1",""), call. = FALSE)
    }
    if(t!=sum(na.rm=TRUE,vt))
      stop("the number of columns in 'X' ", ngettext(dfreq,"minus 1 ",""),
      "is not equal to the total number of capture occasions (sum of the 'vt' components)", 
      call. = FALSE)
  } else { # Uniquement pour closedp.bc, descriptive, closedp.internal et closedpCI.internal
    if (is.null(t)) stop("argument 't' must be given if 'dtype' takes the value \"nbcap\"", call. = FALSE) 
    if (is.infinite(t)) { t <- if (dfreq) max(X[X[,2]!=0, 1]) else max(X) }
  }
  
  # Validation de l'argument X
  if (dtype == "hist") {
    if (any(X[, 1:t] != 1 & X[,1:t] != 0))
      stop(ngettext(dfreq, "every columns of 'X' but the last one ", "'X' "), "must contain only zeros and ones", call. = FALSE)  
  } else { # Uniquement pour closedp.bc, descriptive, closedp.internal et closedpCI.internal
    ncolverif <- if(dfreq) 2 else 1
    if (ncol(X) != ncolverif){
      error <- sprintf(ngettext(ncolverif, "'X' must have %d column",
                                           "'X' must have %d columns"), ncolverif)
      stop(error, call. = FALSE)
    }
    if (any((X[, 1] %% 1) != 0) || any(X[, 1] < 1) || any(X[, 1] > t))
      stop(ngettext(dfreq, "the first column of ",""),
           "'X' must contain only integers between 1 and ", t, call. = FALSE)
  }
  if (dfreq) {
    if (any((X[, ncol(X)] %% 1) != 0) || any(X[, ncol(X)] < 0)) 
      stop("the last column of 'X' must contain capture history frequencies, therefore non-negative integers", call. = FALSE)  
  }
  
  # Pour omettre les lignes avec uniquement des zéros s'il y en a
  zeros <- if(dfreq) apply(X[,-ncol(X),drop=FALSE],1,sum)==0 else apply(X,1,sum)==0
  X <- X[!zeros, , drop=FALSE]
  ### Message d'avis omis car en fait le modèle est ajusté avec les fréquences nulles
  ### C'est pour le bon fonctionnement du code (histfreq) que ces lignes sont omises.
  # if (sum(zeros)>0) warning("the data matrix 'X' contains cases with no capture; these are ignored", call. = FALSE)
  ###
  
  # Avertissement données problématiques    
  if (warn) {
    if (t>20) warning("There is more than 20 capture occasions. This function migth fail.\n  We suggest using the 'periodhist' function to reduce the size of your data set.", immediate.=TRUE, call. = FALSE)   
    ### Message d'avis omis je ne sais plus pourquoi...
    #temp <- if (dfreq) X[X[,t+1]!=0,] else X
    #if (any(colSums(temp)==0)) warning("There is no capture on some occasions. Results can be unstable.\n  We suggest removing from the data set the occasions without captures.", immediate.=TRUE, call. = FALSE)
    ###
  }
  
  # Sortie des résultats
  return(list(X=X,t=t))
}

"valid.t0" <- function(t0, typet, t, tmin=2, tinf=FALSE) {
  if (!typet) {    
    if(is.null(t0)) {
      t0 <- t
    } else {
      if (tinf && is.infinite(t0)) {
        t0 <- t
      } else {
        tmax <- if(tinf) Inf else t
        if ( !(is.numeric(t0) && (t0 %% 1) == 0 && t0 >= tmin && t0 <= tmax && length(t0) == 1) )
          stop("'t0' must be an integer between ", tmin, " and ", tmax, ", the number of capture occasions, inclusively", call. = FALSE)
      }
    }
  } else {
    if (!is.null(t0)) {  
      ## Cette condition peut seulement être rencontrée avec la fonction closedp.bc 
      ## pour un modèle qui n'utilise pas l'Argument t0.
      if (t0 != t) warning("the input argument 't0' could not be used with the requested model", call. = FALSE)
      t0 <- NULL
    }
  }
  return(t0)
}

"valid.vt" <- function(vt)
{
  if (length(vt)==1)
    stop("'vt' must be at least of length 2\n",
         "(to analyze data for only one primary period, please use a 'closedp' function)", 
         call. = FALSE)
  if (!is.numeric(vt) || any((vt %% 1)!=0) || any(vt<=0)) 
    stop("the 'vt' components must be positive integers", call. = FALSE)
}

"valid.vm" <- function(vm,values,vt,typet=TRUE)
{
  valuesMsg <- if(typet) values else values[-grep("t",values)]
  # Note length(vt)==1 signifie que valid.vm est appelée d'une fonction closedp ou de openp
  error <- gettext(ngettext(length(vt),"'m'","'vm' elements")," can only take one of these values: ",
                 paste(dQuote(valuesMsg), collapse=", "))
  if(is.null(vm)) stop(error, call. = FALSE)		  
  if(length(vt)==1) vm <- vm[1] else if (length(vm)==1) vm <- rep(vm,length(vt))
  if(length(vt)!=length(vm)) stop("'vm' must be of length 1 or have the same length than 'vt'", call. = FALSE)
  for (i in 1:length(vm))
  {
    if(!(vm[i] %in% values)) stop(error, call. = FALSE)
    if(vm[i] %in% c("Mh","Mth") && vt[i] < 3) stop("heterogeneous models require at least 3 capture occasions", call. = FALSE)
    if(!typet && vm[i] %in% c("Mt","Mth")) stop("models with a temporal effect cannot be adjusted with a '*.0' function", call. = FALSE)
  }
  if(length(vt)!=1&&(vm[1]=="none"||vm[length(vm)]=="none")) stop("the 'no model' cannot be chosen for the first or the last period", call. = FALSE)
  return(vm)
}

"valid.h" <- function(h, values, m, call)
{
  if(length(h)!=1) h <- h[1]
  msg2 <- gettext("a function or a character string taking one of these values: ", paste(dQuote(values), collapse=", "))
  if(is.null(m)) {  ## m est NULL uniquement si mX est utilisé
    if(!(is.function(h) || h %in% values || is.null(h)))
      stop("'h' must be NULL, ", msg2, call. = FALSE) 
  } else if(m %in% c("Mh", "Mth")) {
    if (is.null(h)) {
      h <- "Chao"  ## valeur par défaut  
    } else if(!(is.function(h) || h %in% values)) {
      stop("'h' must be ", msg2, call. = FALSE)
    }
  } else {
    if (!is.null(h)) {
      h <- NULL
      warning("the input argument 'h' has been ignored since the requested model is not heterogeneous", call. = FALSE)
    }
  }
  hname <- if (is.function(h)) { deparse(call[["h"]]) 
      } else if (is.null(h))   { "none"
      } else                   { h }
  return(list(h=h, hname=hname))
}

"valid.vh" <- function(vh,values,vm)
{
  pos <- which(vm %in% c("Mh", "Mth"))
  nh <- length(pos)
  if(!is.list(vh)) vh <- if(length(vh)==1) list(vh) else as.list(vh) # transforme vh en liste si ça n'en est pas une
  if (nh>0) { 
    if(!(length(vh)==nh||length(vh)==1))
      stop("'vh' must be of length 1 or of length equals to the number of heterogeneous models specified in 'vm'", call. = FALSE)
    for (i in 1:length(vh)){
      if (is.null(vh[[i]])) {
        vh[[i]] <- "Chao"  ## valeur par défaut
      } else {
        if(!(is.function(vh[[i]]) || vh[[i]] %in% values))
          stop("'vh' elements must be a function or a character string taking one of these values:", paste(dQuote(values), collapse=", "), call. = FALSE)
      }
    }
    if(length(vh) == 1) vh <- rep(vh,nh)  ## Si vh est de longueur 1, appliquer ce h à toutes les périodes avec hétérogénéité
  }
  vht<-vector("list",length=length(vm))
  vht[pos]<-vh  ## fixe à NULL la valeur de h pour les périodes sans modèle hétérogène
  return(vht)
}

"valid.vtheta" <- function(vtheta,vh)
{
        pos<-which(vh=="Poisson"|vh=="Gamma")
        nP <- length(pos)
        if (nP>0)
        {
            if(!(length(vtheta)==nP||length(vtheta)==1))
                stop("'vtheta' must be of length 1 or of length equals to the number of Poisson or Gamma heterogeneous models specified in 'vh'", call. = FALSE)
            if (!is.numeric(vtheta)) stop("the 'vtheta' elements must be numerics", call. = FALSE)
            if (length(vtheta)==1) vtheta <- rep(vtheta,nP)
            vthetat<-rep(NA,length(vh))
            vthetat[pos]<-vtheta
            return(vthetat)
        }
}

"valid.mX" <- function(mX, typet, t, t0){
#### Description : validation de mX, appelé uniquement par closedpCI.internal
  if(!is.null(mX)){
    if (class(mX)=="formula"){  ## Si mX est une formule :
      if(!typet) stop("'mX' cannot be a formula for 'closedpCI.0'", call. = FALSE)
      if(length(mX)==3)
        stop("the formula given as 'mX' argument must not contain a response variable", call. = FALSE)
      if(any(!(all.vars(mX) %in% c(".", paste("c", 1:t, sep="")))))
        stop("the only accepted variable names in the formula given as 'mX' argument are ",
             "c1 to c", t, ", which represent the capture occasions", call. = FALSE)
      # L'erreur si l'intercept est enlevé sera généré dans getmX car j'ai besoin de histpos
    } else {  ## Si mX est supposée être une matrice :
      testmX <- try(as.matrix(mX), silent=TRUE)
      if (class(testmX) == "try-error"){
        stop("cannot turn the given 'mX' argument into a matrix", call. = FALSE)  
      } else { mX <- testmX }
      nrows <- if(typet) 2^t-1 else t0
      if (nrow(mX) != nrows) stop("'mX' must have ", ngettext(typet, "2^t-1", "t0"), " rows", call. = FALSE)
      if (any(colSums(mX==1)== nrow(mX)))
        stop("'mX' must not contain a column of ones since an intercept is always added to the model", call. = FALSE)
    }
  }
  return(mX)
}

"valid.mname" <- function(mname, typet=FALSE, m, hname, theta, call){
  if (is.null(mname)) { # Si aucun mname n'a été donné, on le crée
    # partie pour décrire le modèle d'hétérogénéité
    ph <- if (hname == "none") { NULL  
          } else if (hname %in% c("Poisson","Gamma")) { paste(hname, theta, sep="")  
          } else if (hname == "Chao") { paste(hname, " (LB)", sep="")            
          } else { paste(hname, sep="") }
    # création du mname  
    if(!is.null(call[["mX"]])) {
      mname <-  if (length(call[["mX"]]) == 1) deparse(call[["mX"]]) else "mX"
      if (!is.null(ph)) mname <- gsub("\\s","", paste(mname,"+h=", ph, sep=""))
        # Je ne veux aucun espace ici, gsub sert à enlever toutes les white spaces 
    } else { 
      mname <- if(is.null(ph)) m else paste(m, ph, sep=" ")
    }

  } else { # Si un mname a été donné en entrée, on le valide
    valid.one(mname,"character")
  }  
  return(mname)
}

"valid.alpha" <- function(alpha){
  valid.one(alpha,"numeric")
  if (alpha<0|alpha>1) stop("'alpha' must be between 0 and 1", call. = FALSE)
}

"valid.theta" <- function(theta){
  valid.one(theta,"numeric")
  if (theta<=0) stop("'theta' must be a positive number", call. = FALSE)
}

"valid.initsig" <- function(initsig){
  valid.one(initsig,"numeric")
  if (initsig<=0) stop("'initsig' must be positive since it is a standard error parameter", call. = FALSE)
}

"valid.fmaxSupCL" <- function(fmaxSupCL){
  valid.one(fmaxSupCL,"numeric")
  if (fmaxSupCL<3) stop("'fmaxSupCL' must be greater or equal to 3", call. = FALSE)
}



######################
# Version légèrement modifiée de la fonction glm

glmRcapture <- function(formula, family = gaussian, data, weights,
    subset, na.action, start = NULL,
    etastart, mustart, offset,
    control = list(...),
    model = TRUE, method = "glm.fit",
    x = FALSE, y = TRUE,
    contrasts = NULL, ...)
{
  call <- match.call()
  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
          "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame")) return(mf)
  
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  
  mt <- attr(mf, "terms") # allow model.frame to have updated it
  
  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
    stop("negative weights not allowed")
  
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call(if(is.function(method)) "method" else method,
          x = X, y = Y, weights = weights, start = start,
          etastart = etastart, mustart = mustart,
          offset = offset, family = family, control = control,
          intercept = attr(mt, "intercept") > 0L))
  
  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
#  if(length(offset) && attr(mt, "intercept") > 0L) {
#    fit$null.deviance <-
#        eval(call(if(is.function(method)) "method" else method,
#                x = X[, "(Intercept)", drop=FALSE], y = Y,
#                weights = weights, offset = offset, family = family,
#                control = control, intercept = TRUE))$deviance
#  }
## Sophie 18 mai 2012 : J'ai mis ce bout en commentaire car je ne m'intéresse
## pas à la NULL deviance. Je n'utilise pas du tout cette valeur. Alors
## les messages d'erreurs générées par cet appel à glm.fit ne sont pas
## pertient pour un utilisateur de Rcapture. En omettant ce deuxième appel
## à glm.fit, j'omet aussi les warnings.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    fit$null.deviance <- NA
  }
## À la place je donne la valeur NA à fit$null.deviance  
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula,
          terms = mt, data = data,
          offset = offset, control = control, method = method,
          contrasts = attr(X, "contrasts"),
          xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}

