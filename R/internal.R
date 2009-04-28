###################################################################################################
# Fonctions pour des bouts de code apparaissant plus d'une fois dans mes fonctions

"glm.call"<-function(Y,mX.,cst){
     for (i in 1:dim(mX.)[2]) assign(colnames(mX.)[i],mX.[,i])
     glmo <- if (all(cst==0)||is.null(cst)) {
          paste("glm(Y ~ ",paste(colnames(mX.),collapse=" + "),", family=poisson)",sep="")
     } else {
          paste("glm(Y ~ offset(cst) + ",paste(colnames(mX.),collapse=" + "),", family=poisson)",sep="")
     }
     glmo <- suppressWarnings(eval(parse(text=glmo)[[1]]))
     if(!glmo$converge) warning("algorithm did not converge")
     return(glmo)     
}


"Chao.neg" <- function(glmo,nca)
{
    testneg<-any(glmo$coef[-(1:nca)]<0)
    while(testneg) {# Répéter la boucle jusqu'à ce qu'aucun eta ne soit négatif
        # Détermination de la colonne à enlever dans mX
        pos<-nca
        while(glmo$coef[pos+1]>0) pos <- pos + 1
        # Retrait de la bonne colonne de mX et réajustement du modèle
        mX. <- if(is.null(glmo$offset)) glmo$model[,-1] else glmo$model[,-(1:2)]
        mX. <- mX.[,-pos,drop=FALSE]
        glmo <- glm.call(glmo$y,mX.,glmo$offset)                   
        testneg<- ifelse(length(glmo$coef)>nca,any(glmo$coef[-(1:nca)]<0),FALSE)
    }
    return(glmo)
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
       error <- paste("'",deparse(substitute(x)),"' must be a length-one object of type ",type,sep="")
       stop(error)
    } 
}

"valid.X" <- function(X,dfreq,dtype="hist",t=NULL,warn=FALSE)
{
    X <- as.matrix(X)
    if (dtype=="hist") { # cas dtype="hist"         
         if (dfreq) {
             if (is.null(t)) {
                    t <- dim(X)[2]-1
             } else {
                    if(t!=dim(X)[2]-1) stop("'t' is not equal to the number of columns in X minus 1")
             }        
             if (any(X[,1:t]!=1&X[,1:t]!=0)) stop("every columns of 'X' but the last one must contain only zeros and ones")
             if (any((X[,t+1]%%1)!=0)||any(X[,t+1]<0)) 
                    stop("the last column of 'X' must contain capture histories frequencies, therefore non-negative integers")
         } else {
             if (is.null(t)) {
                    t <- dim(X)[2]
             } else {
                    if(t!=dim(X)[2]) stop("'t' is not equal to the number of columns in X")
             }        
             if(any(X!=1&X!=0)) stop("'X' must contain only zeros and ones")
         }
    } else {  # cas dtype="nbcap", pour le moment valide seulement pour populations fermées
         if (is.null(t)) stop("argument 't' must be given if 'dtype' takes the value 'nbcap'") 
         if (dfreq) {
             if (dim(X)[2]!=2) stop("X must have two columns") 
             if (any((X[,1]%%1)!=0)||any(X[,1]<0)||any(X[,1]>t)) stop("the first column of 'X' must contain only integers between 0 and 't'")
             if (any((X[,2]%%1)!=0)||any(X[,2]<0))
                    stop("the last column of 'X' must contain capture histories frequencies, therefore non-negative integers")
         } else {
             if (dim(X)[2]!=1) stop("X must have only one column") 
             if (any((X%%1)!=0)||any(X<0)||any(X>t)) stop("'X' must contain only integers between 0 and 't'")
         }
    }
    
    ## Pour omettre les lignes avec uniquement des zéros s'il y en a
    zeros <- if(dfreq) apply(X[,-dim(X)[2],drop=FALSE],1,sum)==0 else apply(X,1,sum)==0
    # if (sum(zeros)>0) warning("the data matrix 'X' contains cases with no capture; these are ignored")
    X <- X[!zeros,,drop=FALSE]
    
    if (warn) {
        # Avertissement données problématiques    
        if (t>20) warning("Warning: There is more than 20 capture occasions. This function migth fail.\n  We suggest using the 'periodhist' function to reduce the size of your data set.\n",immediate.=TRUE)
        
#        temp <- if (dfreq) X[X[,t+1]!=0,] else X
#        if (any(colSums(temp)==0)) warning("Warning: There is no capture on some occasions. Results can be unstable.\n  We suggest removing from the data set the occasions without captures.\n",immediate.=TRUE)
    }

    return(list(X=X,t=t))
}

"valid.vt" <- function(vt,t)
{
    if(t!=sum(na.rm=TRUE,vt)) stop("the number of columns in 'X', minus 1 if 'dfreq' is TRUE, is not equal to the total number of capture occasions (sum of the 'vt' components)")
    if (any((vt %% 1)!=0)) stop("the 'vt' components must be integers")
}

"valid.vm" <- function(vm,values,vt,typet=TRUE)
{
    if(length(vt)==1) vm <- vm[1] else if (length(vm)==1) vm <- rep(vm,length(vt))
    if(length(vt)!=length(vm)) stop("'vm' must be of length 1 or have the same length than 'vt'")
    for (i in 1:length(vm))
    {
        if(!vm[i]%in%values){
          error <- paste(ifelse(length(vt)==1,"'m'","'vm' elements"),"can only take one of these values: ",paste("'",values,"'",sep="",collapse=" ")) 
          stop(error)
        }
        if(vm[i]=="Mh"&&vt[i]<3) stop("Mh models require at least 3 capture occasions")
        if(vm[i]=="Mth"&&vt[i]<3) stop("Mth models require at least 3 capture occasions")
        if(!typet&&vm[i]%in%c("Mt","Mth")) stop("models with a temporal effect can not be adjusted")
    }
    if(length(vt)!=1&&(vm[1]=="none"||vm[length(vm)]=="none")) stop("the 'no model' cannot be chosen for the first or the last period")
    return(vm)
}

"valid.vh" <- function(vh,values,vm)
{
        pos<-which(vm=="Mh"|vm=="Mth")
        nh <- length(pos)
        if(!is.list(vh)) vh <- if(length(vh)==1) list(vh) else as.list(vh) # transforme vh en liste si ça n'en est pas une
        if (nh>0)
        { 
            if(length(vm)==1&&length(vh)!=1) vh <- list(vh[1]) else if(length(vm)!=1&&length(vh)==1) vh <- rep(vh,nh)
            if(length(vm)!=1&&!(length(vh)==nh||length(vh)==1))
                stop("'vh' must be of length 1 or of length equals the number of heterogenous models specified in 'vm'")
            for (i in 1:length(vh)){
                if(!is.function(vh[[i]])) if(!vh[[i]]%in%values) {
                    error <- paste(ifelse(length(vm)==1,"'h'","'vh' elements"),"must be a function or a charater strings taking one of these values:",
                             paste("'",values,"'",sep="",collapse=" ")) 
                    stop(error)
                }     
            }
        }
        vht<-vector("list",length=length(vm))
        vht[pos]<-vh
        return(vht)
}

"valid.vtheta" <- function(vtheta,vh)
{
        pos<-which(vh=="Poisson"|vh=="Gamma")
        nP <- length(pos)
        if (nP>0)
        {
            if(!(length(vtheta)==nP||length(vtheta)==1))
                stop("'vtheta' must be of length 1 or of length equals the number of Poisson or Gamma heterogenous models specified in 'vh'")
            if (!is.numeric(vtheta)) stop("the 'vtheta' elements must be numerics")
            if (length(vtheta)==1) vtheta <- rep(vtheta,nP)
            vthetat<-rep(NA,length(vh))
            vthetat[pos]<-vtheta
            return(vthetat)
        }
}


####################################################################################################
# Bouts de code communs à plus d'une fonction
valid.0 <- quote({
     typet<-FALSE
     dtype <- dtype[1]
     if(!dtype%in%c("hist","nbcap")) stop("'dtype' must be given the value 'hist' or 'nbcap'")
     if (missing(t)) t <- NULL else valid.one(t,"numeric")
})

valid.t <- quote({
     typet<-TRUE
     dtype<-"hist"
     t <- NULL
})

model.0 <- quote({ # Construction du modèle de type .t
     Y <- histfreq.0(X=X,dfreq=dfreq,dtype=dtype,vt=t)
     histpos <- histpos.0(t)
     cst <- log(choose(t, t:1))
     if (t0!=t) {
          mXomit <- rbind(diag(1,t-t0),matrix(0,t0,t-t0))
          mXomit <- mXomit[,(t-t0):1] # Pour renverser l'ordre des colonnes
          colnames(mXomit) <- paste("omit",(t0+1):t,sep="")
     } else mXomit <- NULL
})

model.t <- quote({ # Construction du modèle de type .0
     Y <- histfreq.t(X,dfreq=dfreq)
     histpos <- histpos.t(t)
     cst <- rep(0,length(Y))
     mXomit <- NULL
})


