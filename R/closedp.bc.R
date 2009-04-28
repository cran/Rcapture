"closedp.bc" <- function(X,dfreq=FALSE,dtype=c("hist","nbcap"),t,t0=t,m=c("M0","Mt","Mh","Mth","Mb","Mbh"), h=c("Chao","Poisson","Darroch","Gamma"), theta=2)
{
##### Validation des arguments fournis en entrée
    valid.one(dfreq,"logical")
    dtype <- dtype[1]
    if (!dtype%in%c("hist","nbcap")) stop("'dtype' must be given the value 'hist' or 'nbcap'")
    if (missing(t)) t <- NULL else valid.one(t,"numeric")
    Xvalid<-valid.X(X,dfreq,dtype,t,warn=FALSE)
         X <- Xvalid$X
         t <- Xvalid$t
         n <- sum(histfreq.0(X=X,dfreq=dfreq,dtype=dtype,vt=t)) 
    if (t<2) stop("the number of capture occasions 't' must be at least 2")    
    valid.one(t0,"numeric")
    if ((t0%%1)!=0||t0>t||t0<2) 
         stop("'t0' must be an integer between t, the number of capture occasion, and 2 inclusively")
    m<-valid.vm(m,c("M0","Mt","Mh","Mth","Mb","Mbh"),t)
    hname <- deparse(substitute(h))
    h<-valid.vh(h,c("Chao","Poisson","Darroch","Gamma"),m)[[1]]
    hname <- if(is.function(h)) hname else h
    valid.one(theta,"numeric")

    # Certains modèles ne peuvent être ajustés par closedp.bc
    if (dtype=="nbcap" && m%in%c("Mt","Mth","Mb","Mbh")) stop("'X' can not be of type 'nbcap' for models 'Mt', 'Mth', 'Mb' and 'Mbh'") 
    if(m=="Mbh"&&t<4) stop("The biais correction cannot be performed for the Mbh model when there is less than 4 capture occasions")
    if(m%in%c("Mh","Mth","Mb")&&t<3) 
          stop("the biais correction cannot be performed for models Mh, Mth or Mb when there is less than 3 capture occasions")
    if(m=="Mth" && !hname%in%c("Chao","Poisson","Darroch"))
          stop("the biais correction cannot be performed for models Mth when 'h' is Gamma or a function")
    if(m=="Mth" && hname=="Poisson" && theta!=2)
          stop("the biais correction can be performed for the Mth Poisson model only with theta equals to 2")
          
    if(m=="Mth" && t>20) warning("Warning: There is more than 20 capture occasions and a Mth model was requested. This function migth fail.\n  We suggest using the 'periodhist' function to reduce the size of your data set.\n",immediate.=TRUE)


##### Calculs exacts   
    if (m=="Mt") {
          converge <- NULL
          ni <- getni(X=X,dfreq=dfreq,t=t)
          if (t==2) { # Seber p.60
               m2 <- ni[1] + ni[2] - n
               N <- (ni[1]+1)*(ni[2]+1)/(m2+1) - 1
               erreurtype <- sqrt((ni[1]+1)*(ni[2]+1)*(ni[1]-m2)*(ni[2]-m2)/((m2+2)*(m2+1)^2))
          } else { # Seber p.131 et 133 avec correction des fréquences selon Rivest & Lévesque 2001
               Nprelim <- suppressWarnings(closedpCI.0(X=X,dfreq=dfreq,m="M0")$results[1,1])
               ni.corr <- ni + (t+2)/(2*t)
               n.corr <- n + t/2
               eqn <- function(N) prod(1-ni.corr/N)-(1-n.corr/N)
               eqnsolve <- uniroot(eqn,c(n,2*Nprelim))
               N.corr <- eqnsolve$root 
               N <- N.corr - t/2
               erreurtype <- sqrt((1/(N.corr-n.corr) + (t-1)/N.corr - sum(1/(N.corr-ni.corr)))^(-1))
          ######## Sans correction pour le biais ########
               ### eqn <- function(N) prod(1-ni/N)-(1-n/N)
               ### eqnsolve <- uniroot(eqn,c(n,2*Nprelim))
               ### N <- eqnsolve$root 
               ### erreurtype <- sqrt((1/(N-n) + (t-1)/N - sum(1/(N-ni)))^(-1))
          }
    } else if (m=="Mh" && hname=="Chao") {
          converge <- NULL
          Y <- histfreq.0(X=X,dfreq=dfreq,dtype=dtype,vt=t)
          fi <- rev(Y)
          N <- n + ((t - 1) * fi[1] * (fi[1] - 1))/(2 * t * (fi[2] + 1))
          erreurtype <- sqrt(((t - 1) * fi[1] * (fi[1] - 1))/(2 * t * (fi[2] + 1)) + 
                             ((t - 1)^2 * fi[1] * (fi[1] - 1) * (fi[1]^2 + 4 * fi[1] * fi[2] + 3 * fi[1] - 6 * fi[2] - 6))/
                              (4 * t^2 * (fi[2] +1)^2 * (fi[2] + 2)))
    } else {
  
     ##### Calculs avec glm
     
          # Élaboration préliminaire du modèle
          if (m%in%c("Mb","Mbh")) {
               Y <- getui(X=X,dfreq=dfreq,t=t)
               if(m=="Mbh") { omitY <- Y[1]; Y <- Y[-1]} 
               mX. <- if(m=="Mb") c(0:(t-1)) else c(0:(t-2))
               cst <- rep(0,length(Y))
          } else {
               if (m%in%c("M0","Mh")) {
                    eval(model.0)
               } else if (m=="Mth") {
                    eval(model.t)
               }
               mX. <- Xclosedp(t=t,m=m,h=h,theta=theta,histpos=histpos)$mat
               mX. <- cbind(mX.,mXomit)
          }
             
          # Ajustement des fréquences
          delta <- rep(0,length(Y))
          if (m%in%c("Mb","Mbh")) {
               delta[1] <- +2
               delta[2] <- -1 
               delta[3] <- -1     
          } else {
               if (m=="M0") {
                    del <- c(-1/2,1, 0)
               } else if (m=="Mh") {
                    hfct <- if (is.function(h)) h else 
                            if (hname=="Poisson") function(x) hP(x,theta=theta) else 
                            if (hname=="Darroch") hD else 
                            if (hname=="Gamma") function(x) hG(x,theta=theta)   
                    cc<-(hfct(2)-2*hfct(1))/(hfct(3)-2*hfct(2)+hfct(1))
                    del <- c(-(1+cc)/2,1+cc,(1-cc)/2)     
               } else if (m=="Mth" && hname=="Chao") {
                    del <- c((t-2)/(2*t),2/(t*(t-1)),0)
               } else if (m=="Mth" && hname=="Poisson") {
                    del <- c((2*t-5)/(4*t),3/(t*(t-1)),3/(2*t*(t-1)*(t-2)))
               } else if (m=="Mth" && hname=="Darroch") {
                    del <- c((t-3)/(2*t),4/(t*(t-1)),0)
               }
               nbcap <- rowSums(histpos)
               delta[nbcap==1] <- del[1]                   
               delta[nbcap==2] <- del[2]                  
               delta[nbcap==3] <- del[3]      
          }
          Yd <- pmax(0,Y + delta)
     
          # Ajustement du modèle
          glmo <- suppressWarnings(glm(Yd ~ offset(cst) + mX.,family=poisson))
                 # on omet les warnings car ils sont nombreux étant donné que les valeurs de Y ne sont pas entières.
          converge <- glmo$converge
          if(!converge) warning("algorithm did not converge")

          # Calcul de N et de son erreur type
          varcov <- summary(glmo)$cov.unscaled
          if (m%in%c("Mb","Mbh")) {
                N <- exp(glmo$coef[1])/(1-exp(glmo$coef[2]))
                v <- N*c(1,exp(glmo$coef[2])/(1-exp(glmo$coef[2])))
                erreurtype <- sqrt((t(v)%*%varcov%*%v) - N)
                if (m=="Mbh") N <- N + omitY
           } else {
                N <- n + exp(glmo$coef[1])
                erreurtype <- sqrt(exp(glmo$coef[1])+(exp(2*glmo$coef[1]))*varcov[1,1])
           }
     
     }


    # Préparation des sorties
    mname <- if (m%in%c("M0","Mt","Mb","Mbh")) m else if (hname%in%c("Poisson","Gamma")) paste(m,paste(hname,theta,sep="")) else paste(m,hname)
    results <- matrix(c(N,erreurtype),1,2)
    dimnames(results) <- list(mname,c("abundance","stderr"))
    ans <- list(n=n,results=results,converge=converge)
    class(ans) <- "closedp.bc"
    ans
}


print.closedp.bc <- function(x, ...) {       
        cat("\nNumber of captured units:",x$n,"\n\n")
        cat("Abundance estimation with bias correction:\n")
        tableau <- round(x$results,1)       
        print.default(tableau, print.gap = 2, quote = FALSE, right=TRUE, ...)
        if (!is.null(x$converge)) { if (!x$converge) cat(paste("\n ** : The model did not converge")) }
        cat("\n")
        invisible(x)
}
