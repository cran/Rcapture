uifit <- function(x.closedp)
{

    #####################################################################################################################################
    # Validation de l'argument fourni en entrée
    if(class(x.closedp)!="closedp") stop("'x.closedp' must be a object of class 'closedp'")
    #####################################################################################################################################

        t <- x.closedp$t

        ifirstcap <- 0
        for (i in 1:t) { ifirstcap <- c(ifirstcap,rep(i,2^(t-i))) }
        ifirstcap <- ifirstcap[-1]

        # Valeurs observées
        uobs <- rep(0,t)
        for ( i in 1:t ) { uobs[i] <- sum(na.rm=TRUE,x.closedp$glmM0$y[ifirstcap==i]) } 
        uobs<-c(uobs,rep(NA,5))      
        
        # Modèle M0
        p <- exp(x.closedp$glmM0$coef[2])/(1+exp(x.closedp$glmM0$coef[2]))
        upredM0 <- x.closedp$results[1,1]*p*(1-p)^(0:(t+4))
        
        # Modèle Mt
        upredMt <- rep(0,t)
        for ( i in 1:t ) { upredMt[i] <- sum(na.rm=TRUE,x.closedp$glmMt$fitted.values[ifirstcap==i]) }
        upredMt <- c(upredMt,rep(NA,5))       

        # Modèle Mh Chao
        upredMhC <- rep(0,t)
        for ( i in 1:t ) { upredMhC[i] <- sum(na.rm=TRUE,x.closedp$glmMhC$fitted.values[ifirstcap==i]) }
        upredMhC <- c(upredMhC,rep(NA,5))       

        # Modèle Mh Poisson
        EprobaP_general <- function(i,beta,tau,a,t,k){ (exp(beta)*(1+exp(beta)*a^i)^(t-k))*(a*tau)^i/(factorial(i)*sum(na.rm=TRUE,choose(t,0:t)*exp(beta*(0:t)+tau*a^(0:t))))}
        value_Eproba <- rep(0,t+5)
        for (j in 1:(t+5))
        {
            EprobaP <- function(i){ EprobaP_general(i,x.closedp$glmMhP$coef[2],x.closedp$glmMhP$coef[3],2,t,j)}
            value_Eproba[j] <- sum(na.rm=TRUE,EprobaP(0:100))       
        }
        upredMhP <- x.closedp$results[4,1]*value_Eproba

        # Modèle Mh Darroch
        if(x.closedp$glmMhD$coef[3]>0)
        {
            EprobaD_general <- function(x,beta,tau,t,k){ exp(-(x^2)/(2*tau))*((1+exp(beta+x))^(t-k))*exp(beta+x)/(sqrt(2*pi*tau)*sum(na.rm=TRUE,choose(t,0:t)*exp(beta*(0:t)+tau*((0:t)^2)/2))) }
            value_Eproba <- rep(0,t+5)
            for (i in 1:(t+5))
            {
                EprobaD <- function(x){EprobaD_general(x,x.closedp$glmMhD$coef[2],x.closedp$glmMhD$coef[3],t,i)}
                value_Eproba[i] <- integrate(EprobaD,-100,100)$value
            }       
            upredMhD <- x.closedp$results[5,1]*value_Eproba
        } else {
            upredMhD <- rep(0,t)
            for ( i in 1:t ) { upredMhD[i] <- sum(na.rm=TRUE,x.closedp$glmMhD$fitted.values[ifirstcap==i]) }
            upredMhD <- c(upredMhD,rep(NA,5))       
        }        
   
        # Modèle Mth Chao
        upredMthC <- rep(0,t)
        for ( i in 1:t ) { upredMthC[i] <- sum(na.rm=TRUE,x.closedp$glmMthC$fitted.values[ifirstcap==i]) }
        upredMthC <- c(upredMthC,rep(NA,5))       

        # Modèle Mth Poisson
        upredMthP <- rep(0,t)
        for ( i in 1:t ) { upredMthP[i] <- sum(na.rm=TRUE,x.closedp$glmMthP$fitted.values[ifirstcap==i]) }
        upredMthP <- c(upredMthP,rep(NA,5))       

        # Modèle Mth Darroch
        upredMthD <- rep(0,t)
        for ( i in 1:t ) { upredMthD[i] <- sum(na.rm=TRUE,x.closedp$glmMthD$fitted.values[ifirstcap==i]) }
        upredMthD <- c(upredMthD,rep(NA,5))       

        # Modèle Mb
        p <- 1-exp(x.closedp$glmMb$coef[2])/(1+exp(x.closedp$glmMb$coef[3]))
        upredMb <- x.closedp$results[9,1]*p*(1-p)^(0:(t+4))
               
        # Modèle Mbh
        upredMbh <- rep(0,t)
        for ( i in 1:t ) { upredMbh[i] <- sum(na.rm=TRUE,x.closedp$glmMbh$fitted.values[ifirstcap==i]) }
        upredMbh <- c(upredMbh,rep(NA,5))       

        # Stat d'ajustement du chi-deux
        statM0 <- sum(na.rm=TRUE,((uobs[1:t]-upredM0[1:t])^2)/upredM0[1:t])
        statMt <- sum(na.rm=TRUE,((uobs[1:t]-upredMt[1:t])^2)/upredMt[1:t])
        statMhC <- sum(na.rm=TRUE,((uobs[1:t]-upredMhC[1:t])^2)/upredMhC[1:t])
        statMhP <- sum(na.rm=TRUE,((uobs[1:t]-upredMhP[1:t])^2)/upredMhP[1:t])
        statMhD <- sum(na.rm=TRUE,((uobs[1:t]-upredMhD[1:t])^2)/upredMhD[1:t])
        statMthC <- sum(na.rm=TRUE,((uobs[1:t]-upredMthC[1:t])^2)/upredMthC[1:t])
        statMthP <- sum(na.rm=TRUE,((uobs[1:t]-upredMthP[1:t])^2)/upredMthP[1:t])
        statMthD <- sum(na.rm=TRUE,((uobs[1:t]-upredMthD[1:t])^2)/upredMthD[1:t])
        statMb <- sum(na.rm=TRUE,((uobs[1:t]-upredMb[1:t])^2)/upredMb[1:t])
        statMbh <- sum(na.rm=TRUE,((uobs[1:t]-upredMbh[1:t])^2)/upredMbh[1:t])
             
        # Présentation des résultats
        tableau <- cbind(uobs,upredM0,upredMt,upredMhC,upredMhP,upredMhD,upredMthC,upredMthP,upredMthD,upredMb,upredMbh)
        unames<-rep(0,t+5)
        for (i in 1:(t+5)) { unames[i] <- paste("u",i,sep = "") }
        dimnames(tableau) <- list(unames,c("observed","M0","Mt","Mh Chao","Mh Poisson2","Mh Darroch","Mth Chao","Mth Poisson2","Mth Darroch","Mb","Mbh"))
        
        stat<-matrix(c(statM0, statMt, statMhC, statMhP, statMhD, statMthC, statMthP, statMthD, statMb, statMbh),ncol=1) 
        dimnames(stat) <- list(c("M0","Mt","Mh Chao","Mh Poisson2","Mh Darroch","Mth Chao","Mth Poisson2","Mth Darroch","Mb","Mbh"),"Chi-square value")

        # Statistiques sur le jour de la première capture
        Mean <- colSums((1:t)*tableau[1:t,])/colSums(tableau[1:t,])
        Variance <- colSums(((1:t)^2)*tableau[1:t,])/colSums(tableau[1:t,]) - Mean^2
        firstcapt <- cbind(Mean,Variance)
        
        # output
        list(predicted=tableau,fit.stat=stat,day.first.capt=firstcapt)        
}
