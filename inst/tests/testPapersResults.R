context("Paper results")

test_that("Rcapture reproduces Rivest (2011) results for the R3 data of Darroch et al. (1993)", {
      histpos <- histpos.t(3)
      DarR3<-cbind(histpos, c(72, 155, 7, 71, 13, 53, 43))

      # saturated loglinear model
      sat <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=~(c1+c2+c3)^2,mname="saturated")
      expect_that(round(sat$results[,"abundance"]), equals(1240))
      expect_that(round(sat$results[,"stderr"]), equals(451))

      # removing the [E1,E3] interaction
      inter2 <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=~.+c1:c2+c2:c3,mname="inter2")
      expect_that(round(inter2$results[,"abundance"]), equals(850))
      expect_that(round(inter2$results[,"stderr"]), equals(186))
      expect_that(round(inter2$results[,"deviance"],1), equals(3.8))
      expect_that(inter2$results[,"df"], equals(1))
      
      # Best fitting LB model : [E1,E2] interaction only
      LB <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=~.+c1:c2,h="LB",mname="LB")
      expect_that(round(LB$results[,"abundance"]), equals(681))
      expect_that(round(LB$results[,"stderr"]), equals(78))
      expect_that(round(LB$results[,"deviance"],2), equals(3.45))
      expect_that(LB$results[,"df"], equals(1))
      
      # normal model with [E1,E2] interaction only
      nor <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=~.+c1:c2,h="Normal",mname="normal")      
      expect_that(round(nor$results[,"abundance"]), equals(1400))
      expect_that(round(nor$results[,"stderr"]), equals(559))
      # There was an error in the paper

      # Darroch model with [E1,E2] interaction only
      Dar <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=~.+c1:c2,h="Darroch",mname="Darroch")      
      expect_that(round(Dar$results[,"abundance"],1), equals(1181.8))
      expect_that(round(Dar$results[,"stderr"],1), equals(404.2))
      # There was an error in the paper
})


test_that("Rcapture reproduces Rivest (2011) results for the diabetes data of Bruno et al. (1994)", {
      ### Total data set
      histpos <- histpos.t(4)
      diabetes<-cbind(histpos,c(58,157,18,104,46,650,12,709,14,20,7,74,8,182,10))

      # saturated loglinear model
      sat <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~(c1+c2+c3+c4)^2,mname="saturated")
      inter <- sat$glm$coefficients[6:11]      
      expect_that(all(inter > 0.15 & inter < 1.94),is_true())
      
      # Bruno et al. (1994) model
      Bruno <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~(c1+c2+c3+c4)^2-c1:c4,mname="Bruno")
      expect_that(round(Bruno$results[,"abundance"]), equals(2771))
      expect_that(round(Bruno$results[,"stderr"]), equals(145))
      expect_that(round(Bruno$results[,"deviance"],2), equals(7.62))
      expect_that(Bruno$results[,"df"], equals(5))
      
      # chosen interaction set in Rivest (2001)
      nohetero <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,mname="nohetero")
      expect_that(round(nohetero$results[,"abundance"]), equals(2472))
      expect_that(round(nohetero$results[,"stderr"]), equals(53))
      expect_that(round(nohetero$results[,"deviance"],2), equals(21.97))
      expect_that(nohetero$results[,"df"], equals(7))
      
      # LB model with the chosen interaction set
      LB <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="LB",mname="LB")
      expect_that(round(LB$results[,"abundance"]), equals(2588))
      expect_that(round(LB$results[,"stderr"]), equals(75))
      expect_that(round(LB$results[,"deviance"],2), equals(3.13))
      expect_that(LB$results[,"df"], equals(6))

      # Poisson2 model with the chosen interaction set
      Pois2 <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="Poisson",mname="Poisson2")
      expect_that(round(Pois2$results[,"abundance"]), equals(2573))
      expect_that(round(Pois2$results[,"stderr"]), equals(76))
      expect_that(round(Pois2$results[,"deviance"],2), equals(12.88))
      expect_that(Pois2$results[,"df"], equals(6))
      
      # Darroch model with the chosen interaction set
      Dar <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="Darroch",mname="Darroch")
      expect_that(round(Dar$results[,"abundance"]), equals(2752))
      expect_that(round(Dar$results[,"stderr"]), equals(133))
      expect_that(round(Dar$results[,"deviance"],2), equals(8.32))
      expect_that(Dar$results[,"df"], equals(6))
      
      # normal model with the chosen interaction set
      nor <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="Normal",mname="normal")
      expect_that(round(nor$results[,"abundance"]), equals(2764))
      expect_that(round(nor$results[,"stderr"]), equals(101))
      expect_that(round(nor$results[,"deviance"],2), equals(8.22))
      expect_that(nor$results[,"df"], equals(6))
      
      # Gamma3.5 model with the chosen interaction set
      Gam3.5 <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="Gamma",mname="Gamma3.5")
      expect_that(round(Gam3.5$results[,"abundance"]), equals(2964))
      expect_that(round(Gam3.5$results[,"stderr"]), equals(219))
      expect_that(round(Gam3.5$results[,"deviance"],2), equals(6.74))
      expect_that(Gam3.5$results[,"df"], equals(6))

      # Gamma0.5 model with the chosen interaction set
      Gam0.5 <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=~.+c1:c3+c2:c4+c3:c4,h="Gamma",h.control=list(theta=0.5),mname="Gamma0.5")
      expect_that(round(Gam0.5$results[,"abundance"]), equals(4238))
      expect_that(round(Gam0.5$results[,"stderr"]), equals(952))
      expect_that(round(Gam0.5$results[,"deviance"],2), equals(5.49))
      expect_that(Gam0.5$results[,"df"], equals(6))
      
      ### Pourrait aussi être fait pour les 3 autres jeux de données      
}) 


test_that("Rcapture reproduces Coull and Agresti (1999) results for the snowshoe hare data set", {
      data(hare)
      # Mth Normal
      fit <- closedpCI.t(X=hare,m="Mth",h="Normal")     
      # Table 2, logistic-normal model
      expect_that(round(fit$results[,"abundance"],digits=1), equals(92.0))
      expect_that(round(fit$fit$parameters["estimate","sigma"],digits=2), equals(0.97))
      
#      hi <- histpos.t(6)
#      reord <- order(hi[,6], hi[,5], hi[,4], hi[,3], hi[,2], hi[,1], decreasing = TRUE)
#      expect_that(round(as.vector(fit$fit$fitted.values[reord]),digits=1), equals(
#        c(0.7,0.7,0.3,0.5,0.5,0.9,0.4,1.0,  0.4,0.6,0.3,0.7,0.4,1.1,0.5,1.9,  0.4,0.7,0.3,0.8,0.5,1.3,0.6,2.3, 
#          0.4,0.9,0.4,1.5,0.6,2.6,1.1,6.8,  0.3,0.4,0.2,0.5,0.3,0.8,0.3,1.3,  0.2,0.5,0.2,0.9,0.4,1.5,0.6,3.9,
#          0.3,0.6,0.3,1.1,0.5,1.8,0.8,4.8,  0.3,1.2,0.5,3.2,0.9,5.4,2.3)))
### Seules deux valeurs prédites ne sont pas égales (éléments 29 et 49 des vecteurs ci-dessus).      

      # Mth Darroch
      fit2 <- closedpCI.t(X=hare,m="Mth",h="Darroch")     
      # Table 2, Homogeneous two-factor association
      expect_that(round(fit2$results[,"abundance"],digits=1), equals(90.5))
    })


test_that("Rcapture reproduces Cormack, Chang and Smith (2000) results in Tables 3a and 3b", {
      # Table 3a : Données avec 4 listes
      cor.dat<-cbind(histpos.t(4),c(72,5,29,53,24,3,12,14,0,0,0,1,0,0,2))
      #liste 1= death liste 2= NTOF  liste 3=WCC liste 4=MOSH
      
      a1 <- closedpCI.t(cor.dat, dfreq = TRUE, m="Mt")
      expect_that(round(a1$results[,"abundance"]), equals(215))
      expect_that(round(a1$results[,"deviance"],digits=1), equals(87.1))
      expect_that(a1$results[,"df"], equals(10))
      
      a2 <- closedpCI.t(cor.dat, dfreq = TRUE, mX = ~ (c1+c2+c3+c4)^2)
      expect_that(a2$results[,"abundance"] > 1000000, is_true())
      expect_that(round(a2$results[,"deviance"],digits=1), equals(1.3))
      expect_that(a2$results[,"df"], equals(4))
      
      a3 <- closedpCI.t(cor.dat, dfreq = TRUE, mX = ~ . + c3:c4)
      expect_that(round(a3$results[,"abundance"]), equals(215))
#      expect_that(round(a3$results[,"deviance"],digits=1), equals(12.7))
      expect_that(a3$results[,"df"], equals(9))
      
      a4 <- closedpCI.t(cor.dat, dfreq = TRUE, mX = ~ . + c3:c4 + c1:c3)
#      expect_that(round(a4$results[,"abundance"]), equals(215))
      expect_that(round(a4$results[,"deviance"],digits=1), equals(7.7))
      expect_that(a4$results[,"df"], equals(8))
      
      a5 <- closedpCI.t(cor.dat, dfreq = TRUE, mX = ~ . + c1:c2 + c1:c3 + c3:c4)
      expect_that(round(a5$results[,"abundance"]), equals(218))
      expect_that(round(a5$results[,"deviance"],digits=1), equals(2.6))
      expect_that(a5$results[,"df"], equals(7))
      
      a6 <- closedpCI.t(cor.dat, dfreq = TRUE, mX = ~ . + c1:c2 + c1:c3 + c1:c4 + c3:c4)
      expect_that(a6$results[,"abundance"] > 1000000, is_true())
#      expect_that(round(a6$results[,"deviance"],digits=1), equals(1.4))
      expect_that(a6$results[,"df"], equals(6))
      
      
      
      ### Rcapture ne donne pas les mêmes bornes d'intervalle que celles de l'article.
      
      # Table 3b : Données avec 3 listes
      cor.dat3<-cbind(histpos.t(3),c(72,5,29,54,24,3,14))

      b1 <- closedpCI.t(cor.dat3, dfreq = TRUE, m="Mt")
      expect_that(round(b1$results[,"abundance"]), equals(209))
      expect_that(round(b1$results[,"deviance"],digits=1), equals(75.1))
      expect_that(b1$results[,"df"], equals(3))
      expect_that(round(b1$CI[,"InfCL"]), equals(203))
      expect_that(round(b1$CI[,"SupCL"]), equals(217))
      
      b2 <- closedpCI.t(cor.dat3, dfreq = TRUE, mX= ~ . + c2:c3)
#      expect_that(round(b2$results[,"abundance"]), equals(221))
      expect_that(round(b2$results[,"deviance"],digits=1), equals(1.2))
      expect_that(b2$results[,"df"], equals(2))
      expect_that(round(b2$CI[,"InfCL"]), equals(210))
#      expect_that(round(b2$CI[,"SupCL"]), equals(237))
      
      b3 <- closedpCI.t(cor.dat3, dfreq = TRUE, mX= ~ . + c1:c2 + c2:c3)
#      expect_that(round(b3$results[,"abundance"]), equals(223))
      expect_that(round(b3$results[,"deviance"],digits=1), equals(0.6))
      expect_that(b3$results[,"df"], equals(1))
#      expect_that(round(b3$CI[,"InfCL"]), equals(208))
#      expect_that(round(b3$CI[,"SupCL"]), equals(246))
      
      b4 <- closedpCI.t(cor.dat3, dfreq = TRUE, mX= ~ . + c1:c3 + c2:c3)
      expect_that(round(b4$results[,"abundance"]), equals(233))
#      expect_that(round(b4$results[,"deviance"],digits=1), equals(0.3))
      expect_that(b4$results[,"df"], equals(1))
      expect_that(round(b4$CI[,"InfCL"]), equals(204))
#      expect_that(round(b4$CI[,"SupCL"]), equals(313))
      
      b5 <- closedpCI.t(cor.dat3, dfreq = TRUE, mX= ~ (c1+c2+c3)^2)
#      expect_that(round(b5$results[,"abundance"]), equals(247))
      expect_that(round(b5$results[,"deviance"],digits=1), equals(0))
      expect_that(b5$results[,"df"], equals(0))
#      expect_that(round(b5$CI[,"InfCL"]), equals(204))
#      expect_that(round(b5$CI[,"SupCL"]), equals(395))
      
      
# J'arrive à des résultats très proches mais pas égaux à ceux de l'article pour les lignes en commentaires

# Je me demande pourquoi les IC que donne closedpCI.t sont si loin de ceux de la table 3a alors
# qu'ils sont pratiquement égaux à ceux de la table 3b.
    })

    
test_that("Rcapture reproduces Dorazio and Royle (2003) results in Table 3 for the normal model", {
      data(hare)
      xMhN <- closedpCI.t(hare,  m="Mh", h="Normal")
      expect_that(round(xMhN$results[,"abundance"],digits=1), equals(91.7))
      expect_that(round(xMhN$results[,"deviance"],digits=1), equals(61.7))
      expect_that(xMhN$results[,"df"], equals(60))
      
      xMthN <- closedpCI.t(hare,  m="Mth", h="Normal")
#      expect_that(round(xMthN$results[,"abundance"],digits=1), equals(91.9))
#      expect_that(round(xMthN$results[,"deviance"],digits=1), equals(52.8))
      expect_that(xMthN$results[,"df"], equals(55))
# J'arrive à des résultats très proches mais pas égaux à ceux de l'article pour les lignes en commentaires
      
      
      data(BBS2001)
      xx01 <- closedpCI.0(X=BBS2001, dfreq=TRUE, dtype="nbcap", t=50, m="Mh", h="Normal")
      expect_that(round(xx01$results[,"abundance"],digits=1), equals(82.2))
      expect_that(round(xx01$results[,"deviance"],digits=1), equals(35.3))
      expect_that(xx01$results[,"df"], equals(47))
      
      BBS1997 <- cbind(1:27,c(14,10,6,3,3,6,2,3,2,1,1,1,0,3,2,0,3,1,1,1,0,0,1,2,0,0,1))
      xx97 <- closedpCI.0(X=BBS1997, dfreq=TRUE, dtype="nbcap", t=50, m="Mh", h="Normal")
      expect_that(round(xx97$results[,"abundance"],digits=1), equals(78.0))
      expect_that(round(xx97$results[,"deviance"],digits=1), equals(29.0))
      expect_that(xx97$results[,"df"], equals(47))      
      
      BBS1998 <- cbind(1:30,c(14,9,7,4,2,4,4,0,4,2,0,1,0,3,0,3,1,0,0,1,1,0,1,1,2,0,1,0,1,1))
      xx98 <- closedpCI.0(X=BBS1998, dfreq=TRUE, dtype="nbcap", t=50, m="Mh", h="Normal")
      expect_that(round(xx98$results[,"abundance"],digits=1), equals(80.0))
      expect_that(round(xx98$results[,"deviance"],digits=1), equals(37.9))
      expect_that(xx98$results[,"df"], equals(47))      
   
      BBS1999 <- cbind(1:31,c(11,12,10,4,4,1,4,2,3,3,0,2,4,1,1,0,1,1,1,2,0,1,0,0,0,0,0,0,0,1,3))
      xx99 <- closedpCI.0(X=BBS1999, dfreq=TRUE, dtype="nbcap", t=50, m="Mh", h="Normal")
      expect_that(round(xx99$results[,"abundance"],digits=1), equals(85.2))
      expect_that(round(xx99$results[,"deviance"],digits=1), equals(38.9))
      expect_that(xx99$results[,"df"], equals(47))      

      BBS2000 <- cbind(1:35,c(14,7,5,5,5,4,3,4,1,1,2,1,1,3,3,1,1,0,2,0,0,0,1,0,1,1,0,1,0,0,1,0,0,0,1))
      xx00 <- closedpCI.0(X=BBS2000, dfreq=TRUE, dtype="nbcap", t=50, m="Mh", h="Normal")
      expect_that(round(xx00$results[,"abundance"],digits=1), equals(82.0))
      expect_that(round(xx00$results[,"deviance"],digits=1), equals(24.1))
      expect_that(xx00$results[,"df"], equals(47))
    })    

    
test_that("Rcapture reproduces some of Cormack (1993) results in Tables 3 and 4", {
      data(bunting)
      m1 <- closedpCI.t(bunting, dfreq=TRUE, m="Mt")
      expect_that(round(m1$results[,"deviance"],digits=1), equals(854.2))
      expect_that(m1$results[,"df"], equals(246))
      
      m2 <- openp(bunting, dfreq=TRUE)
      expect_that(round(m2$model.fit[,"deviance"],digits=2), equals(219.41))
      expect_that(m2$model.fit[,"    df"], equals(234))
      expect_that(round(as.vector(m2$survivals[1:6,"estimate"]),digits=3), equals(c(0.488,0.223,0.551,0.324,0.408,0.333)))
      expect_that(round(as.vector(m2$capture.prob[2:7,"estimate"]),digits=3), equals(c(0.378,0.23,0.374,0.472,0.451,0.619)))
      
      m3 <- openp(bunting, dfreq=TRUE, m="ep")     
      expect_that(round(m3$model.fit[,"deviance"],digits=2), equals(233.85))
      expect_that(m3$model.fit[,"    df"], equals(239))
      
      keep2<-apply(histpos.t(8),1,sum)>1
      m5<- suppressWarnings(openp(bunting,dfreq=TRUE,keep=keep2))
      expect_that(round(m5$model.fit[,"deviance"],digits=2), equals(125.18))
      expect_that(m5$model.fit[,"    df"], equals(228))
      expect_that(round(as.vector(m5$survivals[2:6,"estimate"]),digits=3), equals(c(0.485,0.674,0.729,0.518,0.556)))
      expect_that(round(as.vector(m5$capture.prob[2:7,"estimate"]),digits=3), equals(c(0.687,0.382,0.633,0.606,0.612,0.714)))
      
    })    

test_that("Rcapture reproduces some of Royle (2006) results in Table 2", {
      data(catb)
      xx <- closedp.0(catb, dfreq = TRUE, dtype = "nbcap", t = 11)
      expect_that(round(xx$results["M0","deviance"],digits=2), equals(8.05))
      expect_that(11-xx$results["M0","df"], equals(2))
      expect_that(round(xx$parameters$M0[,"p"],digits=3), equals(0.219))
      
      MhN <- closedpCI.0(catb, dfreq = TRUE, dtype = "nbcap", t = 11, m="Mh", h="Normal")
      expect_that(round(MhN$results[,"deviance"],digits=2), equals(6.02))
      expect_that(11-MhN$results[,"df"], equals(3))
      

      JAY <- cbind(1:5,c(9,11,6,5,2))
      xx <- closedp.0(JAY, dfreq = TRUE, dtype = "nbcap", t = 11)
      expect_that(round(xx$results["M0","deviance"],digits=2), equals(1.88))
      expect_that(11-xx$results["M0","df"], equals(2))
      expect_that(round(xx$parameters$M0[,"p"],digits=3), equals(0.199))
      
      MhN <- closedpCI.0(JAY, dfreq = TRUE, dtype = "nbcap", t = 11, m="Mh", h="Normal")
      expect_that(round(MhN$results[,"deviance"],digits=2), equals(1.87))
      expect_that(11-MhN$results[,"df"], equals(3))

      
      CYT <- cbind(1:9,c(6,7,5,3,1,4,5,3,2))
      xx <- closedp.0(CYT, dfreq = TRUE, dtype = "nbcap", t = 11)
      expect_that(round(xx$results["M0","deviance"],digits=2), equals(41.22))
      expect_that(11-xx$results["M0","df"], equals(2))
      expect_that(round(xx$parameters$M0[,"p"],digits=3), equals(0.385))
      
      MhN <- closedpCI.0(CYT, dfreq = TRUE, dtype = "nbcap", t = 11, m="Mh", h="Normal")
      expect_that(round(MhN$results[,"deviance"],digits=2), equals(8.93))
      expect_that(11-MhN$results[,"df"], equals(3))
      
      
      SOSP <- cbind(1:11,c(5,1, 5, 1, 5, 3, 1, 3, 1, 0, 1))
      xx <- closedp.0(SOSP, dfreq = TRUE, dtype = "nbcap", t = 11)
      expect_that(round(xx$results["M0","deviance"],digits=2), equals(37.89))
      expect_that(11-xx$results["M0","df"], equals(2))
      expect_that(round(xx$parameters$M0[,"p"],digits=3), equals(0.419))
      
      MhN <- closedpCI.0(SOSP, dfreq = TRUE, dtype = "nbcap", t = 11, m="Mh", h="Normal")
#      expect_that(round(MhN$results[,"deviance"],digits=2), equals(11.24))
      expect_that(11-MhN$results[,"df"], equals(3))

      
      TREE <- cbind(1:10,c(7, 8, 1, 3, 3, 1, 2, 3, 0, 1))
      xx <- closedp.0(TREE, dfreq = TRUE, dtype = "nbcap", t = 11)
      expect_that(round(xx$results["M0","deviance"],digits=2), equals(43.69))
      expect_that(11-xx$results["M0","df"], equals(2))
      expect_that(round(xx$parameters$M0[,"p"],digits=3), equals(0.331))
      
      MhN <- closedpCI.0(TREE, dfreq = TRUE, dtype = "nbcap", t = 11, m="Mh", h="Normal")
#      expect_that(round(MhN$results[,"deviance"],digits=2), equals(9.70))
      expect_that(11-MhN$results[,"df"], equals(3))

# J'arrive à des résultats très proches mais pas égaux à ceux de l'article pour les lignes en commentaires
      
})    


test_that("Rcapture reproduces some of Cormack (1989) results in Tables 4, 5 and 6", {
      #### Table 4
      data(hare)
      xx <- closedp.t(hare)
      expect_that(round(xx$results["M0","abundance"],digits=1), equals(75.4))
      expect_that(round(xx$results["Mt","abundance"],digits=1), equals(75.1))
      expect_that(round(xx$results["M0","deviance"]-xx$results["Mt","deviance"],digits=1), equals(10.2))
      
      ## hare5 <- hare[rowSums(hare)<=5,]
      ## closedp.t(hare5)
      ## Cette façon de faire ne fonctionne pas car ça ne supprime pas une ligne de la matrice de design X et une donnée
      ## de Y. Ça donne plutôt la valeur 0 à la fréquence de l'historique 111111.       
      mat <- histpos.t(6)
      col <- rep(0, 2^6-1)
      col[rowSums(mat) == 6] <- 1
      M05 <- closedpCI.t(hare, mX=cbind(rowSums(mat), col), mname="M0 without 111111")
      # On obtient les mêmes résultats pour N avec closedp.0(hare, t0=5), mais pas les mêmes déviances et ddl.    
      Mt5 <- closedpCI.t(hare, mX=cbind(mat, col), mname="Mt without 111111")
      expect_that(round(M05$results[,"abundance"],digits=1), equals(77.2))
      expect_that(round(Mt5$results[,"abundance"],digits=1), equals(76.8))
      expect_that(round(M05$results[,"deviance"]-Mt5$results[,"deviance"],digits=1), equals(10.7))
      
      cols <- matrix(0, nrow=2^6-1, ncol=7)
      pos1 <- (1:(2^6-1))[rowSums(mat) >= 5]
      diag(cols[pos1,]) <- 1
      M04 <- closedpCI.t(hare, mX=cbind(rowSums(mat), cols), mname="M0 without 5 capt")
      Mt4 <- closedpCI.t(hare, mX=cbind(mat, cols), mname="Mt without 5 capt")
      expect_that(round(M04$results[,"abundance"],digits=1), equals(77.5))
      expect_that(round(Mt4$results[,"abundance"],digits=1), equals(77.1))
      expect_that(round(M04$results[,"deviance"]-Mt4$results[,"deviance"],digits=1), equals(9.8))
      
      cols <- matrix(0, nrow=2^6-1, ncol=22)
      pos1 <- (1:(2^6-1))[rowSums(mat) >= 4]
      diag(cols[pos1,]) <- 1
      M03 <- closedpCI.t(hare, mX=cbind(rowSums(mat), cols), mname="M0 without 4 capt")
      Mt3 <- closedpCI.t(hare, mX=cbind(mat, cols), mname="Mt without 4 capt")
      expect_that(round(M03$results[,"abundance"],digits=1), equals(78.4))
#      expect_that(round(Mt3$results[,"abundance"],digits=1), equals(78.2))
      expect_that(round(M03$results[,"deviance"]-Mt3$results[,"deviance"],digits=1), equals(6.3))

      
      # Pour reproduire les résultats de la première partie de la table 4.
      expect_that(round(xx$results["M0","deviance"]-M05$results[,"deviance"],digits=1), equals(9.9))
      expect_that(round(M05$results[,"deviance"]-M04$results[,"deviance"],digits=1), equals(3.7)) # résultats inversés sur la ligne dans table de Cormack
      expect_that(round(M04$results[,"deviance"]-M03$results[,"deviance"],digits=1), equals(14.1))

      expect_that(round(xx$results["Mt","deviance"]-Mt5$results[,"deviance"],digits=1), equals(10.4))
      expect_that(round(Mt5$results[,"deviance"]-Mt4$results[,"deviance"],digits=1), equals(2.8)) # résultats inversés sur la ligne dans table de Cormack
      expect_that(round(Mt4$results[,"deviance"]-Mt3$results[,"deviance"],digits=1), equals(10.6))


      #### Table 5
      data(duck)
      
      cp <- closedpCI.t(duck, dfreq=TRUE, m="Mt")
      expect_that(round(cp$results[,"deviance"],digits=1), equals(655.2))
      expect_that(cp$results[,"df"], equals(56))
      
      upF <- suppressWarnings(openp(duck, dfreq=TRUE, neg=FALSE))
      expect_that(round(upF$model.fit[,"deviance"],digits=1), equals(83.1))
      expect_that(upF$model.fit[,"    df"], equals(48))
      
      up <- openp(duck, dfreq=TRUE)
      expect_that(round(up$model.fit[,"deviance"],digits=1), equals(83.4))
      expect_that(up$model.fit[,"    df"], equals(49))
      expect_that(round(as.vector(up$N[2:5,"estimate"])), equals(c(379,455,357,456)))
#      expect_that(round(as.vector(up$capture.prob[2:5,"estimate"]),digits=2), equals(c(0.47,0.46,0.45,0.42)))
      expect_that(round(as.vector(up$survivals[2:4,"estimate"]),digits=2), equals(c(0.96,0.71,1)))
      expect_that(round(as.vector(up$birth[2:4,"estimate"])), equals(c(93,32,99)))
      
      epF <- suppressWarnings(openp(duck, dfreq=TRUE, m="ep", neg=FALSE))
      expect_that(round(epF$model.fit[,"deviance"],digits=1), equals(84.8))
      expect_that(epF$model.fit[,"    df"], equals(51))
      expect_that(round(as.vector(epF$N[2:5,"estimate"])), equals(c(385,464,354,441)))
      expect_that(round(as.vector(epF$capture.prob[2:5,"estimate"]),digits=2), equals(rep(0.45,4)))
      expect_that(round(as.vector(epF$survivals[2:4,"estimate"]),digits=2), equals(c(0.96,0.70,0.99)))
      expect_that(round(as.vector(epF$birth[2:4,"estimate"])), equals(c(93,29,90)))
      
      
      #### Table 6
      data(duck)
      
      all <- openp(duck, dfreq=TRUE)
      expect_that(round(all$model.fit[,"deviance"],digits=1), equals(83.4))
      expect_that(all$model.fit[,"    df"], equals(49))
      expect_that(round(as.vector(all$N[2:5,"estimate"])), equals(c(379,455,357,456)))
      pos <- c(1,33,17,9,5,3,2)
      expect_that(round(as.vector(fitted.values(all$glm)[pos]),digits=1), equals(c(3.9,1.7,4.5,4.5,4.8,5.4,3.4)))
      
      keepNo6 <- rowSums(histpos.t(6)) != 6
      No6 <- openp(duck, dfreq=TRUE, keep=keepNo6)
      expect_that(round(No6$model.fit[,"deviance"],digits=1), equals(67.3))
      expect_that(No6$model.fit[,"    df"], equals(48))
      expect_that(round(as.vector(No6$N[2:5,"estimate"])), equals(c(388,470,369,471)))
#      pos <- pos[-1]-1
#      expect_that(round(as.vector(fitted.values(No6$glm)[pos]),digits=1), equals(c(1.3,3.6,3.6,4.7,4.0,2.6)))
      
      keepNo5 <- rowSums(histpos.t(6)) < 5
      No5 <- openp(duck, dfreq=TRUE, keep=keepNo5)
      expect_that(round(No5$model.fit[,"deviance"],digits=1), equals(56.8))
      expect_that(No5$model.fit[,"    df"], equals(42))
#      expect_that(round(as.vector(No5$N[2:5,"estimate"])), equals(c(393,481,385,491)))

# J'arrive à des résultats très proches mais pas égaux à ceux de l'article pour les lignes en commentaires

    })    

    
test_that("Rcapture reproduces some of Agresti (1994) results in Table 2", {
      data(hare)
      xx <- closedpCI.t(hare, m="Mth", h="Darroch")
      expect_that(round(xx$results[,"abundance"],digits=1), equals(90.5))
      expect_that(round(xx$results[,"deviance"],digits=1), equals(50.7))
      expect_that(xx$results[,"df"], equals(55))
    })    


test_that("Rcapture reproduces some of Agresti (1994) results in Table 2", {
      data(hare)
      M0 <- closedpCI.t(hare)
      expect_that(round(M0$results[,"abundance"]), equals(75))
      expect_that(round(M0$CI[,"InfCL"]), equals(68+1))
#      expect_that(round(M0$CI[,"SupCL"]), equals(68+16))
      
      Mt <- closedpCI.t(hare, m="Mt")
      expect_that(round(Mt$results[,"abundance"]), equals(75))
      expect_that(round(Mt$CI[,"InfCL"]), equals(68+1))
#      expect_that(round(Mt$CI[,"SupCL"]), equals(68+15))

# J'arrive à des résultats très proches mais pas égaux à ceux de l'article pour les lignes en commentaires

# On devrait se poser des questions sur notre calcul d'IC profile. Ça ne marche pas très bien.
    })    

### Récapitulatif par jeu de données du package :

## BBS2001
# Pas dans article JSS
# J'ai essayé autant que possible de reproduire les résultats de Dorazio et Royle (2003)

## bunting
# Article JSS p.18 -> on reproduit des résultats de Cormack 1993 p.46

## catb
# Pas dans article JSS
# J'ai essayé autant que possible de reproduire les résultats de Royle (2006)

## duck
# Article JSS p.21 -> on reproduit des résultats de Cormack 1989 tables 5 et 6

## hare
# Article JSS p.8 -> On reproduit les résultats de Agresti (1994) table 2
#                    + ceux de Cormack 1989 table4
# A comparer avec Coull&Agresti1999 : N=92 et sigma=0.97 : on a les mêmes résultats
# À comparer avec Dorazio&Royle2003 : Table 1, résultats prochent de logistic-normal
# À comparer avec Cormack 1992 pour les IC profile

## HIV
# Pour le jeu de données HIV, on obtient des résultats proches de ceux de Abeni et al. (1994)
# mais pas égaux car pas même méthode d'estimation (mentionné dans notre article JSS p.11)

## mvole
# À venir : article de Louis-Paul et moi 2007

## rbvole
# À venir : article de Louis-Paul 2004

