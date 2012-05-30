context("Mathematical validation")
      
test_that("the .t and .0 'closedp' and 'closedpCI' functions give the same results for the same models", {
      data(hare)
      
      res.t <- closedp.t(X=hare, dfreq=FALSE)
      res.0 <- closedp.0(X=hare, dfreq=FALSE)
      
      fct.t <- closedpCI.t(X=hare,dfreq=FALSE,m="Mh",h="Normal")
      fct.0 <- closedpCI.0(X=hare,dfreq=FALSE,dtype="hist",m="Mh",h="Normal")
      
      psi <- function(x) { 0.5^x - 1 }
      matX.t <- rowSums(histpos.t(6))
      Mh.t <- closedpCI.t(X=hare,dfreq=FALSE,mX=matX.t,h=psi)
      matX.0 <- histpos.0(6)
      Mh.0 <- closedpCI.0(X=hare,dfreq=FALSE,mX=matX.0,h=psi)
      
      expect_that(res.0$results["M0","abundance"], equals(res.t$results["M0","abundance"]))
      expect_that(res.0$results["M0","stderr"], 
                  equals(res.t$results["M0","stderr"], tolerance=0.0001))
      expect_that(res.0$results["M0","df"], equals(res.t$results["M0","df"]-(2^6-1-6)))
              
      expect_that(res.0$results["Mh Chao (LB)","abundance"], 
                  equals(res.t$results["Mh Chao (LB)","abundance"]))
      expect_that(res.0$results["Mh Chao (LB)","stderr"], 
                  equals(res.t$results["Mh Chao (LB)","stderr"], tolerance=0.0001))
      expect_that(res.0$results["Mh Chao (LB)","df"], 
                  equals(res.t$results["Mh Chao (LB)","df"]-(2^6-1-6)))
              
      expect_that(res.0$results["Mh Poisson2","abundance"], 
                  equals(res.t$results["Mh Poisson2","abundance"]))
      expect_that(res.0$results["Mh Poisson2","stderr"], 
                  equals(res.t$results["Mh Poisson2","stderr"], tolerance=0.0001))
      expect_that(res.0$results["Mh Poisson2","df"], 
                  equals(res.t$results["Mh Poisson2","df"]-(2^6-1-6)))
              
      expect_that(res.0$results["Mh Darroch","abundance"], 
                  equals(res.t$results["Mh Darroch","abundance"]))
      expect_that(res.0$results["Mh Darroch","stderr"], 
                  equals(res.t$results["Mh Darroch","stderr"], tolerance=0.0001))
      expect_that(res.0$results["Mh Darroch","df"], 
                  equals(res.t$results["Mh Darroch","df"]-(2^6-1-6)))
              
      expect_that(res.0$results["Mh Gamma3.5","abundance"], 
                  equals(res.t$results["Mh Gamma3.5","abundance"]))
      expect_that(res.0$results["Mh Gamma3.5","stderr"], 
                  equals(res.t$results["Mh Gamma3.5","stderr"], tolerance=0.0001))
      expect_that(res.0$results["Mh Gamma3.5","df"], 
                  equals(res.t$results["Mh Gamma3.5","df"]-(2^6-1-6)))
              
      expect_that(fct.0$results[,"abundance"], equals(fct.t$results[,"abundance"]))
      expect_that(fct.0$results[,"stderr"], equals(fct.t$results[,"stderr"]))
      expect_that(fct.0$results[,"InfCL"], equals(fct.t$results[,"InfCL"]))
      expect_that(fct.0$results[,"SupCL"], equals(fct.t$results[,"SupCL"]))
      expect_that(fct.0$results[,"df"], equals(fct.t$results[,"df"]-(2^6-1-6)))
      
      expect_that(Mh.0$results[,"abundance"], equals(Mh.t$results[,"abundance"]))
      expect_that(Mh.0$results[,"stderr"], equals(Mh.t$results[,"stderr"], tolerance=0.0001))
      expect_that(Mh.0$results[,"df"], equals(Mh.t$results[,"df"]-(2^6-1-6)))
      expect_that(Mh.0$CI[,"abundance"], equals(Mh.t$CI[,"abundance"]))
      expect_that(Mh.0$CI[,"InfCL"], equals(Mh.t$CI[,"InfCL"]))
      expect_that(Mh.0$CI[,"SupCL"], equals(Mh.t$CI[,"SupCL"]))
    })

test_that("'closedpCI.t' and 'closedp.t' give the same results for the same models", {
      data(hare)
          
      res <- closedp.t(X=hare)
      resCI <- vector(mode="list")
      resCI[[1]] <- closedpCI.t(X=hare,dfreq=FALSE,m="M0")
      resCI[[2]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mt")
      resCI[[3]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mh",h="Chao")
      resCI[[4]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mh",h="Poisson")
      resCI[[5]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mh",h="Darroch")
      resCI[[6]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mh",h="Gamma")
      resCI[[7]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mth",h="Chao")
      resCI[[8]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mth",h="Poisson")
      resCI[[9]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mth",h="Darroch")
      resCI[[10]] <- closedpCI.t(X=hare,dfreq=FALSE,m="Mth",h="Gamma")
      
      for (i in 1:10)
        expect_that(res$results[i,,drop=FALSE], is_identical_to(resCI[[i]]$results))

    })      
    
test_that("the degrees of freedom are good", {
      data(BBS2001)
      m1 <- closedpCI.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=50,m="Mh",h="Normal")
      m2 <- closedpCI.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=50,t0=20,m="Mh",h="Normal")
      m3 <- closedpCI.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=Inf,m="Mh",h="Normal") 
      m4 <- closedpCI.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=Inf,t0=20,m="Mh",h="Normal")
      tobs <- max(BBS2001[BBS2001[,2]!=0, 1])
      
      expect_that(m1$results[,"df"], equals(50-3))
      expect_that(m2$results[,"df"], equals(20-3))
      expect_that(m3$results[,"df"], equals(tobs-3))
      expect_that(m4$results[,"df"], equals(20-3))
    })

test_that("the mX + h arguments works correctly", {
      histpos <- histpos.t(3)
      DarR3 <- cbind(histpos, c(72, 155, 7, 71, 13, 53, 43))

      # Example avec h="Darroch"
      matX <- cbind(histpos,histpos[,1]*histpos[,2],(rowSums(histpos)^2)/2)
      rmX <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=matX,mname="Darroch")
      matX <- cbind(histpos,histpos[,1]*histpos[,2])
      rmXh <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=matX,h="Darroch",mname="Darroch") 
      
      expect_that(rmX$results[,"abundance"], equals(rmXh$results[,"abundance"]))
      expect_that(rmX$results[,"stderr"], equals(rmXh$results[,"stderr"]))
      expect_that(rmX$results[,"deviance"], equals(rmXh$results[,"deviance"]))
      expect_that(rmX$results[,"df"], equals(rmXh$results[,"df"]))

      # Example avec h="Chao", mais sans eta négatif fixés à zéro
      matX <- cbind(histpos,histpos[,1]*histpos[,2],c(1,rep(0,6)))
      rmX <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=matX,mname="LB")
      matX <- cbind(histpos,histpos[,1]*histpos[,2])
      rmXh <- closedpCI.t(X=DarR3,dfreq=TRUE,mX=matX,h="Chao",mname="LB")
      
      expect_that(rmX$results[,"abundance"], equals(rmXh$results[,"abundance"]))
      expect_that(rmX$results[,"stderr"], equals(rmXh$results[,"stderr"]))
      expect_that(rmX$results[,"deviance"], equals(rmXh$results[,"deviance"]))
      expect_that(rmX$results[,"df"], equals(rmXh$results[,"df"]))

      # Example avec h="Chao", avec eta négatif fixés à zéro
      histpos <- histpos.t(4)
      diabetes<-cbind(histpos,c(58,157,18,104,46,650,12,709,14,20,7,74,8,182,10))
      matX <- cbind(histpos,histpos[,1]*histpos[,3],histpos[,2]*histpos[,4],histpos[,3]*histpos[,4])
      nbcap <- rowSums(histpos)
      matX_LB <- cbind(matX, pmax(nbcap-2,0)) # pmax(nbcap-3,0) enlevé car eta négatif
      rmX <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=matX_LB,mname="LB")
      matX_LB <- cbind(matX)
      rmXh <- closedpCI.t(X=diabetes,dfreq=TRUE,mX=matX_LB,h="Chao",mname="LB")
      
      expect_that(rmX$results[,"abundance"], equals(rmXh$results[,"abundance"]))
      expect_that(rmX$results[,"stderr"], equals(rmXh$results[,"stderr"]))
      expect_that(rmX$results[,"deviance"], equals(rmXh$results[,"deviance"]))
      expect_that(rmX$results[,"df"], equals(rmXh$results[,"df"]))
    })

