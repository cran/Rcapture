context("Arguments validation and default values")

test_that("meaningful errors are printed if 'X' is not adequate", {
      expect_that(closedp(c(1,2)), throws_error("'X' must contain at least 2 capture occasions"))
      expect_that(closedp(closedp(matrix(c(1,1,2,0),2,2))), 
          throws_error("'X' must contain only zeros and ones"))
      
      data(hare)      
      for (fct in c("closedpCI.t","closedpCI.0","closedp.bc"))
        expect_that(do.call(fct,args=list(X=hare[,1:2],m="Mh")),
            throws_error("heterogeneous models require at least 3 capture occasions"))      
            
      data(BBS2001)
      expect_that(closedpCI.0(BBS2001, dfreq=TRUE, dtype="nbcap", t=20),
          throws_error("the first column of 'X' must contain only integers between 1 and"))
      expect_that(closedpCI.0(BBS2001, dtype="nbcap", t=20),
          throws_error("'X' must have 1 column"))
      expect_that(closedpCI.0(BBS2001[,1], dfreq=TRUE, dtype="nbcap", t=20),
          throws_error("'X' must have 2 columns"))
    })

test_that("meaningful errors are printed if 't' or 't0' is not adequate", {
      data(hare)
      expect_that(closedp.bc(X=hare, t=Inf),
          throws_error("'t' can take the value 'Inf' only with the functions 'closedp.0', 'closedpCI.0' and 'descriptive'"))
      expect_that(closedpCI.0(hare[,1:3],t0=2,m="Mh",h="Normal"),
          throws_error("the given t0 argument is too small for the requested model"))
      for (t0 in c(1,7,3.5))
        expect_that(closedpCI.0(hare,t0=t0,m="Mh",h="Normal"),
            throws_error("'t0' must be an integer between 2 and"))
      
      expect_that(closedpCI.0(hare, t=5),
          throws_error("'t' is not equal to the number of columns in 'X'"))
      for (t in c(1,3.5))
        expect_that(closedpCI.0(hare, t=t),
            throws_error("if not Null, 't' must be an integer greater or equal to 2 or take the value 'Inf'"))
      expect_that(closedpCI.t(hare, t=5),
          throws_error("unused argument"))

      for (m in c("Mt","Mth","Mb","Mbh"))
          expect_that(closedp.bc(hare, t0=5, m=m),
              gives_warning("the input argument 't0' could not be used with the requested model"))
      
      data(HIV)
      expect_that(closedpCI.0(HIV, dfreq=TRUE, t=5),
          throws_error("'t' is not equal to the number of columns in 'X' minus 1"))
      
      data(BBS2001)
      expect_that(closedpCI.0(BBS2001, dfreq=TRUE, dtype="nbcap"),
          throws_error("argument 't' must be given if 'dtype' takes the value \"nbcap\""))
    })

test_that("meaningful errors are printed if 'mX' is not adequate", {
      data(hare)
      for (fct in c("closedpCI.t","closedpCI.0"))
        expect_that(do.call(fct,args=list(X=hare,mX=function(x)x^2)),
            throws_error("cannot turn the given 'mX' argument into a matrix"))
      expect_that(closedpCI.0(hare,mX=histpos.t(6)),
          throws_error("'mX' must have t0 rows"))
      expect_that(closedpCI.0(hare,mX=~.),
          throws_error("'mX' cannot be a formula for 'closedpCI.0'"))
      expect_that(closedpCI.t(hare,mX=~c1+c10),
          throws_error("the only accepted variable names in the formula"))
      expect_that(closedpCI.t(hare,mX=~0+.),
          throws_error("the intercept must not be removed"))
      expect_that(closedpCI.t(hare,mX=~.-1),
          throws_error("the intercept must not be removed"))
    })

test_that("meaningful errors are printed if 'm', 'h' or 'theta' (sometimes in 'h.control') is not adequate", {
      data(hare)
      psi <- function(x) { 0.5^x - 1 }
      
      for (fct in c("closedpCI.t","closedpCI.0"))
        expect_that(do.call(fct,args=list(X=hare, m="Mh", h="test")),
            throws_error("'h' must be a function or a character string taking one of these values: \"Chao\", \"LB\", \"Poisson\", \"Darroch\", \"Gamma\", \"Normal\""))
      for (h in c("Poisson", "Darroch")) {
        for (fct in c("closedpCI.t","closedpCI.0"))
          expect_that(do.call(fct,args=list(X=hare, m="Mh", h=h, h.control=list(theta=0))),
              throws_error("'theta' must be a positive number"))      
        expect_that(closedp.bc(X=hare, m="Mh", h=h, theta=0),
            throws_error("'theta' must be a positive number"))      
      }
      
      # Vérifications avec closedpCI.t
      expect_that(closedpCI.t(hare,m="Mb"), 
          throws_error("'m' can only take one of these values: \"M0\", \"Mt\", \"Mh\", \"Mth\""))

      # Vérifications avec closedpCI.0
      expect_that(closedpCI.0(hare,m="Mt"), 
          throws_error("models with a temporal effect cannot be adjusted"))
      expect_that(closedpCI.0(hare,m="Mb"), 
          throws_error("'m' can only take one of these values: \"M0\", \"Mh\""))
      
      # Vérifications avec closedp.bc
      expect_that(closedp.bc(hare,m="Mth",h="Poisson",theta=3), 
          throws_error("the biais correction can be performed for model Mth Poisson only with theta equals to 2"))
      expect_that(closedp.bc(hare,m="Mhh"), 
          throws_error("'m' can only take one of these values: \"M0\", \"Mt\", \"Mh\", \"Mth\", \"Mb\", \"Mbh\""))
      expect_that(closedp.bc(hare,m="Mth",h="Gamma"), 
          throws_error("the biais correction cannot be performed for model Mth when 'h' is \"Gamma\" or a function"))
      expect_that(closedp.bc(hare,m="Mth",h=psi), 
          throws_error("the biais correction cannot be performed for model Mth when 'h' is \"Gamma\" or a function"))
      expect_that(closedp.bc(hare[,1:3],m="Mbh"), 
          throws_error("the biais correction cannot be performed for model Mbh with less than 4 capture occasions"))
      expect_that(closedp.bc(X=hare, m="Mh", h="Normal"),
          throws_error("'h' must be a function or a character string taking one of these values: \"Chao\", \"LB\", \"Poisson\", \"Darroch\", \"Gamma\""))
      
    })

test_that("'m', 'h', 'theta' (sometimes in 'h.control') and 'mname' default values are correct", {
      data(hare)
      psi <- function(x) { 0.5^x - 1 }
      
      # Vérifications avec closedpCI.t
      expect_that(rownames(closedpCI.t(hare)$results), equals("M0"))
      expect_that(rownames(closedpCI.t(hare,m="M0")$results), equals("M0"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh")$results), equals("Mh Chao (LB)"))
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Chao")$results), equals("Mh Chao (LB)"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="LB")$results), equals("Mh LB"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Darroch")$results), equals("Mh Darroch"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Poisson")$results), equals("Mh Poisson2"))
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Poisson",h.control=list(theta=3))$results), equals("Mh Poisson3"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Gamma")$results), equals("Mh Gamma3.5"))
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Gamma",h.control=list(theta=2))$results), equals("Mh Gamma2"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h="Normal")$results), equals("Mh Normal"))      
      expect_that(rownames(closedpCI.t(hare,m="Mh",h=psi)$results), equals("Mh psi"))
      
      expect_that(rownames(closedpCI.t(hare,m="Mt")$results), equals("Mt"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth")$results), equals("Mth Chao (LB)"))
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Chao")$results), equals("Mth Chao (LB)"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="LB")$results), equals("Mth LB"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Darroch")$results), equals("Mth Darroch"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Poisson")$results), equals("Mth Poisson2"))
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Poisson",h.control=list(theta=3))$results), equals("Mth Poisson3"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Gamma")$results), equals("Mth Gamma3.5"))
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Gamma",h.control=list(theta=2))$results), equals("Mth Gamma2"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h="Normal")$results), equals("Mth Normal"))      
      expect_that(rownames(closedpCI.t(hare,m="Mth",h=psi)$results), equals("Mth psi"))
              
      matX <- histpos.t(6)
      expect_that(rownames(closedpCI.t(hare,mX=matX)$results), equals("matX"))
      expect_that(rownames(closedpCI.t(hare,mX=histpos.t(6))$results), equals("mX"))
      expect_that(rownames(closedpCI.t(hare,mX=~.)$results), equals("mX"))
      expect_that(rownames(closedpCI.t(hare,mX=matX,h="Poisson")$results), equals("matX+h=Poisson2"))
      expect_that(rownames(closedpCI.t(hare,mX=histpos.t(6),h="Gamma")$results), equals("mX+h=Gamma3.5"))
      expect_that(rownames(closedpCI.t(hare,mX=~.,h="Normal")$results), equals("mX+h=Normal"))
      
      
      # Vérifications avec closedpCI.0
      expect_that(rownames(closedpCI.0(hare)$results), equals("M0"))
      expect_that(rownames(closedpCI.0(hare,m="M0")$results), equals("M0"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh")$results), equals("Mh Chao (LB)"))
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Chao")$results), equals("Mh Chao (LB)"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="LB")$results), equals("Mh LB"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Darroch")$results), equals("Mh Darroch"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Poisson")$results), equals("Mh Poisson2"))
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Poisson",h.control=list(theta=3))$results), equals("Mh Poisson3"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Gamma")$results), equals("Mh Gamma3.5"))
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Gamma",h.control=list(theta=2))$results), equals("Mh Gamma2"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h="Normal")$results), equals("Mh Normal"))      
      expect_that(rownames(closedpCI.0(hare,m="Mh",h=psi)$results), equals("Mh psi"))

      matX <- histpos.0(6)
      expect_that(rownames(closedpCI.0(hare,mX=matX)$results), equals("matX"))
      expect_that(rownames(closedpCI.0(hare,mX=histpos.0(6))$results), equals("mX"))
      expect_that(rownames(closedpCI.0(hare,mX=matX,h="Chao")$results), equals("matX+h=Chao(LB)"))
      expect_that(rownames(closedpCI.0(hare,mX=histpos.0(6),h="LB")$results), equals("mX+h=LB"))

      # Vérifications avec closedp.bc
      expect_that(rownames(closedp.bc(X=hare)$results), equals("M0"))
      expect_that(rownames(closedp.bc(X=hare, m="M0")$results), equals("M0"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh")$results), equals("Mh Chao (LB)"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Chao")$results), equals("Mh Chao (LB)"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="LB")$results), equals("Mh LB"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Darroch")$results), equals("Mh Darroch"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Poisson")$results), equals("Mh Poisson2"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Poisson", theta=3)$results), equals("Mh Poisson3"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Gamma")$results), equals("Mh Gamma3.5"))
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h="Gamma", theta=2)$results), equals("Mh Gamma2"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mh", h=psi)$results), equals("Mh psi"))
      
      expect_that(rownames(closedp.bc(X=hare, m="Mt")$results), equals("Mt"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mth")$results), equals("Mth Chao (LB)"))
      expect_that(rownames(closedp.bc(X=hare, m="Mth", h="Chao")$results), equals("Mth Chao (LB)"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mth", h="LB")$results), equals("Mth LB"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mth", h="Darroch")$results), equals("Mth Darroch"))      
      expect_that(rownames(closedp.bc(X=hare, m="Mth", h="Poisson")$results), equals("Mth Poisson2"))
      
    })


### Je pourrais passer en revue tous les appels à la fonction stop
### et faire un test pour chaque erreur.