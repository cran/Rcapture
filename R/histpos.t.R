"histpos.t" <- function(t)
{

    #####################################################################################################################################
    # Validation de l'argument fourni en entrée
    if (any((t %% 1)!=0)) stop("'t' must be an integer")
    #####################################################################################################################################

        X0 <- cbind(c(1, 1, 0), c(1, 0, 1))
        if (isTRUE(all.equal(t,2)))  # il y n'y a que 2 occasions de capture
        { 
            X0 
        } else    # il y a plus de 2 occasions de capture
                {
                  for(i in (3:t))
                        {
                                X1 <- cbind( c(rep(1, 2^(i-1)), rep(0,((2^(i-1))-1))), rbind(X0,rep(0,(i-1)),X0))
                                X0 <- X1
                        }
                X1
                }

}
