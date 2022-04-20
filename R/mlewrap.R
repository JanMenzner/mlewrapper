#' @export
mlewrap <-
  function(y, # n:1-Matrix (DV-Values)
           X, # n:k-matrix (includes Int. & IV-Values
           Z, # Predictor Variables for Variance
           theta, # Parameter-Vector
           par = if(ll=="linear"){
             if(hetero){
               rep(0, (ncol(X) + ncol(Z)))
             } else{
               rep(0, ncol(X) + 1)
             }
           } else{rep(0, ncol(X))}, # Creates vector of 0s as start values
           par_scobit = c(res_logit$par, 0), # Sets scobit start values per default as parameters estimated by logistic regression
           control, # Defines control for optim()
           method, # Defines method of optim()
           lower=-Inf, # Bounds on the variables for the "L-BFGS-B" method
           seed = F, # If T, seed = 1312
           digits, # Defines decimal points in output
           ci=95, # Define CI Level
           bootstrap=F, # If T, bootstrapped SE, p- and t-values included
           samplesize = nrow(X), # Defines Size of Bootstrap Samples
           bs.no = 1000, # Defines No. of Bootstrap Samples
           IV = c(paste("Var", 1:(ncol(X) - 1))), # Character-Vector of IV-Names
           ZV = c(paste("Z", 1:(ncol(Z) - 1))), # Character-Vector of Variable-Names predicting sigma2
           ll="linear", # Specify log likelihood-function (default = linear)
           hetero = F, # Specify whether heteroskedastic regression or not
           parallel = F # If specified, test proportional odds assumption
  ) {

    # Making the function smarter
    if(ll!="ologit" & ll!="oprobit"){
      if (!all(X[, 1] == 1)) {
      X <- cbind(1, X)
      }
    }
    if(hetero){
      if (!all(Z[, 1] == 1)) {
        Z <- cbind(1, Z)
      }}

    ##### Likelihood Function ---
    if (ll=="linear"){
      if(!hetero){
        lik_func <- function(y, X, theta) {# Linear-Regression, Not Heteroskedastic
          n <- nrow(X)                                # No. of observations
          k <- ncol(X)                                # No. of betas
          beta <- theta[1:k]                      # Subset the parameter vector theta
          gamma <- theta[k+1]
          sigma2 <- exp(gamma)                # Ensure that sigma^2 is positive
          e <- y - X %*% beta                     # Calculate the residuals
          llik <- - (n/2)* log(sigma2) -          # Estimate Log-Likelihood
            (t (e) %*% (e) / (2*sigma2) )
          return(llik)                                # Return Result for optim()
        }
      } else{
        lik_func <- function(theta, y, X, Z) {# Linear-Regression, Heteroskedastic
          n <- nrow(X)                                # No. of observations
          k <- ncol(X)                                # No. of betas
          l <- ncol(Z)                                # No. of gammas
          beta <- theta[1:k]                      # Subset the parameter vector theta
          gamma <- theta[(k+1):(k+l)]
          sigma2 <- exp(Z %*% gamma)          # Parameterize sigma squared
          e <- y - X %*% beta                     # Calculate the residuals
          ll <- -1/2*log(sigma2) -              # Estimate Log-Likelihood
            1/2*(e^2/(sigma2))
          ll <- sum(ll)
          return(ll)                                # Return Result for optim()
        }
      }
    } else if(ll=="logit" | ll=="scobit"){# Logit

      if(length(unique(y)) > 2){stop("Dependent Variable not Binary")}

      lik_func <- function(theta, y, X) {
        # theta consists merely of beta (dim is ncol(X))
        beta <- theta[1:ncol(X)]
        # linear predictor; make sure that X is stored as.matrix
        mu <- X %*% beta
        # link function
        p <- 1/(1 + exp(-mu))
        # log-likelihood
        ll <- y * log(p) + (1 - y) * log(1 - p)
        # sum
        ll <- sum(ll)
        return(ll)
      }
    } else if(ll=="ologit" | ll=="oprobit"){# Ordered Logit or probit
      cats <- sort(unique(y))  # Different categories in dep.var
      K <- length(unique(y))  # Number of categories
      L <- matrix(NA, nrow = length(y), ncol = K)  # Empty indicator matrix
      for (j in 1:K) { # Fills indicator matrix
        L[, j] <- y == cats[j]
      }
      lik_func <- function(theta, L, X) {# Log-Likelihood Function
        k <- ncol(X)
        J <- ncol(L)

        beta <- theta[1:k]                                 # Subset theta
        tau <- theta[(k + 1):(k + J - 1)]                  # Subset theta
        U <- X %*% beta                                    # Linear Predictor
        probs <- matrix(nrow = length(U), ncol = J)        # Probabilities in a matrix
        # Calculate Probabilities for the different tau values
        if(ll=="ologit"){
        probs[, 1] <- 1 / (1 + exp(-(tau[1] - U)))         # First Category
        for (j in 2:(J - 1)) {                             # Between categories
          probs[, j] <- 1 / (1 + exp(-(tau[j] - U))) -
            1 / (1 + exp(-(tau[j - 1] - U)))
        }
        probs[, J] <- 1 - 1 / (1 + exp(-(tau[J - 1] - U))) # Final Category
        }
        if(ll=="oprobit"){
        probs[, 1] <- pnorm(tau[1] - U)                    # First Category
        for (j in 2:(J - 1)) {                             # Between categories
          probs[, j] <- pnorm(tau[j] - U) -
            pnorm(tau[j - 1] - U)
        }
        probs[, J] <- 1 - pnorm(tau[J - 1] - U)            # Final Category
        }
        ll <- sum(log(probs[L]))                           # sum over probabilities
        return(ll)
      }
    } else{ # Probit
      lik_func <- function(theta, y, X) {
        # theta consists merely of beta (dim is ncol(X))
        beta <- theta[1:ncol(X)]
        # linear predictor; make sure that X is stored as.matrix
        mu <- X %*% beta
        # link function
        p <- pnorm(mu)
        # log-likelihood
        ll <- y * log(p) + (1 - y) * log(1-p)
        # sum
        ll <- sum(ll)
        return(ll)
      }
    }

    ##### Optimization ---
    if(hetero){
      res <-  stats::optim(
        par = par,                                # Define Startvalues for Parameters
        fn = lik_func,                            # Input function from above
        y = y,                                    # Input DV
        X = X,
        Z = Z,                                    # Input Intercept and IVs
        control = control,                        # Define Control-Option
        method = method,                          # Define Climbing-Option
        hessian = T,                              # Include Hessian Matrix in Return
        lower = lower)
    } else if (ll=="ologit" | ll=="oprobit"){
      if(missing(par)){#  Sets 0s for betas and increasing values for taus
        par = c(rep(0,ncol(X)), seq(-1,0,length.out = K-1))
      }
      suppressWarnings(
        res <- stats::optim(
        par = par,                                # Define Startvalues for Parameters
        fn = lik_func,                            # Input function from above
        L = L,                                    # Input DV
        X = X,                                    # Input Intercept and IVs
        control = control,                        # Define Control-Option
        method = method,                          # Define Climbing-Option
        hessian = T,                              # Include Hessian Matrix in Return
        lower = lower)
      )
    } else {
      res <- stats::optim(
        par = par,                                # Define Startvalues for Parameters
        fn = lik_func,                            # Input function from above
        y = y,                                    # Input DV
        X = X,                                    # Input Intercept and IVs
        control = control,                        # Define Control-Option
        method = method,                          # Define Climbing-Option
        hessian = T,                              # Include Hessian Matrix in Return
        lower = lower)
    }

    if(ll=="scobit"){ # Scobit
      res_logit <- res

      lik_func <- function(theta, y, X) {

        beta <- theta[1:ncol(X)]
        gamma <- theta[ncol(X)+1]

        mu <- X %*% beta

        alpha <- exp(gamma)

        p <- 1/((1+exp(-mu))^alpha)

        ll <- y * log(p) + (1 - y) * log(1 - p)
        ll <- sum(ll)
        return(ll)
      }
      res <- stats::optim(
        par = par_scobit,                         # Define Startvalues for Parameters (Defualt: logit parameters)
        fn = lik_func,                            # Input function from above
        y = y,                                    # Input DV
        X = X,                                    # Input Intercept and IVs
        control = control,                        # Define Control-Option
        method = method,                          # Define Climbing-Option
        hessian = T,                              # Include Hessian Matrix in Return
        lower = lower)
    }

    ##### Adding information ---
    if(ll=="linear" & !hetero | ll=="scobit"){
      res$par[length(res$par)] <-                 # Transform gamma to sigma2/alpha
        exp(res$par[length(res$par)])}
    res$se <- sqrt(diag(solve(-(res$hessian))))   # Standard errors
    t <- res$par / res$se                         # t-values
    n <- nrow(X)                                  # No. of Observations
    df <- n - length(res$par)                     # Degrees of freedom
    p <- 2 * pt(-abs(t), df)                      # p-values
    low_ci <- res$par - qnorm((100-ci)/2/100,     # lower bound ci
                              lower.tail=FALSE)*res$se
    up_ci <- res$par + qnorm((100-ci)/2/100,      # upper bound ci
                             lower.tail=FALSE)*res$se

    res$table <- cbind(Est. = res$par,            # Combine results in table format
                       St.Err. =res$se,
                       lowci = low_ci,
                       upci = up_ci,
                       T.Val = t,
                       P.Val = p)

    res$betas <- res$table[1:ncol(X),] # Take only Beta-Coefficients
    if(ll=="ologit" | ll=="oprobit"){
      if(missing(IV)){IV = c(paste("Var", 1:(ncol(X))))}
      rownames(res$table) <- c(IV, paste("Tau", 1:(ncol(L)-1)))
      rownames(res$betas) <- c(IV)
      }
    else{
    rownames(res$betas) <-  c("Int.", IV)         # Give Table Row-Names
    }

    colnames(res$betas)[3] <- colnames(res$table)[3] <-  # Label CIs correctly
      paste("low", paste0(ci,"%"), "CI")
    colnames(res$betas)[4] <- colnames(res$table)[4] <-
      paste("up", paste0(ci,"%"), "CI")

    if(ll=="logit" | ll=="probit"){
      res$aic = -2 * res$value + 2 * (ncol(X)+1)
      res$bic = -2 * res$value + log(n) * (ncol(X)+1)
    }

    if(ll=="ologit" | ll=="oprobit"){
      res$aic = -2 * res$value + 2 * (ncol(X)+length(unique(y))-1)
      res$bic = -2 * res$value + log(n) * (ncol(X)+length(unique(y))-1)
    }

    if(hetero){
      res$gammas <- res$table[(ncol(X)+1):(ncol(X)+ncol(Z)),] # Take only Gamma-Coefficients
      rownames(res$gammas) <-  c("Int.", ZV)         # Give Table Row-Names
      colnames(res$gammas)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly
      res$meanvar <- apply(exp(Z %*% res$par[(ncol(X)+1):(ncol(X)+ncol(Z))]),2,mean)
    }

    if(ll=="scobit"){
      # Compute log-likelihood ratio between scobit and logit models
      R <- 2 * (res$value - res_logit$value)
      llratio <- pchisq(R,
                        df = 1,
                        lower.tail = FALSE)
      # --> Significant Result indicates that the scobit model is better
    }

    ##### Optional: Add bootstrapped SE, p- and t-values ---
    if(bootstrap){

      if (seed==T){set.seed(1312)}                # Set Seed for replicability

      if(hetero){
        mle_bs <- function(){
          data <- cbind(y,X,Z)                      # Prepare Data Set
          bs_samp  <-                               # Draw Bootstrap Sample
            data[sample(1:samplesize, replace = TRUE),]
          y_bs <- bs_samp[,1]                       # Index Bootstrap DV
          X_bs <- bs_samp[,2:(ncol(X)+1)]           # Index Int. and Bootstrap IVs
          Z_bs <- bs_samp[,-c(1:(ncol(X)+1))]       # Index Sigma2 IVs

          res_bs  <- stats::optim(                  # Running similar optimization as above:
            par = par,                              # Define Startvalues for Parameters
            fn = lik_func,                          # Input function from above
            y = y_bs,                               # Input Bootstrap DV
            X = X_bs,                               # Input Bootstrap Intercept and IVs
            Z = Z_bs,                                # Input Bootstrap Sigma2 IVs
            control = control,                      # Define Control-Option
            method = method,                         # Define Climbing-Option
            lower = lower)
          return(res_bs$par)}
      } else if(ll=="ologit" | ll=="oprobit"){
        mle_bs <- function(){
          data <- cbind(y,X  )                      # Prepare Data Set
          bs_samp  <-                               # Draw Bootstrap Sample
            data[sample(1:samplesize, replace = TRUE),]
          y_bs <- bs_samp[,1]                       # Index Bootstrap DV
          X_bs <- bs_samp[,2:(ncol(X)+1)]           # Index Bootstrap IVs
          cats_bs <- sort(unique(y_bs))  # Different categories in dep.var
          K_bs <- length(unique(y_bs))  # Number of categories
          L_bs <- matrix(NA, nrow = length(y_bs), ncol = K_bs)  # Empty indicator matrix
          for (j in 1:K_bs) { # Fills indicator matrix
            L_bs[, j] <- y_bs == cats[j]
          }

          res_bs  <- stats::optim(                  # Running similar optimization as above:
            par = par,                              # Define Startvalues for Parameters
            fn = lik_func,                          # Input function from above
            L = L_bs,                               # Input Bootstrap DV
            X = X_bs,                               # Input Bootstrap Intercept and IVs
            control = control,                      # Define Control-Option
            method = method,                         # Define Climbing-Option
            lower = lower)
          return(res_bs$par)}
      } else {
        mle_bs <- function(){
          data <- cbind(y, X)                       # Prepare Data Set
          bs_samp  <-                               # Draw Bootstrap Sample
            data[sample(1:samplesize, replace = TRUE),]
          y_bs <- bs_samp[,1]                       # Index Bootstrap DV
          X_bs <- bs_samp[,2:ncol(bs_samp)]         # Index Int. and Bootstrap IVs

          suppressWarnings(
            res_bs  <- stats::optim(                  # Running similar optimization as above:
            if(ll=="scobit"){par = par_scobit}
            else{par = par},                        # Define Startvalues for Parameters
            fn = lik_func,                          # Input function from above
            y = y_bs,                               # Input Bootstrap DV
            X = X_bs,                               # Input Bootstrap Intercept and IVs
            control = control,                      # Define Control-Option
            method = method,                        # Define Climbing-Option
            lower = lower) )
          return(res_bs$par)
          }

      }

      suppressWarnings(
        bs_res <- t(replicate(bs.no, mle_bs()))     # Repeat this procedure bs.no times
      )
      res$bs_se <- sqrt(apply(bs_res, 2, var))    # Extract the Bootstrap standard errors
      t_bs <- res$par / res$bs_se                 # Bootstrapped t-values
      p_bs <- 2 * pt(-abs(t_bs), df)              # Bootstrapped p-values
      low_ci_bs <- res$par - qnorm((100-ci)/2/100, # lower bound ci
                                   lower.tail=FALSE)*res$bs_se
      up_ci_bs <- res$par + qnorm((100-ci)/2/100,  # upper bound ci
                                  lower.tail=FALSE)*res$bs_se

      res$table_bs <- cbind(Est. = res$par,       # Combine results in table format
                            Boot.SE =res$bs_se,
                            lowci = low_ci_bs,
                            upci = up_ci_bs,
                            Boot.T.Val = t_bs,
                            Boot.P.Val = p_bs)

      if(ll=="ologit" | ll=="oprobit"){
        res$betas_bs <- res$table_bs
        rownames(res$table_bs) <- rownames(res$table)
      } else {res$betas_bs <- res$table_bs[1:ncol(X),] # Exclude Sigma^2-Row
      rownames(res$betas_bs) <- rownames(res$betas) # Give Table Row-Names
      }
      colnames(res$betas_bs)[3:4] <- colnames(res$table_bs)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly

      if(hetero){
        res$gammas_bs <- res$table_bs[(ncol(X)+1):(ncol(X)+ncol(Z)),] # Take only Gamma-Coefficients
        rownames(res$gammas_bs) <- c("Int.", ZV)         # Give Table Row-Names
        colnames(res$gammas_bs)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly
        res$meanvar <- apply(exp(Z %*% res$par[(ncol(X)+1):(ncol(X)+ncol(Z))]),2,mean)
      }
    }



    # Add Test for proportional odds assumption (logit)
    if(parallel){
      M <- sort(unique(y))                   # Creates Ordered Vector with all values in y
      M <- M[c(-1)]                          # Excludes first value
      Y <- matrix(NA, nrow = length(y), ncol = length(M))
      X_test <- cbind(1, X)                  # Add intercept
      for(i in M){
        Y[,i] <- ifelse(as.numeric(y) < i, 0, 1)
      }

      test_lik <- function(theta, y, X) {
        # theta consists merely of beta (dim is ncol(X))
        beta <- theta[1:ncol(X)]
        # linear predictor; make sure that X is stored as.matrix
        mu <- X %*% beta
        # link function
        if(ll=="ologit"){p <- 1/(1 + exp(-mu))}
        if(ll=="oprobit"){p <- pnorm(mu)}
        # log-likelihood
        ll <- y * log(p) + (1 - y) * log(1 - p)
        # sum
        ll <- sum(ll)
        return(ll)
      }

      prop <- list()
      res$prop_beta <- matrix(NA, nrow = (1+ncol(X)+length(M)), ncol = (1 + length(M)))
      res$prop_se <- matrix(NA, nrow = (1+ncol(X)+length(M)), ncol = (1 + length(M)))
      if(missing(IV)){IV = c(paste("Var", 1:(ncol(X))))}
      rownames(res$prop_se) <- rownames(res$prop_beta) <- c("Int.", IV, paste("Tau", 1:(length(M))))
      res$prop_beta[,1] <- c(NA, res$par)
      res$prop_se[,1] <- c(NA, res$se)
      for(i in 1:ncol(Y)){
        prop[[i]] <- stats::optim(
          par = c(0, par[1:ncol(X)]),               # Define Startvalues for Parameters
          fn = test_lik,                            # Input function from above
          y = Y[,i],                                # Input DV
          X = X_test,                               # Input Intercept and IVs
          control = control,                        # Define Control-Option
          method = method,                          # Define Climbing-Option
          hessian = T,
          lower = lower)

        prop[[i]]$se <- sqrt(diag(solve(-(prop[[i]]$hessian))))   # Standard error

        res$prop_beta[,1+i] <- c(prop[[i]]$par, rep(NA, length(M)))
        res$prop_se[,1+i] <- c(prop[[i]]$se, rep(NA, length(M)))
      }

      colnames(res$prop_beta) <- colnames(res$prop_se) <-
        c("Ordered", paste(paste0("<", min(M):max(M)), "= 0"))
    }


    ##### Format Output ---
    if (missing(digits)) {                        # Prints Coefficient-Table,
      cat("Beta-Coefficients \n")                 # Log.Likelihood and Variance &
      if(ll=="ologit" | ll=="oprobit"){print(res$table)} else {print(res$betas)} # saves output
      cat("\n")
      if(hetero){
        cat("Gamma-Coefficients \n")
        print(res$gammas)
        cat("\n")}
      if(bootstrap){
        cat("Bootstrapping \n")
        if(ll=="ologit" | ll=="oprobit"){print(res$table_bs)} else {print(res$betas_bs)}
        cat("\n")
        if(hetero){
          print(res$gammas_bs)
          cat("\n")}}
      cat("Log. Likelihood:", res$value, "\n")
      if(ll=="linear"){
        if(hetero){
          cat("Mean Var.", res$meanvar)
        } else{ cat("Variance:", res$par[length(res$par)])}
      } else if(ll=="scobit"){
        cat("Alpha:", res$par[length(res$par)])
        cat("\n")
        cat("Likelihood-ratio test of alpha=1: Prob > chi2 =", llratio)
        if(llratio <= 0.05){
          cat("\n")
          cat(" --> Scobit-Model should be prefered")
        } else{
          cat("\n")
          cat(" --> Logit-Model should be prefered")
        }
      } else {
        cat("AIC:", res$aic)
        cat("\n")
        cat("BIC:", res$bic)
      }
      if(parallel == T){
        cat("\n")
        cat("\n", "Coefficients Proportional?", "\n")
        print(res$prop_beta, na.print = "")
        cat("Standard Errors Proportional?", "\n")
        print(res$prop_se, na.print = "")
      }
      invisible(res)
    }


    else{
      cat("Beta-Coefficients \n")                 # Prints Coefficient-Table,
      if(ll=="ologit" | ll=="oprobit"){print(apply(res$table, c(1,2),round,digits)) # Log.Likelihood and Variance
        } else {print(apply(res$betas, c(1,2),round,digits))} # rounded to digits-argument &
      cat("\n")                                               # saves output
      if(hetero){
        cat("Gamma-Coefficients \n")
        res$gammas <- apply(res$gammas, c(1,2),round,digits)
        print(res$gammas)
        cat("\n")}
      if(bootstrap){
        cat("Bootstrapping \n")
        if(ll=="ologit" | ll=="oprobit"){print(apply(res$table_bs, c(1,2),round,digits))
        } else {print(apply(res$betas_bs, c(1,2),round,digits))}
        cat("\n")
        if(hetero){
          print(apply(res$gammas_bs, c(1,2),round,digits))
          cat("\n")}}
      cat("Log. Likelihood:", round(res$value,digits), "\n")
      if(ll=="linear"){
        if(hetero){
          cat("Mean Var.", round(res$meanvar, digits))
        } else{cat("Variance:", round(res$par[length(res$par)], digits))}
      } else if(ll=="scobit"){
        cat("Alpha:", round(res$par[length(res$par)], digits))
        cat("\n")
        cat("Likelihood-ratio test of alpha=1: Prob > chi2 =", round(llratio, digits))
        if(llratio <= 0.05){
          cat("\n")
          cat(" --> Scobit-Model should be prefered")
        } else{
          cat("\n")
          cat(" --> Logit-Model should be prefered")
        }
      } else {
        cat("AIC:", round(res$aic, digits))
        cat("\n")
        cat("BIC:", round(res$bic, digits))
      }
      if(parallel == T){
        cat("\n")
        cat("\n", "Coefficients Proportional?", "\n")
        print(apply(res$prop_beta, c(1,2),round,digits), na.print = "--")
        cat("Standard Errors Proportional?", "\n")
        print(apply(res$prop_se, c(1,2),round,digits), na.print = "--")
      }
      invisible(res)
    }
  }


