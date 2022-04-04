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
           hetero = F # Specify whether heteroskedastic regression or not
  ) {

    # Making the function smarter
    if (!all(X[, 1] == 1)) {
      X <- cbind(1, X)
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
          beta_hat <- theta[1:k]                      # Subset the parameter vector theta
          gamma_hat <- theta[k+1]
          sigma2_hat <- exp(gamma_hat)                # Ensure that sigma^2 is positive
          e <- y - X %*% beta_hat                     # Calculate the residuals
          llik <- - (n/2)* log(sigma2_hat) -          # Estimate Log-Likelihood
            (t (e) %*% (e) / (2*sigma2_hat) )
          return(llik)                                # Return Result for optim()
        }
      } else{
        lik_func <- function(theta, y, X, Z) {# Linear-Regression, Heteroskedastic
          n <- nrow(X)                                # No. of observations
          k <- ncol(X)                                # No. of betas
          l <- ncol(Z)                                # No. of gammas
          beta_hat <- theta[1:k]                      # Subset the parameter vector theta
          gamma_hat <- theta[(k+1):(k+l)]
          sigma2_hat <- exp(Z %*% gamma_hat)          # Parameterize sigma squared
          e <- y - X %*% beta_hat                     # Calculate the residuals
          llik <- -1/2*log(sigma2_hat) -              # Estimate Log-Likelihood
            1/2*(e^2/(sigma2_hat))
          llik <- sum(llik)
          return(llik)                                # Return Result for optim()
        }
      }
    } else if(ll=="logit"){# Logit
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
    ##### Adding information ---
    if(ll=="linear" & !hetero){
      res$par[length(res$par)] <-                 # Transform gamma_hat to sigma^2_hat
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
    res$betas <- res$table[1:ncol(X),]            # Take only Beta-Coefficients
    rownames(res$betas) <-  c("Int.", IV)         # Give Table Row-Names
    colnames(res$betas)[3] <-                     # Label CIs correctly
      paste("low", paste0(ci,"%"), "CI")
    colnames(res$betas)[4] <-
      paste("up", paste0(ci,"%"), "CI")

    if(ll=="logit" | ll=="probit"){
      res$aic = -2 * res$value + 2 * (ncol(X)+1)
      res$bic = -2 * res$value + log(n) * (ncol(X)+1)
    }

    if(hetero){
      res$gammas <- res$table[(ncol(X)+1):(ncol(X)+ncol(Z)),] # Take only Gamma-Coefficients
      rownames(res$gammas) <-  c("Int.", ZV)         # Give Table Row-Names
      colnames(res$gammas)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly
      res$meanvar <- apply(exp(Z %*% res$par[(ncol(X)+1):(ncol(X)+ncol(Z))]),2,mean)
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
      } else {
        mle_bs <- function(){
          data <- cbind(y, X)                       # Prepare Data Set
          bs_samp  <-                               # Draw Bootstrap Sample
            data[sample(1:samplesize, replace = TRUE),]
          y_bs <- bs_samp[,1]                       # Index Bootstrap DV
          X_bs <- bs_samp[,2:ncol(bs_samp)]         # Index Int. and Bootstrap IVs

          res_bs  <- stats::optim(                         # Running similar optimization as above:
            par = par,                              # Define Startvalues for Parameters
            fn = lik_func,                          # Input function from above
            y = y_bs,                               # Input Bootstrap DV
            X = X_bs,                               # Input Bootstrap Intercept and IVs
            control = control,                      # Define Control-Option
            method = method,                        # Define Climbing-Option
            lower = lower)
          return(res_bs$par)}
      }

      bs_res <- t(replicate(bs.no, mle_bs()))     # Repeat this procedure bs.no times
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

      res$betas_bs <- res$table_bs[1:ncol(X),]      # Exclude Sigma^2-Row
      rownames(res$betas_bs) <- rownames(res$betas) # Give Table Row-Names
      colnames(res$betas_bs)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly

      if(hetero){
        res$gammas_bs <- res$table_bs[(ncol(X)+1):(ncol(X)+ncol(Z)),] # Take only Gamma-Coefficients
        rownames(res$gammas_bs) <- c("Int.", ZV)         # Give Table Row-Names
        colnames(res$gammas_bs)[3:4] <- colnames(res$betas)[3:4] # Label CIs correctly
        res$meanvar <- apply(exp(Z %*% res$par[(ncol(X)+1):(ncol(X)+ncol(Z))]),2,mean)
      }
    }

    ##### Format Output ---
    if (missing(digits)) {                        # Prints Coefficient-Table,
      cat("Beta-Coefficients \n")                 # Log.Likelihood and Variance &
      print(res$betas)                            # saves output
      cat("\n")
      if(hetero){
        cat("Gamma-Coefficients \n")
        print(res$gammas)
        cat("\n")}
      if(bootstrap){
        cat("Bootstrapping \n")
        print(res$betas_bs)
        cat("\n")
        if(hetero){
          print(res$gammas_bs)
          cat("\n")}}
      cat("Log. Likelihood:", res$value, "\n")
      if(ll=="linear"){
        if(hetero){
          cat("Mean Var.", res$meanvar)
        } else{ cat("Variance:", res$par[length(res$par)])}
      } else {
        cat("AIC:", res$aic)
        cat("\n")
        cat("BIC:", res$bic)
      }
      invisible(res)
    }


    else{
      cat("Beta-Coefficients \n")                 # Prints Coefficient-Table,
      print(apply(res$betas, c(1,2),round,digits)) # Log.Likelihood and Variance
      cat("\n")                                   # rounded to digits-argument &
      if(hetero){                                 # saves output
        cat("Gamma-Coefficients \n")
        res$gammas <- apply(res$gammas, c(1,2),round,digits)
        print(res$gammas)
        cat("\n")}
      if(bootstrap){
        cat("Bootstrapping \n")
        print(apply(res$betas_bs, c(1,2),round,digits))
        cat("\n")
        if(hetero){
          print(apply(res$gammas_bs, c(1,2),round,digits))
          cat("\n")}}
      cat("Log. Likelihood:", round(res$value,digits), "\n")
      if(ll=="linear"){
        if(hetero){
          cat("Mean Var.", round(res$meanvar, digits))
        } else{cat("Variance:", round(res$par[length(res$par)], digits))}
      }
      else {
        cat("AIC:", round(res$aic, digits))
        cat("\n")
        cat("BIC:", round(res$bic, digits))
      }
      invisible(res)
    }
  }

