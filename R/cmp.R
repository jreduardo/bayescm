# @title COM-Poisson Model (mean parametrization) in JAGS Language
# @author Eduardo Jr
# @description A string with defined model in a JAGS language
jags_string_cmp <- function() {
    "
 data {
    for (i in 1:n) {
      ze[i] <- 0
    }
  }
  model {
    # Constant to zero trick (large enough to ensure la[i] > 0)
    C <- 1e5
    # Log-verossimilhanca
    for (i in 1:n) {
      lambda[i] <- (mu[i] + (nu - 1)/(2 * nu))^nu
      for (j in 1:sumto) {
        Z[i, j] <- exp(j * log(lambda[i]) - nu * logfact(j))
      }
      la[i] <- -ll[i] + C
      ll[i] <- y[i] * log(lambda[i]) - nu * logfact(y[i]) - log(sum(Z[i, ]))
      mu[i] <- exp(inprod(X[i, ], theta[]))
      ze[i] ~ dpois(la[i])
    }
    # Prioris
    for (j in 1:p) {
      theta[j] ~ dnorm(0.0, 1e-3)
    }
    phi ~ dnorm(0.0, 0.1)
    nu <- exp(phi)
  }
"
}

#' @title Fitting COM-Poisson (mean parametrization) Model Based on MCMC
#' @author Eduardo Jr
#' @description Build matrices and organize data to use
#'     \code{\link[rjags]{jags.model}}
#'     and\code{\link[rjags]{coda.samples}} for analysis of Bayesian
#'     models using Markov Chain Monte Carlo (MCMC).
#' @param formula A formula to define fixed effects.
#' @param data The data frame.
#' @param sumto Upper bound to infinite sum of the normalization
#'     constant
#' @param .control_model See \link{control_model}.
#' @param .control_samples See \link{control_samples}.
#' @export
#' @examples
#'
#' \donttest{
#'
#' # Simulate data
#' beta <- c(3, -1)
#' X <- cbind(1, runif(50))
#' y <- rpois(50, lambda = exp(X %*% beta))
#'
#' # Sampling and summarise posterior
#' model <- bayescmp(y ~ X - 1, sumto = 100L)
#' vapply(model$samples, function(x) apply(x, 2L, mean), double(3))
#'
#' }
#'
bayescmp <- function(formula, data, sumto,
                     .control_model = control_model(),
                     .control_samples = control_samples()) {
    #-------------------------------------------
    # Build matrices
    if (missing(data)) data <- environment(formula)
    frame <- stats::model.frame(formula, data)
    terms <- attr(frame, "terms")
    y <- stats::model.response(frame)
    X <- stats::model.matrix(terms, frame)
    yy <- y
    if (any(y == 0)) yy[y == 0] <- 0.01
    #-------------------------------------------
    # Define data
    data_jags <- list("n" = nrow(X),
                      "p" = ncol(X),
                      "X" = X,
                      "y" = yy,
                      "sumto" = sumto)
    #-------------------------------------------
    # Settings for MCMC
    if (is.null(.control_model$inits)) {
        m0 <- stats::glm.fit(x = X, y = y, family = stats::poisson())
        coefs <- m0$coefficients
        .control_model$inits <-  list("theta" = coefs, "phi" = 0)
    }
    if (is.null(.control_samples$variable.names)) {
        .control_samples$variable.names <- names(.control_model$inits)
    }
    #-------------------------------------------
    # Initialize model
    model <- rjags::jags.model(
        file     = base::textConnection(jags_string_cmp()),
        data     = data_jags,
        inits    = .control_model$inits,
        n.chains = .control_model$n.chains,
        n.adapt  = .control_model$n.adapt,
        quiet    = .control_model$quiet)
    #-------------------------------------------
    # MCMC Sampling
    samples <- rjags::coda.samples(
        model          = model,
        variable.names = .control_samples$variable.names,
        n.iter         = .control_samples$n.iter,
        thin           = .control_samples$thin,
        na.rm          = .control_samples$na.rm)
    samples <- coda::as.mcmc.list(
        lapply(samples, function(x) {
            pnames <- colnames(x)
            colnames(x)[grepl("theta", pnames)] <- colnames(X)
            coda::as.mcmc(x)
        }))
    #-------------------------------------------
    # Output
    output <- list(
        "model" = "COM-Poisson",
        "formula" = formula,
        "model" = model,
        "posterior_samples" = samples
    )
    class(output) <- "bayescm"
    return(output)
}
