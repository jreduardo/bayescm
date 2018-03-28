#' @name bayescm_control
#' @author Eduardo Jr
#' @title Controls Initialization of JAGS Model and MCMC Sampling
#' @description Parameters to control initialization of JAGS model,
#'     \code{\link[rjags]{jags.model}} and MCMC sampling,
#'     \code{\link[rjags]{jags.model}}.
#' @param ... Arguments passed to \code{\link[rjags]{jags.model}} and
#'     \code{\link[rjags]{coda.samples}}.See arguments of both function.
#' @return An object of class \code{bayescm_control} to use in
#'     \link{bayesgct} and \link{bayescmp}.
#'
NULL

#' rdname bayescm_control
#' @param ... Arguments passed to \code{\link[rjags]{jags.model}}. See
#'     the function arguments.
#' @export
control_model <- function(...) {
    default <- list(inits    = NULL,
                    n.chains = 3L,
                    n.adapt  = 1000L,
                    quiet    = FALSE)
    userval <- list(...)
    if (!all(names(userval) %in% names(default))) {
        txt <- "Incorrect values. Values must be options for rjags::jags.model."
        stop(txt)
    }
    output <- utils::modifyList(default, userval)
    attr(output, "jagsfun") <- "jags.model"
    class(output) <- "bayescm_control"
    return(output)
}

#' rdname bayescm_control
#' @param ... Arguments passed to \code{\link[rjags]{coda.samples}}. See
#'     the function arguments.
#' @export
control_samples <- function(...) {
    default <- list(variable.names = NULL,
                    n.iter         = 10000,
                    thin           = 1,
                    na.rm          = TRUE)
    userval <- list(...)
    if (!all(names(userval) %in% names(default))) {
        txt <- "Incorrect values. Values must be options for rjags::jags.model."
        stop(txt)
    }
    output <- utils::modifyList(default, userval)
    attr(output, "jagsfun") <- "coda.samples"
    class(output) <- "bayescm_control"
    return(output)
}

# @title Print method for MCMC controls
print.bayescm_control <- function(x, ...) {
    x[vapply(x, is.null, logical(1))] <- "'default'"
    title <- sprintf(
        "MCMC controls (arguments for rjags::%s() function)\n",
        attr(x, "jagsfun"))
    out <- sprintf("  %s = %s",
                   format(names(x), width = 15),
                   format(unlist(x), justify = "centre"))
    cat(title, "\n")
    cat(out, sep = "\n")
    invisible()
}
