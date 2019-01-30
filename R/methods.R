
da <- read.table("~/GitProjects/article-reparcmp/data/nitrofen.txt",
                 header = TRUE)
str(da)

model0 <- bayesgct(novos ~ dose, data = da,
                   .control_model = control_model(n.chains = 1),
                   .control_samples = control_samples(n.iter = 1000))

model1 <- bayesgct(novos ~ dose, data = da,
                   .control_samples = control_samples(n.iter = 1000))

model2 <- bayescmp(novos ~ dose, data = da, sumto = 100L,
                   .control_samples = control_samples(n.iter = 1000))

model3 <- bayescmp0(novos ~ dose, data = da, sumto = 100L,
                    .control_samples = control_samples(n.iter = 1000))


bmark <- function(...) microbenchmark::microbenchmark(...)
glm(novos ~ dose, data = da)

str(model1$samples)

library(rjags)
summary(model1$samples)

x <- model1
str(x)



p0 <- dpois(0, 1/7)
p1 <- 1 - dpois(0, 1/7)

p0^2 + p1^2 + 2*p0*p1

p0^4*p1^3+
p1^4*p0^3


#-------------------------------------------
# Functions
all_comb <- function(x, n) {
    out <- do.call(expand.grid, lapply(1:n, function(i) x))
    attr(out, "out.attrs") <- NULL
    names(out) <- paste0("y", 1:n)
    out
}
verify_consecutive <- function(x, n = length(x)) {
    for (i in 2:n) {
        if((x[i] + x[i - 1]) == 0) return(TRUE)
    }
    return(FALSE)
}
simulate_problem <- function(rate = 1/7, ndays = 7) {
    x <- rpois(ndays, rate)
    verify_consecutive(x, ndays)
}

#-------------------------------------------
# Build all possibilities
da <- all_comb(0:1, 7)

# Compute probabilities for each possibility
da$prob <- apply(da, 1, function(x) {
    prod(dbinom(x, size = 1, prob = 1 - dpois(0, 1/7)))
})

# Number of days where Y_i=0
da$z <- apply(da[, 1:7], 1, sum)

# Verify the sequences that has consecutive 0 values
da$consecutive <- apply(da, 1, verify_consecutive)

# Compute the probability
res <- 1 - with(da, sum(prob[consecutive]))
res

#-------------------------------------------
# Simulate the problem to confirm
consegue <- replicate(1e6, simulate_problem())
sim <- 1 - sum(consegue)/length(consegue)
sim


#=======================================================================

round(aggregate(prob ~ z, data = da, unique), 6)
round(cbind(0:7, dbinom(0:7, 7, 1 - dpois(0, 1/7))), 6)

table(subset(da, consecutive == 1)$z)

arr <- function(n, k) {
    exp(lfactorial(n) - lfactorial(n - k))
}

arr(7, 1)
arr(6, 2)
arr(7, 3)

choose(7, 1)
choose(7, 2)

choose(7, 5)


sum(da$prob)

xx <- all_comb(0:1, 2)
xx$Z <- apply(xx, 1, sum)
xx$py <- apply(xx, 1, function(i) prod(dbinom(i, 1, 0.3)))

aggregate(py ~ Z, data = xx, sum)
cbind(0:max(xx$Z), dbinom(0:max(xx$Z), max(xx$Z), 0.3))

dbinom(0, 1000, 0.0001)
dbinom(0, 1000, 0.0001)




1 - dbinom(0, 10, 0.05)
dbinom(10, 10, 0.05)
10*0.05
10*0.05*0.95

bmark(dpois(0, 1000 * 0.0001), dbinom(0, 1000, 0.0001))

# Print
x <- model0

# Summary
summary.bayescm <- function(object, ...) {
    samples <- do.call("rbind", object$posterior_samples)
    vstats <- apply(samples, 2, function(x) {
        c("Mean" = mean(x), "SD" = sd(x),
          quantile(x, 0.025), quantile(x, 0.975))
    })
    out <- list("statistics" = t(vstats),
                "formula" = object$formula,
                "model" = object$model,
                "chains" = length(object$posterior),
                "samples" = nrow(object$posterior[[1L]]))
    class(out) <- "summary.bayescm"
    return(out)
}

# Print summary
print.summary.bayescm <- function(x, ...) {
    npars <- nrow(x$statistics)
    cat("\nModel: ", sprintf("%s (%s parameters)", x$model, npars),
        "\n", sep = "")
    cat("  Call     => ", deparse(x$formula), "\n", sep = "")
    cat("  Samples  => ", x$samples, " x ", x$chains, "\n", sep = "")
    cat("  Burnin   => ", "NEED", "\n", sep = "")
    cat("  Thinning => ", "NEED", "\n\n", sep = "")
    print(x$statistics)
    invisible(x)
}



model0

x <- summary.bayescm(model0)
x

object <- model0
str(object)





coda:::summary.mcmc.list(x$posterior_samples)$statistics
coda:::summary.mcmc.list(model0$posterior_samples)
coda:::print.summary.mcmc

x$formula

nchains <- length(x$samples)
niter <- nrow(x$samples[[1]])



nchain

start(x)

coda:::summary.mcmc.list
coda:::print.summary.mcmc

#=======================================================================

# Simulate Mixed Poisson Model
set.seed(10209770)

n <- 1000
x <- seq(0, 1, length.out = n)
g <- factor(rep(1:3, length.out = n))
da <- data.frame(x = x, g = g)

beta <- rbind(1, -1.5)
sigma <- 1

X <- model.matrix(~x)
Z <- model.matrix(~g - 1)
b <- rnorm(nlevels(g), sd = sigma)

Xb <- X %*% beta
Zb <- Z %*% b
mu <- exp(Xb + Zb)

y <- rpois(n, lambda = mu)
da <- data.frame(y = y, x = x, g = g)

library(lattice)
xyplot(y ~ x, type = c("p", "smooth"), groups = g, data = da)



#=======================================================================
library(tibble)
set.seed(10209770)

# Covariable and grouping factor
n <- 24
x <- seq(0, 1, length.out = n)
g <- factor(rep(LETTERS[1:3], length.out = n))
da <- tibble(x = x, g = g)
da


# Parameters
beta <- rbind(1, -1.5)
sigma_0 <- 3
sigma_1 <- 1
r <- 0.5
covar <- r * sigma_0 * sigma_1
Sigma <- matrix(c(sigma_0^2, covar, covar, sigma_1^2), ncol = 2)

# Generate matrices
form <- ~x + (x | g)
mats <- build_matrices(form, data = da)
X <- mats$X
Z <- mats$Z


b <- c(mvrnorm(3, c(0, 0), Sigma))


Z %*% b

b <- rnorm(nlevels(g), sd = sigma)

Xb <- X %*% beta
Zb <- Z %*% b
mu <- exp(Xb + Zb)

y <- rpois(n, lambda = mu)
da <- data.frame(y = y, x = x, g = g)

library(lattice)
xyplot(y ~ x, type = c("p", "smooth"), groups = g, data = da)

#=======================================================================
poisRE0 <-
" model {
    # Log-verossimilhanca
    for (i in 1:n) {
      mu[i] <- exp(X[i, ] %*% beta[] + Z[i, ] %*% b[])
      y[i] ~ dpois(mu[i])
    }
    # Random effects (uncorrelated random-effects)
    for (i in 1:n_re) {
      for (j in 1:n_gr) {
      b[k] ~ dnorm(0.0, tau)
      }
    }
    # Prioris
    tau ~ dgamma(0.001, 0.001)
    for (j in 1:p) {
      beta[j] ~ dnorm(0.0, 1e-3)
    }
  }
"



data0 <- with(
    da,
    list("y" = y,
         "n" = length(y),
         "m" = length(unique(g)),
         "ind" = as.numeric(cult))
)

jagsmodel0 <- jags.model(
    textConnection(poisRE0),
    data = data0,
    n.chains = 3,
    n.adapt = 1000
)

amostra0 <- coda.samples(
    jagsmodel0, c("b0", "sigma", "u", "mu"),
    n.iter = 10000, thin = 10,
    n.chains = 3,
    n.adapt = 1000)


#=======================================================================

library(lme4)
m0 <- glmer(y ~ x + (x | g), family = poisson, data = da)

ranef(m0)
plot(ranef(m0)$g[, 1], ranef(m0)$g[, 2])

anova(m0, m1)
coef(m0)

mform <- y ~ x + (x | g)
mform <- y ~ x + (x | g) + (1 | rep(1:2, 500))


(bar.f <- findbars(mform)) # list with 3 terms
mf <- stats::model.frame(subbars(mform), data = da)
rt <- mkReTrms(bar.f,mf)
names(rt)
rt$theta
rt$Lambdat
lme4::nobars(mform)


build_matrices <- function(formula, data) {
    sform <- lme4::subbars(term = formula)
    rform <- lme4::findbars(term = formula)
    fform <- lme4::nobars(term = formula)
    mframe <- stats::model.frame(formula = sform, data = data)
    relist <- lme4::mkReTrms(bars = rform, fr = mframe)
    if (length(relist$flist) > 1) {
        stop("Models with more than one grouping factor have not yet been implemented.")
    }
    X <- stats::model.matrix(object = fform, data = data)
    Z <- t(relist$Zt)
    out <- list(X = X, Z = Z)
    return(out)
}


split(model.response(mframe), relist$flist[[1]])




mats <- build_matrices(formula = y ~ x + (x || g), data = da)
head(mats$X)
head(mats$Z)

colnames(t(rt$Zt))


m1 <- glmer(y ~ x + (x || g), family = poisson, data = da)
getME(m0, "Z")
getME(m1, "theta")

t(KhatriRao(t(Z), t(X)))

names(lmod)
str(lmod$fr)
names(lmod$reTrms)

lmod <- lFormula(y ~ x + (x || g), da)
lmod <- lFormula(y ~ x + (x | g), da)
lmod <- lFormula(y ~ x + (x | g) + (1 | rep(1:2, 500)), da)

head(t(lmod$reTrms$Zt))
lmod$reTrms$theta
lmod$reTrms$Lind
lmod$reTrms$Gp
lmod$reTrms$lower
head(t(lmod$reTrms$Lambdat))
length(lmod$reTrms$flist)
lmod$reTrms$cnms

lmod$reTrms$flist[[1]]


#-------------------------------------------

library(MASS)
library(lme4)

generate_data = function(
	n # number of units
	, k # number of trials within each condition within each unit
	, noise # measurement noise variance
	, I # population intercept
	, vI # across-units variance of intercepts
	, A # population A effect
	, vA # across-units variance of A effects
	, rIA # across-units correlation between intercepts and A effects
){
	Sigma = c(
		vI , sqrt(vI*vA)*rIA
		, sqrt(vI*vA)*rIA , vA
	)
	Sigma = matrix(Sigma,2,2)
	means = mvrnorm(n,c(I,A),Sigma)
	temp = expand.grid(A=c('a1','a2'),value=0)
	temp$A = factor(temp$A)
	contrasts(temp$A) = contr.sum
	from_terms = terms(value~A)
	mm = model.matrix(from_terms,temp)
	data = expand.grid(A=c('a1','a2'),unit=1:n,trial=1:k,
                           KEEP.OUT.ATTRS = FALSE)
	for(i in 1:n){
		data$value[data$unit==i] = as.numeric(mm %*% means[i,]) + rnorm(k*2,0,sqrt(noise))
	}
	data$unit = factor(data$unit)
	data$A = factor(data$A)
	contrasts(data$A) = contr.sum
	return(data)
}

this_data = generate_data(
	n = 100 # number of units
	, k = 100 # number of trials within each condition within each unit
	, noise = 1 # measrurement noise variance
	, I = 2 # population intercept
	, vI = 3 # across-units variance of intercepts
	, A = 4 # population A effect
	, vA = 5 # across-units variance of A effects
	, rIA = .6 # across-units correlation between intercepts and A effects
)

fit = lmer(
	data = this_data
	, formula = value ~ (1+A|unit) + A
)

plot(ranef(fit)[[1]])
print(fit)

#-------------------------------------------
