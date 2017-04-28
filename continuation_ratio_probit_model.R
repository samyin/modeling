chemo <- read.csv("chemotherapy.csv")

d <- length(unique(chemo$Response))
q <- d - 1 + ncol(chemo) - 1
n <- nrow(chemo)
N <- sum(chemo$Response) - sum(chemo$Response == max(chemo$Response))

y <- vector(length = N)
W <- matrix(nrow = N, ncol = d - 1)
x <- matrix(nrow = N, ncol = ncol(chemo))
s <- 0
for (i in 1:n) {
  if (chemo$Response[i] == 1) {
    y[s + 1] <- 1
    W[s + 1, 1] <- 1
    W[s + 1, -1] <- rep(0, d - 2)
    x[s + 1, 1] <- i
    s <- s + 1
  } else if (chemo$Response[i] == d) {
    y[(s + 1):(s + d - 1)] <- 0
    W[(s + 1):(s + d - 1), ] <- diag(d - 1)
    x[(s + 1):(s + d - 1), 1] <- rep(i, d - 1)
    s <- s + d - 1
  } else {
    r <- chemo$Response[i]
    y[(s + 1):(s + r - 1)] <- 0
    y[s + r] <- 1
    W[(s + 1):(s + r), ] <- diag(d - 1)[1:r, ]
    x[(s + 1):(s + r), 1] <- rep(i, r)
    s <- s + r
  }
}
for (j in 1:N) {
  for (k in 1:(ncol(chemo) - 1)) {
    x[j, k + 1] <- chemo[x[j, 1], k]
  }
}
data_trans <- as.data.frame(cbind(W, x[, -1], y))
colnames(data_trans) <- c(paste0("ind_", 1:(d - 1)), 
                          "Therapy", "Gender", "Response")
rm(W, x, y, i, j, k, r, s)

# Frequentist GLM
model_1 <- glm(Response ~ -1 + ., data = data_trans, 
               family = binomial(link = "probit"))
summary(model_1)

# Gibbs Sampling
library(truncnorm)
library(MASS)
library(coda)

x <- as.matrix(data_trans[, -ncol(data_trans)])
y <- data_trans[, ncol(data_trans)]
t <- 10000

z <- matrix(nrow = nrow(data_trans), ncol = t)
b <- matrix(nrow = ncol(data_trans) - 1, ncol = t)
Sigma0 <- diag(nrow(b))
Sigma <- solve(solve(Sigma0) + t(x) %*% x)
mu0 <- as.vector(rep(0, nrow(b)))
mu <- matrix(nrow = nrow(b), ncol = t)

z[, 1] <- rep(0, nrow(z))
b[, 1] <- rep(0, nrow(b))
mu[, 1] <- rep(0, nrow(mu))

for (i in 1:(t - 1)) {
  for (j in 1:nrow(data_trans)) {
    if (y[j] == 0) {
      z[j, i + 1] <- rtruncnorm(1, a = -Inf, b = 0, 
                                mean = x[j, ] %*% b[, i], sd = 1)
    } else {
      z[j, i + 1] <- rtruncnorm(1, a = 0, b = Inf, 
                                mean = x[j, ] %*% b[, i], sd = 1)
    }
  }
  mu[, i + 1] <- Sigma %*% (t(x) %*% z[, i + 1])
  b[, i + 1] <- mvrnorm(1, mu[, i + 1], Sigma)
}

theta_sample <- as.mcmc(t(b[, 2001:10000]))
colnames(theta_sample) <- c(paste0("ind_", 1:(d - 1)), "Therapy", "Gender")
summary(theta_sample)