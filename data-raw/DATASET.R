## code to prepare `DATASET` dataset goes here

set.seed(3)
n <- 500

X_a <- rnorm(n, 1, 1.5)
Q_a <- rnorm(n, 1, 1.5)
Y_0_a <- (1/4)*sqrt(abs(X_a)) + (1/4)*sqrt(abs(Q_a))
tau_a <- (1/6)*X_a + (1/3)*(exp(X_a))^(1/4) + (1/6)*Q_a
D_a <- rbinom(n=n, size=1, prob=1/(1+1/exp(-0.3+0.8*tau_a)))


Y_a <- Y_0_a + tau_a*D_a

X_b <- rnorm(n, 0, 1.5)
Q_b <- rnorm(n, 0, 1.5)
Y_0_b <- (1/4)*sqrt(abs(X_b)) + (1/4)*sqrt(abs(Q_b))
tau_b <- (1/10)*X_b + (exp(X_b))^(1/8) + (1/10)*Q_b
D_b <- rbinom(n=n, size=1, prob=1/(1+1/exp(0.2-0.8*tau_b)))


Y_b <- Y_0_b + tau_b*D_b

exp_data <- as.data.frame(cbind(c(Y_a,Y_b),c(D_a,D_b),c(X_a,X_b),c(Q_a,Q_b)))
colnames(exp_data) <- c("outcome","treatment","confounder","Q")
exp_data$group_a <- c(rep(1,n),rep(0,n))

#mean(Y_a)-mean(Y_b)
#mean(Y_0_a)-mean(Y_0_b)
#mean(tau_b)*(mean(D_a)-mean(D_b))
#mean(D_a)*(mean(tau_a)-mean(tau_b))
#cov(D_a, tau_a)-cov(D_b, tau_b)

exp_data <- exp_data[sample(1:(2*n), 2*n, replace = FALSE),]

head(exp_data)

usethis::use_data(exp_data, overwrite = TRUE, internal = FALSE, version=3)
