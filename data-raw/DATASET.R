## code to prepare `DATASET` dataset goes here

set.seed(100)
X_a <- rnorm(500, 2, 1)
Q_a <- rnorm(500, 2, 1)
Y_0_a <- (1/4)*sqrt(abs(X_a)) + (1/4)*sqrt(abs(Q_a))
tau_a <- (1/6)*X_a + (1/3)*(exp(X_a))^(1/4) + (1/6)*Q_a
D_a <- rbinom(n=500, size=1, prob=1/(1+1/exp(-10+12*tau_a)))

Y_a <- Y_0_a + tau_a*D_a

X_b <- rnorm(500, 0, 1)
Q_b <- rnorm(500, 0, 1)
Y_0_b <- (1/4)*sqrt(abs(X_b)) + (1/4)*sqrt(abs(Q_b))
tau_b <- (1/10)*X_b + (exp(X_b))^(1/8) + (1/10)*Q_b
D_b <- rbinom(n=500, size=1, prob=1/(1+1/exp(7-10*tau_b)))
Y_b <- Y_0_b + tau_b*D_b

exp_data <- as.data.frame(cbind(c(Y_a,Y_b),c(D_a,D_b),c(X_a,X_b),c(Q_a,Q_b)))
colnames(exp_data) <- c("outcome","treatment","confounder","Q")
exp_data$group_a <- c(rep(1,500),rep(0,500))
exp_data$group_b <- c(rep(0,500),rep(1,500))

usethis::use_data(exp_data, overwrite = TRUE, internal = FALSE)
