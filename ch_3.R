set.seed(123456789)

get_data <- function() {
  N <- 1000
  a <- 2
  b <- 3
  c <- 2
  e <- 3
  f <- -1
  d <- 4
  z <- runif(N)
  u_1 <- rnorm(N, mean = 0, sd = 3)
  u_2 <- rnorm(N, mean = 0, sd = 1)
  
  x <- f + d * z + u_2 + c * u_1
  y <- a + b * x + e * u_1
  data.frame(x, y, z)
}

df <- get_data()

lm1 <- lm(y ~ x, data = df)

lm1

X <- cbind(1, df$x)
Z <- cbind(1, df$z)
y <- df$y

beta_hat_ols <- solve(t(X) %*% X) %*% t(X) %*% y
beta_hat_iv <- solve(t(Z) %*% X) %*% t(Z) %*% y

beta_hat_ols
beta_hat_iv

lm_iv <- function(y, X_in, Z_in = X_in, reps = 100,
                  min_in = 0.05, max_in = 0.95) {
  
  set.seed(123456789)
  X <- cbind(1, X_in)
  Z <- cbind(1, Z_in)
  
  bs_mat <- matrix(NA, reps, dim(X)[2])
  
  N <- length(y)
  for (r in 1:reps) {
    # sample(1:N, size = N, replace = TRUE)
    index_bs <- round(runif(N, min = 1, max = N))
    y_bs <- y[index_bs]
    X_bs <- X[index_bs, ]
    Z_bs <- Z[index_bs, ]
    bs_mat[r, ] <- solve(t(Z_bs) %*% X_bs) %*% t(Z_bs) %*% y_bs
  }

  # Present results
  tab_res <- matrix(NA, dim(X)[2], 4)
  tab_res[, 1] <- colMeans(bs_mat)
  for (j in 1:dim(X)[2]) {
    tab_res[j, 2] <- sd(bs_mat[ , j])
    tab_res[j, 3] <- quantile(bs_mat[ , j], min_in)
    tab_res[j, 4] <- quantile(bs_mat[ , j], max_in)
  }
  colnames(tab_res) <- c("coef", "sd", as.character(min_in),
                         as.character(max_in))
  return(tab_res)
  
}

print(lm_iv(df$y, df$x), digits = 3)
print(lm_iv(df$y, df$x, df$z), digits = 3)

# 3.4: Returns to schooling ----

