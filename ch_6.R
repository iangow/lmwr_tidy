# 6.2.2. Simulation of censored data ----

set.seed(123456789)

N <- 500
a <- 2
b <- -3
x <- runif(N)
u <- rnorm(N)
y_star <- a + b * x + u
# y <- pmax(y_star, 0)
y <- ifelse(y_star > 0, y_star, 0)
lm1 <- lm(y ~ x)
lm1$coefficients[2]

lm(y[x < 0.6] ~ x[x < 0.6])$coefficients[2]

length(x[x < 0.6])
sum(x < 0.6)

glm(y > 0 ~ x, family = binomial(link = "probit"))$coefficients

# 6.2.5 Tobit estimator in R ----
f_tobit <- function(par, y, X) {
  X <- cbind(1, X)
  sigma <- exp(par[1]) # Use exp() to keep value positive
  beta <- par[-1]
  
  z <- (y - X %*% beta)/sigma
  
  log_lik <- sum((y == 0) * log(pnorm(z)) +
                    (y > 0) * (log(dnorm(z)) - log(sigma)))
  
  # Return negative because we are minimizing
  return(-log_lik)
}
    
par1 <- c(0, lm1$coefficients)
a1 <- optim(par = par1, fn = f_tobit, y = y, X = x)

# sigma
exp(a1$par[1])

a1$par[-1]

# Check against AER
library(AER)
tobit(y ~ x)

# 6.4.2 Simulation of a selection model ----

library(mvtnorm)
set.seed(123456789)
N <- 100
a <- 6
b <- -3
c <- 4
d <- -5

mu <- c(0, 0)
Sigma <- rbind(c(1, 0.5), c(0.5, 1))


df <- 
  tibble(x = runif(N),
         z = runif(N),
         u = rmvnorm(N, mean = mu, sigma = Sigma),
         y = ifelse(c + d*z + u[, 1] > 0, a + b * x + u[, 2], 0))

# OLS model
lm1 <- lm(y ~ x, data = df)
lm2 <- update(fm1, subset = z < 0.6)
glm1 <- glm(y > 0 ~ z, family = binomial(link = "probit"))
stargazer(list(lm1, lm2, glm1), type = "text", keep.stat = c("n", "rsq"))

# 6.4.5 Heckman estimator in R ----

f_heckman <- function(par, y, X_in, Z_in = X_in) {
  
  X <- cbind(1, X_in)
  Z <- cbind(1, Z_in)
  rho <- exp(par[1])/(1 + exp(par[1]))
  
  # This is the sigmoid function.
  # Note that in fact rho is between -1 and 1
  k <- dim(X)[2]
  beta <- par[2:(2 + k -1)]
  gamma <- par[(2 + k):length(par)]
  
  Xb <- X %*% beta
  Zg <- Z %*% gamma
  
  Zg_adj <- (Zg + rho * (y - Xb))/((1 - rho^2)^(.5))
  log_lik <- (y == 0) * log(pnorm(-Zg)) +
    (y > 0) * (log(pnorm(Zg_adj)) + log(dnorm(y - Xb)))
  
  return(-sum(log_lik))
}

f_heckman_alt <- function(par, y, X_in, Z_in = X_in) {
  
  X <- cbind(1, X_in)
  Z <- cbind(1, Z_in)
  
  # This is the sigmoid function.
  # Note that in fact rho is between -1 and 1
  rho <- exp(par[1])/(1 + exp(par[1]))
  sigma <- exp(par[2])
  
  k <- dim(X)[2]
  beta <- par[3:(3 + k - 1)]
  gamma <- par[(3 + k):length(par)]
  
  Zg <- Z %*% gamma
  Xb <- X %*% beta
  u2 <- y - Xb

  B <- (Zg + rho/sigma*u2)/sqrt(1 - rho^2)
  
  log_lik <- ifelse(y == 0, pnorm(-Zg, log.p=TRUE),
                    -1/2*log(2*pi) - log(sigma) +
                      (pnorm(B, log.p=TRUE) - 0.5*(u2/sigma)^2))
  
  return(-sum(log_lik))
}

f_heckman_sigma_1 <- function(par, y, X_in, Z_in) {
  par_alt <- c(par[1], 0, par[-1])
  f_heckman_alt(par_alt, y, X_in, Z_in)
}

par1 <- c(0, lm1$coefficients, glm1$coefficients)
a1 <- optim(par = par1, fn = f_heckman, y = df$y, X = df$x, Z = df$z)
exp(a1$par[1])/(1 + exp(a1$par[1]))
a1$par[-2]

a3 <- optim(par = par1, fn = f_heckman_sigma_1,
            y = df$y, X = df$x, Z = df$z)
exp(a3$par[1])/(1 + exp(a3$par[1]))
a3$par[-1]

par2 <- c(0, 1, lm1$coefficients, glm1$coefficients)
a2 <- optim(par = par2, fn = f_heckman_alt,
            y = df$y, X = df$x, Z = df$z)
exp(a2$par[1])/(1 + exp(a2$par[1]))
a2$par[-1]



par2 <- c(0, 1, lm1$coefficients, glm1$coefficients)
a2 <- optim(par = par2, fn = f_heckman_alt, y = df$y, X = df$x, Z = df$z,
            control = list(maxit = 1000, reltol = 1e-15))
exp(a2$par[1])/(1 + exp(a2$par[1]))
a2$par[-1]

library(sampleSelection)
heckit <- selection((y > 0) ~ z, y ~ x, data = df, method = "ml",
                    maxMethod = "NM")
summary(heckit)
heckit$estimate

# 6.5.1 NLSY97 data ----
library(googlesheets4)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(stargazer)

nlsy97_raw <- read_sheet("1TccRf1Cbc4h9fu31WhKhzzleEKKKgiRkLofcebRQQnE",
                         na = "NA")

nlsy97 <-
  nlsy97_raw %>%
  mutate(wage = if_else(CVC_HOURS_WK_YR_ALL.07_XRND > 0,
                        `YINC-1700_2007`/CVC_HOURS_WK_YR_ALL.07_XRND, 0),
        lwage = if_else(wage > 1, log(wage), 0),
        fulltime = CVC_HOURS_WK_YR_ALL.07_XRND > 1750,
        lftwage = if_else(lwage > 0 & fulltime, lwage, 0),
        female = KEY_SEX_1997==2,
        black = KEY_RACE_ETHNICITY_1997==1,
        age = 2007 - KEY_BDATE_Y_1997,
        age2 = age^2,
        college = CV_HIGHEST_DEGREE_0708_2007 >= 3,
        south = CV_CENSUS_REGION_2007==3,
        urban = `CV_URBAN-RURAL_2007`==1,
        msa = CV_MSA_2007 > 1 & CV_MSA_2007 < 5,
        married = CV_MARSTAT_COLLAPSED_2007==2,
        children = CV_BIO_CHILD_HH_2007 > 0) %>%
  select(black, lftwage, age, age2, msa, urban, south, college, female,
         married, children) %>%
  na.omit()

nlsy97 %>%
  filter(lftwage > 0) %>%
  ggplot(aes(x = lftwage, color = female)) +
  geom_density()

lm3 <- lm(lftwage ~ female + age + age2 + black + college + south + msa,
          data = nlsy97)
lm1 <- update(lm3, subset = !female)
lm2 <- update(lm3, subset = female)

stargazer(list(lm1, lm2, lm3), type = "text", keep.stat = c("n", "rsq"))


