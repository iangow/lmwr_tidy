library(googlesheets4)
library(dplyr, warn.conflicts = FALSE)
library(stargazer)

#nlsm_raw <- sheets_read("116KwP9SJUQ8LZHeIH-sTvlFFwCgkKIsz86exk1I7O1s",
#                 na = ".")

url <- "https://davidcard.berkeley.edu/data_sets/proximity.zip"
t <- tempfile(fileext = ".zip")
download.file(url, t)

nlsm <-
  nlsm %>%
  mutate(exp = age76 - ed76 - 6,
         exp2 = exp^2/100) %>%
  filter(!is.na(lwage76))



fm <- list()
fm[[1]] <- lm(lwage76 ~ ed76, data = nlsm)
fm[[2]] <- update(fm[[1]], ~ . + exp + exp2)
fm[[3]] <- update(fm[[2]], ~ . + black + reg76r)
fm[[4]] <- update(fm[[3]], ~ . + smsa76r + smsa66r + reg662 + reg663 +
                    reg664 + reg665 + reg666 + reg667 + reg668 + reg669)
fm[[5]] <- update(fm[[4]], ~ . - ed76 + nearc4)
fm[[6]] <- update(fm[[5]], ed76 ~ .)

stargazer(fm, type = "text", keep.stat = c("n", "rsq"),
          omit = "^(reg66|smsa)")
fm[[5]]$coefficients[2]/fm[[6]]$coefficients[2]
fm[[5]]$coefficients[["nearc4"]]/fm[[6]]$coefficients[["nearc4"]]

iv_fit <- function(y_in, X_in, Z_in = X_in) {
  # Convert supplied data to matrices
  y <- as.matrix(y_in, ncol=1)
  X <- cbind(1, as.matrix(X_in))
  Z <- cbind(1, as.matrix(Z_in))

  res <- solve(t(Z) %*% X) %*% t(Z) %*% y
  rownames(res) <- c("intercept", colnames(X_in))
  return(res)
}

iv_bs <- function(i, y, X_in, Z_in) {

  y_vec <- as.matrix(y, ncol=1)
  N <- length(y_vec)

  index_bs <- round(runif(N, min = 1, max = N))
  y_bs <- y_vec[index_bs, ]
  X_bs <- X_in[index_bs, ]
  Z_bs <- Z_in[index_bs, ]

  as_tibble(t(iv_fit(y_bs, X_bs, Z_bs)))
}

lm_iv <- function(y, X_in, Z_in = X_in, reps = 100) {

  set.seed(123456789)
  bind_rows(lapply(1:reps, iv_bs, y, X_in, Z_in))
}

get_summ_stats <- function(x, p_lower, p_upper) {
  df_temp <- tibble(mean(x),
                    sd(x),
                    quantile(x, p_lower),
                    quantile(x, p_upper))

  names(df_temp) <- c("coef", "sd",
                         paste0("p", as.character(p_lower * 100)),
                         paste0("p", as.character(p_upper * 100)))
  return(df_temp)
}

get_tab_res <- function(bs_mat, p_lower = 0.05, p_upper = 0.95) {
  bind_rows(lapply(bs_mat, get_summ_stats, p_lower, p_upper),
            .id = "var")
}

y <-
  nlsm %>%
  select(lwage76)

X <-
  nlsm %>%
  select(ed76, exp, exp2, black, reg76r,
                     smsa76r, smsa66r, reg662:reg669)
Z1 <-
  nlsm %>%
  mutate(age2 = age76^2) %>%
  select(nearc4, age76, age2, black, reg76r,
         smsa76r, smsa66r, reg662:reg669)

Z2 <-
  nlsm %>%
  mutate(age2 = age76^2) %>%
  select(momdad14, age76, age2, black, reg76r,
         smsa76r, smsa66r, reg662:reg669)

res_Z1 <- lm_iv(y, X, Z1, reps = 1000)
get_tab_res(res_Z1)

res_Z2 <- lm_iv(y, X, Z2, reps = 1000)
bs_diff <- (res_Z1 - res_Z2)[, "ed76"]
summary(bs_diff)

# 3.4.5 Concerns with Distance to College
summarise_by <- function(.data, .by, ..., .groups = NULL) {
  .data %>%
    summarise(across(c(ed76, exp, black, south66,
                       smsa66r, reg76r, smsa76r), mean),
              .groups = "keep") %>%
    pivot_longer(cols = everything())
    pivot_wider(names_from = .by, values_from = value)
}

nlsm %>%
  mutate(nearc4 = case_when(nearc4 == 0 ~ "Not near college",
                            nearc4 == 1 ~ "Near college")) %>%
  group_by(nearc4) %>%
  summarise_by(.by = "nearc4",
               across(c(ed76, exp, black, south66,
                     smsa66r, reg76r, smsa76r), mean))


