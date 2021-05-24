library(dplyr, warn.conflicts = FALSE)
library(haven)

url <- paste0("https://wps.pearsoned.com/wps/media/objects/11422/",
              "11696965/datasets3e/datasets/hmda_aer.dta")

hdma <- 
  read_dta("~/Downloads/hmda_aer.dta") %>%
  mutate(deny = case_when(s7 == 3 ~ TRUE,
                          s7 %in% c(1, 2) ~ FALSE,
                          TRUE ~ NA),
         black = s13 == 3,
         lwage = if_else(s31a > 0 & s31a < 999999, log(s31a), NA_real_),
         emp = if_else(s25a > 1000, NA_real_, s25a),
         married = s23a == "M",
         across(c(s11, s31c, s35, s6, s50, s33, school, s57, s48, s39,
                  chval, s20),
                ~ if_else(. > 999999, NA_real_, .))) %>%
  rename(mhist = s42, chist = s43, phist = s44)
         

lm(deny ~ black, data = hdma)

hdma_cc <- 
  hdma %>%
  select(deny, black, lwage, chist, mhist, phist, emp) %>%
  na.omit()

Y2 <- 
  hdma_cc %>% 
  select(deny) 

X2 <- 
  hdma_cc %>%
  select(black) %>%
  mutate(constant = 1)
  
W2 <- 
  hdma_cc %>% 
  select(-deny, -black) %>% 
  mutate(constant = 1)
  
get_coefs <- function(X, Y) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  solve(t(X) %*% X) %*% t(X) %*% Y
}
 
beta_tilde_hat <- get_coefs(X2, Y2)
Delta_hat <- get_coefs(X2, W2) 
gamma_hat <- get_coefs(W2, Y2)
beta_tilde_hat - Delta_hat %*% gamma_hat
