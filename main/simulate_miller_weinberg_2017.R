# initialize -------------------------------------------------------------------

rm(list = ls())
devtools::load_all(".")
library(foreach)
library(magrittr)
library(codetools)
library(doParallel)
registerDoParallel()
set.seed(1)

# set dimension --------------------------------------------------------------- 

constant <- set_constant()

# set parameter ---------------------------------------------------------------

parameter <- 
  set_parameter(
    constant = constant
  )

# draw exogenous variable -----------------------------------------------------
set_product_characteristics(
  constant = constant
  )

set_market_product(
  constant = constant
)

set_demographic(
  constant = constant
)
exogenous <- 
  set_exogenous(
    constant = constant
  )

# draw shock -------------------------------------------------------------------

shock <- 
  set_shock(
    constant = constant
  )

# set endogenous variable ------------------------------------------------------
delta <- 
  compute_delta(
    M = exogenous$M,
    X = exogenous$X,
    shock = shock,
    parameter = parameter
  )

mu <- 
  compute_mu(
    M = exogenous$M,
    X = exogenous$X,
    D = exogenous$D,
    parameter = parameter
  )

u_df <- 
  compute_indirect_utility(
    delta = delta,
    mu = mu
  )


u_df <- 
  classify_group(u_df)


u_df <-
  adjust_utility(
    u_df = u_df,
    rho = parameter$rho
    )

u_df <-
  compute_group_inclusive_value(
    u_df = u_df,
    rho = parameter$rho
  )


u_df <-
  compute_market_inclusive_value(
    u_df = u_df
    )


u_df <- 
  compute_share_df(
    u_df = u_df,
    rho = parameter$rho
  )

share_df <-
  compute_market_share_df(
    u_df = u_df
  )


cost <-
  compute_marginal_cost(
    M = exogenous$M,
    shock = shock,
    parameter = parameter
  )


omega <- 
  set_ownership(
    M = exogenous$M,
    T_pre = constant$T_pre,
    kappa = parameter$kappa
  )


jacobian <-
  compute_jacobian_all(
    u_df = u_df,
    D = exogenous$D,
    constant = constant,
    parameter = parameter
  )


p <- cost[cost$j>0, "p"]
logp <- log(rep(1, dim(p)[1]))

p_new <- 
  update_price(
    logp = logp,
    M = exogenous$M,
    X = exogenous$X,
    D = exogenous$D,
    parameter = parameter,
    constant = constant,
    shock = shock,
    omega = omega
  )


# set equilibrium --------------------------------------------------------------
distance <- 10000
lambda = 1e-6
while (distance > lambda) {
  p_old <- p_new
  p_new <- 
    update_price(
      logp = log(p_old), 
      M = exogenous$M,
      X = exogenous$X,
      D = exogenous$D,
      parameter = parameter,
      constant = constant,
      shock = shock,
      omega = omega
    )
  distance <- max(abs(p_new - p_old))
  print(distance)
}
p_equilibrium <- p_new

if (!dir.exists("output/simulate")) {
    dir.create("output/simulate", 
    recursive = TRUE)
}
saveRDS(
  p_equilibrium, 
  file = "output/simulate/p_equilibrium.rds"
  )
p_equilibrium <-
  readRDS(file = "output/simulate/p_equilibrium.rds")

M <- exogenous$M %>%
  dplyr::arrange(r, t, j)  # ensure proper ordering if needed

# Update prices for inside goods only
M[M$j > 0, "p"] <- p_equilibrium

delta <- compute_delta(
  M = M,
  X = exogenous$X,
  shock = shock,
  parameter = parameter
)

mu <- compute_mu(
  M = M,
  X = exogenous$X,
  D = exogenous$D,
  parameter = parameter
)

u_df <- compute_indirect_utility(
  delta = delta,
  mu = mu
)

u_df <- adjust_utility(
  u_df = u_df,
  rho = parameter$rho
)

u_df <- compute_group_inclusive_value(
  u_df = u_df,
  rho = parameter$rho
)

u_df <- compute_market_inclusive_value(
  u_df = u_df
)

u_df <- compute_share_df(
  u_df = u_df,
  rho = parameter$rho
)

s_df <- compute_market_share_df(
  u_df = u_df
)

summary_df <- 
  M %>%
  dplyr::filter(j > 0) %>%
  dplyr::select(j, r, t, p) %>%
  dplyr::left_join(
    s_df, 
    by = c("j", "r", "t")
  ) %>%
  dplyr::arrange(r, t, j)

summary_df <- summary_df %>%
  dplyr::mutate(
    merger_status = dplyr::if_else(
      t <= constant$T_pre, 
      "pre", 
      "post"
    )
  )

summary_table <- 
  summary_df %>%
  dplyr::mutate(
    merger_status = dplyr::if_else(
      t <= constant$T_pre,
      "pre",
      "post"
    )
  ) %>%
  dplyr::group_by(j, merger_status) %>%
  dplyr::summarise(
    mean_price = mean(p, na.rm = TRUE),
    sd_price   = sd(p, na.rm = TRUE),
    mean_share = mean(s, na.rm = TRUE),
    sd_share   = sd(s, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(j, merger_status)
summary_table