# initialize -------------------------------------------------------------------

rm(list = ls())
devtools::load_all(".")
library(foreach)
library(magrittr)
library(codetools)
set.seed(1)

# set dimension --------------------------------------------------------------- 

constant <- set_constant()
findGlobals(
  set_constant
)

# set parameter ---------------------------------------------------------------

parameter <- set_parameter()
findGlobals(
  set_parameter
)

# draw exogenous variable -----------------------------------------------------
exogenous <- 
  set_exogenous(
    constant = constant
  )
findGlobals(
  set_exogenous
)

# draw shock -------------------------------------------------------------------

shock <- 
  set_shock(
    constant = constant
  )
findGlobals(
  set_shock
)

# set endogenous variable ------------------------------------------------------
x_rt <- exogenous$x[[1]][[1]]
initial_price <- 
  foreach (
    t = seq_len(
      constant$num_period_before + 
      constant$num_period_after
    )
  ) %do% {
    price_t <- 
      foreach (
        r = seq_len(
          constant$num_market
        )
      ) %do% {
        matrix(
          0,
          nrow = constant$num_firm,
          ncol = 1
        )
      }
    return(price_t)
  }

p_rt <- initial_price[[1]][[1]]
xi_rt <- shock$demand$xi[[1]][[1]]
sigma_d <- shock$demand$sigma_d   # already J x 1
tau_d_t <- shock$demand$tau_d[1]
d_rt <- exogenous$d[[1]][[1]] 
rho <- 0.8

delta_rt <- 
    compute_delta_rt(
        x_rt = x_rt,
        p_rt = p_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d = tau_d_t,
        parameter = parameter
    )
delta_rt

mu_irt <- 
    compute_mu_irt(
        x_rt = x_rt,
        p_rt = p_rt,
        d_rt = d_rt,
        num_consumer = constant$num_consumer,
        parameter = parameter
    )
mu_irt

u_irt <- 
    compute_u_irt(
        delta_rt = delta_rt,
        mu_irt = mu_irt,
        num_consumer = constant$num_consumer
    )
u_irt

u_irt/ (1 - parameter$rho) 
I_i1rt <- compute_inclusive_value_i1rt(
    u_irt = u_irt,
    rho = rho
)
I_i1rt

I_irt <- 
    compute_inclusive_value_irt(
        inclusive_value_i1rt = I_i1rt,
        rho = rho
    )
I_irt

s_ijrt <- 
    compute_share_ijrt(
        u_irt = u_irt,
        inclusive_value_i1rt = I_i1rt,
        inclusive_value_irt = I_irt,
        rho = rho
    )
s_ijrt

s_jrt <- 
    compute_market_share_jrt(
        share_ijrt = s_ijrt
    )
s_jrt

s_ijrt_wrapped <- 
    compute_share_ijrt_wrapper(
        x_rt = x_rt,
        p_rt = p_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d_t = tau_d_t,
        d_rt = d_rt,
        num_consumer = constant$num_consumer,
        parameter = parameter
    )
s_ijrt_wrapped

s_jrt_wrapped <- 
    compute_market_share_jrt(
        share_ijrt = s_ijrt_wrapped
    )
s_jrt_wrapped

jacobian_rt <- 
    compute_jacobian_rt(
        s_ijrt = s_ijrt_wrapped,
        d_rt = d_rt,
        parameter = parameter,
        constant = constant
    )
jacobian_rt

w_rt <- exogenous$w[[1]][[1]]
sigma_s <- shock$cost$sigma_s
tau_s_t <- shock$cost$tau_s[1]
mu_s_r <- shock$cost$mu_s[1]
eta_rt <- shock$cost$eta[[1]][[1]]
eta_rt

mc_rt <- 
    compute_marginal_cost_rt(
        w_rt = w_rt,
        sigma_s = sigma_s,
        tau_s_t = tau_s_t,
        mu_s_r = mu_s_r,
        eta_rt = eta_rt,
        parameter = parameter
    )
mc_rt

owner_rt <- exogenous$owner[[1]][[1]]
owner_rt

p_rt_new <- 
    update_price_rt(
        x_rt = x_rt,
        w_rt = w_rt,
        d_rt = d_rt,
        owner_rt = owner_rt,
        sigma_s = sigma_s,
        tau_s_t = tau_s_t,
        mu_s_r = mu_s_r,
        eta_rt = eta_rt,
        s_ijrt = s_ijrt_wrapped,
        parameter = parameter,
        constant = constant 
    )
p_rt_new

initial_endogenous <- 
    initialize_endogenous_price(
        constant = constant
    )
initial_endogenous$price[[1]][[1]]

endogenous <- 
    set_endogenous(
        constant = constant,
        parameter = parameter,
        exogenous = exogenous,
        shock = shock,
        endogenous = initial_endogenous
    )
endogenous$share[[1]][[1]]
endogenous$price[[1]][[1]]


findGlobals(set_endogenous)
# set equilibrium --------------------------------------------------------------
equilibrium <- 
    set_equilibrium(
        constant = constant
    )

if (!dir.exists("output/simulate")) {
    dir.create("output/simulate", 
    recursive = TRUE)
}
saveRDS(
  equilibrium, 
  file = "output/simulate/equilibrium.rds"
  )
equilibrium <-
  readRDS(file = "output/simulate/equilibrium.rds")

equilibrium$endogenous$share[[1]][[1]]
equilibrium$endogenous$price[[1]][[1]]

findGlobals(set_equilibrium)