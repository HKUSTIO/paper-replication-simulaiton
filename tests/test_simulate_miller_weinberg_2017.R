# initialize -------------------------------------------------------------------

rm(list = ls())
devtools::load_all(".")
library(foreach)
library(magrittr)
library(codetools)
set.seed(1)

# set constant ----------------------------------------------------------------

t <- 1
r <- 1

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

# draw ezxogenous variable -----------------------------------------------------
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

endogenous <-
  set_endogenous(
    constant = constant
  )
findGlobals(
  set_endogenous
)

# set equilibrium --------------------------------------------------------------

equilibrium <- 
  set_equilibrium(
    constant = constant
  )

# solve endogenous variable ----------------------------------------------------

delta_rt <- 
  compute_delta_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    p_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d = equilibrium$shock$demand$tau_d[[t]],
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    intercept = equilibrium$parameter$demand$intercept
  )
findGlobals(
  compute_delta_rt
)
delta_rt

mu_irt <- 
  compute_mu_irt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    p_rt = equilibrium$endogenous$p[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    num_consumer = constant$num_consumer,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta
  )
findGlobals(
  compute_mu_irt
)
mu_irt

u_irt <- 
  compute_u_irt(
    delta_rt = delta_rt,
    mu_irt = mu_irt,
    num_consumer = constant$num_consumer
  )
findGlobals(
  compute_u_irt
)
u_irt


u_irt / (1 - equilibrium$parameter$demand$rho) 

I_i1rt <- 
  compute_inclusive_value_i1rt(
    u_irt = u_irt,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_inclusive_value_i1rt
)
I_i1rt

I_irt <- 
  compute_inclusive_value_irt(
    inclusive_value_i1rt = I_i1rt,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_inclusive_value_irt
)
I_irt

s_individual_rt <- 
  compute_share_individual_rt(
    u_irt = u_irt,
    inclusive_value_i1rt = I_i1rt,
    inclusive_value_irt = I_irt,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_share_individual_rt
)
s_individual_rt

s_rt <- 
  compute_market_share_rt(
    share_individual_rt = s_individual_rt
  )
findGlobals(
  compute_market_share_rt
)
s_rt

s_individual_rt_wrapped <- 
  compute_share_individual_rt_wrapper(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    p_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    intercept = equilibrium$parameter$demand$intercept,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_share_individual_rt_wrapper
)
s_individual_rt_wrapped

s_rt_wrapped <- 
  compute_market_share_rt(
    share_individual_rt = s_individual_rt_wrapped
  )
findGlobals(
  compute_market_share_rt
)
s_rt_wrapped

jacobian_rt <- 
  compute_jacobian_rt(
    s_individual_rt = s_individual_rt_wrapped,
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    alpha = equilibrium$parameter$demand$alpha,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_jacobian_rt
)
jacobian_rt

mc_rt <- 
  compute_marginal_cost_rt(
    w_rt = equilibrium$exogenous$w[[t]][[r]],
    sigma_s = equilibrium$shock$cost$sigma_s,
    tau_s_t = equilibrium$shock$cost$tau_s[[t]],
    mu_s_r = equilibrium$shock$cost$mu_s[[t]],
    eta_rt = equilibrium$shock$cost$eta[[t]][[r]],
    gamma = equilibrium$parameter$cost$gamma
  )
findGlobals(
  compute_marginal_cost_rt
)
mc_rt

owner_rt_adj <- 
  adjust_owner_rt_with_kappa(
    owner_rt = equilibrium$exogenous$owner[[t]][[r]],
    kappa = equilibrium$parameter$conduct$kappa
  )

p_rt_new <- 
  add_markup_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    w_rt = equilibrium$exogenous$w[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    owner_rt = equilibrium$exogenous$owner[[1]][[1]],
    sigma_s = equilibrium$shock$cost$sigma_s,
    tau_s_t = equilibrium$shock$cost$tau_s[[t]],
    mu_s_r = equilibrium$shock$cost$mu_s[[t]],
    eta_rt = equilibrium$shock$cost$eta[[t]][[r]],
    s_rt = equilibrium$endogenous$share[[t]][[r]],
    jacobian_rt = jacobian_rt,
    alpha = equilibrium$parameter$demand$alpha,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    rho = equilibrium$parameter$demand$rho,
    gamma = equilibrium$parameter$cost$gamma
  )
findGlobals(
  add_markup_rt
)
p_rt_new

p_rt_new <-
  update_price_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    w_rt = equilibrium$exogenous$w[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    owner_rt = equilibrium$exogenous$owner[[t]][[r]],
    p_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    sigma_s = equilibrium$shock$cost$sigma_s,
    tau_s_t = equilibrium$shock$cost$tau_s[[t]],
    mu_s_r = equilibrium$shock$cost$mu_s[[t]],
    eta_rt = equilibrium$shock$cost$eta[[t]][[r]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    intercept = equilibrium$parameter$demand$intercept,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho,
    gamma = equilibrium$parameter$cost$gamma
  )
findGlobals(
  update_price_rt
)
p_rt_new

endogenous <- 
  update_endogenous(
    equilibrium = equilibrium
  )
findGlobals(
  update_endogenous
)
endogenous$share[[1]][[1]]
endogenous$price[[1]][[1]]

findGlobals(set_endogenous)

# set equilibrium --------------------------------------------------------------

solve_endogenous_rt <-
  function(
    x_rt,
    w_rt,
    d_rt,
    owner_rt,
    p_rt,
    xi_rt,
    sigma_d,
    tau_d_t,
    sigma_s,
    tau_s_t,
    mu_s_r,
    eta_rt,
    num_consumer,
    alpha,
    beta,
    intercept,
    pi_alpha,
    pi_beta,
    rho,
    gamma
  ) {
    fn <-
      function(
        p_rt 
      ) {
        p_rt_new <-
          update_price_rt(
            x_rt = x_rt,
            w_rt = w_rt,
            d_rt = d_rt,
            owner_rt = owner_rt,
            p_rt = p_rt,
            xi_rt = xi_rt,
            sigma_d = sigma_d,
            tau_d_t = tau_d_t,
            sigma_s = sigma_s,
            tau_s_t = tau_s_t,
            mu_s_r = mu_s_r,
            eta_rt = eta_rt,
            num_consumer = num_consumer,
            alpha = alpha,
            beta = beta,
            intercept = intercept,
            pi_alpha = pi_alpha,
            pi_beta = pi_beta,
            rho = rho,
            gamma = gamma
          )
        return(
          p_rt - p_rt_new
        )
      }
    solution <-
      nleqslv::nleqslv(
        x = p_rt,
        fn = fn,
        lower = 0,
        method = "L-BFGS-B"
      )
    p_rt <- 
      solution$x %>%
      as.matrix()
    share_individual_rt <- 
      compute_share_individual_rt_wrapper(
        x_rt = x_rt,
        p_rt = p_rt,
        
      )
    return()
  }

equilibrium <- 
  solve_equilibrium(
    equilibrium = equilibrium
  )
findGLobals(
  solve_equilibrium
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