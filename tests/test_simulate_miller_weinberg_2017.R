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
    price_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d = equilibrium$shock$demand$tau_d[[t]],
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta
  )
findGlobals(
  compute_delta_rt
)
delta_rt

mu_irt <- 
  compute_mu_irt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    price_rt = equilibrium$endogenous$p[[t]][[r]],
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

inclusive_i1rt <- 
  compute_inclusive_value_i1rt(
    u_irt = u_irt,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_inclusive_value_i1rt
)
inclusive_i1rt

inclusive_irt <- 
  compute_inclusive_value_irt(
    inclusive_value_i1rt = inclusive_i1rt
  )
findGlobals(
  compute_inclusive_value_irt
)
inclusive_irt

s_irt <- 
  compute_share_irt(
    u_irt = u_irt,
    inclusive_value_i1rt = inclusive_i1rt,
    inclusive_value_irt = inclusive_irt,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_share_irt
)
s_irt

share_rt <- 
  compute_share_rt(
    share_irt = s_irt
  )
findGlobals(
  compute_share_rt
)
share_rt

s_irt_wrapped <- 
  compute_share_irt_wrapper(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    price_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho
  )
findGlobals(
  compute_share_irt_wrapper
)
s_irt_wrapped

share_rt <-
  compute_share_rt_wrapper(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    price_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho
  ) 
findGlobals(
  compute_share_rt_wrapper
)
share_rt

jacobian_rt <- 
  compute_jacobian_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    price_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho
  )
# Define a function to numerically calculate the derivatives
jacobian_rt_numerical <- 
  compute_jacobian_numerical_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    price_rt = equilibrium$endogenous$p[[t]][[r]],
    xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
    sigma_d = equilibrium$shock$demand$sigma_d,
    tau_d_t = equilibrium$shock$demand$tau_d[[t]],
    num_consumer = constant$num_consumer,
    alpha = equilibrium$parameter$demand$alpha,
    beta = equilibrium$parameter$demand$beta,
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho
  )

findGlobals(
  compute_jacobian_rt
)
findGlobals(
  compute_jacobian_numerical_rt
)
jacobian_rt
jacobian_rt_numerical

jac_long <-
  foreach::foreach(
    tt = seq_along(equilibrium$exogenous$x),
    .combine = "rbind"
  ) %do% {

    foreach::foreach(
      rr = seq_along(equilibrium$exogenous$x[[tt]]),
      .combine = "rbind"
    ) %do% {
      jac_ana <-
        compute_jacobian_rt(
          x_rt = equilibrium$exogenous$x[[tt]][[rr]],
          d_rt = equilibrium$exogenous$d[[tt]][[rr]],
          price_rt = equilibrium$endogenous$p[[tt]][[rr]],
          xi_rt = equilibrium$shock$demand$xi[[tt]][[rr]],
          sigma_d = equilibrium$shock$demand$sigma_d,
          tau_d_t = equilibrium$shock$demand$tau_d[[tt]],
          num_consumer = constant$num_consumer,
          alpha = equilibrium$parameter$demand$alpha,
          beta = equilibrium$parameter$demand$beta,
          pi_alpha = equilibrium$parameter$demand$pi_alpha,
          pi_beta = equilibrium$parameter$demand$pi_beta,
          rho = equilibrium$parameter$demand$rho
        )

      jac_num <-
        compute_jacobian_numerical_rt(
          x_rt = equilibrium$exogenous$x[[tt]][[rr]],
          d_rt = equilibrium$exogenous$d[[tt]][[rr]],
          price_rt = equilibrium$endogenous$p[[tt]][[rr]],
          xi_rt = equilibrium$shock$demand$xi[[tt]][[rr]],
          sigma_d = equilibrium$shock$demand$sigma_d,
          tau_d_t = equilibrium$shock$demand$tau_d[[tt]],
          num_consumer = constant$num_consumer,
          alpha = equilibrium$parameter$demand$alpha,
          beta = equilibrium$parameter$demand$beta,
          pi_alpha = equilibrium$parameter$demand$pi_alpha,
          pi_beta = equilibrium$parameter$demand$pi_beta,
          rho = equilibrium$parameter$demand$rho
        )

      data.frame(
        Numerical  = as.vector(jac_num),
        Analytical = as.vector(jac_ana),
        period     = tt,    # kept for reference / colouring
        market     = rr
      )
    }
  }

ggplot2::ggplot(
  jac_long,
  ggplot2::aes(
    x      = Numerical,
    y      = Analytical
  )
) +
  ggplot2::geom_point(alpha = 0.55) +
  ggplot2::geom_abline(
    slope     = 1,
    intercept = 0,
    colour    = "red",
    linewidth = 1
  ) +
  ggplot2::coord_equal() +
  ggplot2::guides(colour = "none") +     
  ggplot2::labs(
    title = "Analytical vs numerical Jacobian (all markets & periods)",
    x     = "Numerical derivative",
    y     = "Analytical derivative"
  ) +
  ggplot2::theme_minimal()

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
    owner_rt = equilibrium$exogenous$owner[[t]][[r]],
    sigma_s = equilibrium$shock$cost$sigma_s,
    tau_s_t = equilibrium$shock$cost$tau_s[[t]],
    mu_s_r = equilibrium$shock$cost$mu_s[[t]],
    eta_rt = equilibrium$shock$cost$eta[[t]][[r]],
    share_rt = equilibrium$endogenous$share[[t]][[r]],
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
    price_rt = equilibrium$endogenous$p[[t]][[r]],
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
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho,
    gamma = equilibrium$parameter$cost$gamma
  )
findGlobals(
  update_price_rt
)
p_rt_new

# set equilibrium --------------------------------------------------------------

endogenous_rt <- 
  solve_endogenous_rt(
    x_rt = equilibrium$exogenous$x[[t]][[r]],
    w_rt = equilibrium$exogenous$w[[t]][[r]],
    d_rt = equilibrium$exogenous$d[[t]][[r]],
    owner_rt = equilibrium$exogenous$owner[[t]][[r]],
    price_rt = equilibrium$endogenous$price[[t]][[r]],
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
    pi_alpha = equilibrium$parameter$demand$pi_alpha,
    pi_beta = equilibrium$parameter$demand$pi_beta,
    rho = equilibrium$parameter$demand$rho,
    gamma = equilibrium$parameter$cost$gamma
  ) 
findGlobals(
  solve_endogenous_rt
)
endogenous_rt$share[[t]][[r]]
endogenous_rt$price[[t]][[r]]

equilibrium <- 
  solve_equilibrium(
    equilibrium = equilibrium
  )
findGlobals(
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
  readRDS(
    file = "output/simulate/equilibrium.rds"
  )

equilibrium$endogenous$share[[1]][[1]]
equilibrium$endogenous$price[[1]][[1]]

findGlobals(set_equilibrium)