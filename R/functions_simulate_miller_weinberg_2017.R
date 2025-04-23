set_constant <-
  function(

  ) {
    num_period_before <- 10
    num_period_after <- 10
    num_market <- 10
    num_firm <- 4
    num_consumer <- 1000
    merge <-
      c(
        2,
        3
      )

    return(
      list(
        num_period_before = num_period_before,
        num_period_after = num_period_after,
        num_market = num_market,
        num_firm = num_firm,
        num_consumer = num_consumer,
        merge = merge
      )
    )
  }
  
set_parameter <-
  function(

  ) {
    alpha <- -0.1
    beta <- 1
    rho <- 0.1
    pi_alpha <- 0.01
    pi_beta <- 0.01
    gamma <- 1
    kappa <- 1

    demand <-
      list(
        alpha = alpha,
        beta = beta,
        rho = rho,
        pi_alpha = pi_alpha,
        pi_beta = pi_beta
      )

    cost <-
      list(
        gamma = gamma
      )

    conduct <-
      list(
        kappa = kappa
      )

    parameter <-
      list(
        demand = demand,
        cost = cost,
        conduct = conduct
      )

    return(parameter)
  }

set_exogenous <-
  function(
    constant
  ) {
    x <- 
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        x_t <- 
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            x_rt <-
              rnorm(
                constant$num_firm
              ) %>%
              as.matrix()
            return(x_rt)
          }
        return(x_t)
      }

    w <-
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        w_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            w_rt <-
              rnorm(
                constant$num_firm
              ) %>%
              as.matrix()
            return(w_rt)
          }
        return(w_t)
      }

    d <- 
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        d_t <- 
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            d_rt <-
              rnorm(
                1
              ) %>%
              as.matrix()
            return(d_rt)
          }
        return(d_t)
      }

    owner_before <-
      foreach (
        t = seq_len(
          constant$num_period_before
        )
      ) %do% {
        owner_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            owner_rt <-
              rep(
                1,
                constant$num_firm
              ) %>%
              diag()
            return(owner_rt)
          }
        return(owner_t)
      }

    owner_after <-
      foreach (
        t = seq_len(
          constant$num_period_after
        )
      ) %do% {
        owner_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            owner_rt <-
              rep(
                1,
                constant$num_firm
              ) %>%
              diag()
            owner_rt[
              constant$merge,
              constant$merge
            ] <- 1
            return(owner_rt)
          }
        return(owner_t)
      }

    owner <-
      c(
        owner_before,
        owner_after
      )

    return(
      list(
        x = x,
        w = w,
        d = d,
        owner = owner
      )
    )
  }

set_shock <-
  function(
    constant
  ) {
    sigma_d <-
      rnorm(
        constant$num_firm
      ) %>%
      as.matrix()

    tau_d <-
      rnorm(
        constant$num_period_before + 
        constant$num_period_after
      ) %>%
      as.matrix()

    xi <-
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        xi_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            xi_rt <-
              rnorm(
                constant$num_firm
              ) %>%
              as.matrix()
            return(xi_rt)
          }
        return(xi_t)
      }

    sigma_s <-
      rnorm(
        constant$num_firm
      ) %>%
      as.matrix()

    tau_s <-
      rnorm(
        constant$num_period_before + 
        constant$num_period_after
      ) %>%
      as.matrix()

    mu_s <-
      rnorm(
        constant$num_market
      ) %>%
      as.matrix()

    eta <-
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        eta_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            eta_rt <-
              rnorm(
                constant$num_firm
              ) %>%
              as.matrix()
            return(eta_rt)
          }
        return(eta_t)
      }

    demand <-
      list (
        sigma_d = sigma_d,
        tau_d = tau_d,
        xi = xi
      )

    cost <-
      list(
        sigma_s = sigma_s,
        tau_s = tau_s,
        mu_s = mu_s,
        eta = eta
      )

    return(
      list(
        demand = demand,
        cost = cost
      )
    )
  }

compute_delta_rt <- 
  function(
    x_rt,
    price_rt,
    xi_rt,
    sigma_d,
    tau_d_t,
    alpha,
    beta
  ) {
    delta <- 
      x_rt * beta +
      price_rt * alpha +
      xi_rt +
      sigma_d +
      matrix(
        tau_d_t,
        nrow = nrow(x_rt),
        ncol = 1
      )
    return(
      delta
    )
  }

compute_mu_irt <- 
  function(
    x_rt,
    d_rt,
    price_rt,
    num_consumer,
    pi_alpha,
    pi_beta
  ) {

    d_vec <- 
      matrix(
        d_rt,
        nrow = num_consumer,
        ncol = 1
      )

    x_term <- 
      x_rt %*% t(d_vec) * pi_beta    # J x N

    p_term <- 
      price_rt %*% t(d_vec) * pi_alpha   # J x N

    mu_irt <- 
      x_term + p_term               # J x N

    return(mu_irt)
  }

compute_u_irt <- 
  function(
    delta_rt,
    mu_irt,
    num_consumer
  ) {
    delta_mat <- 
      delta_rt %*% 
      matrix(
        1,
        nrow = 1,
        ncol = num_consumer
      )  # J x N

    u_irt <- 
      delta_mat + 
      mu_irt      # J x N

    return(u_irt)
  }

compute_inclusive_value_i1rt <- 
  function(
    u_irt,
    rho
  ) {
    u_adj <- 
      u_irt / (1 - rho)  # J x 

    inclusive_value_i1rt <- 
      (1 - rho) * 
      log(
        colSums(
          exp(u_adj)
        )
      )  # 1 x N

    return(inclusive_value_i1rt)
  }

compute_inclusive_value_irt <- 
  function(
    inclusive_value_i1rt
  ) {
    inclusive_value_irt <- 
      log(
        1 + 
          exp(
            inclusive_value_i1rt
          )
      )

    return(inclusive_value_irt)  # 1 x N vector
  }

compute_share_irt <- 
  function(
    u_irt,
    inclusive_value_i1rt,
    inclusive_value_irt,
    rho
  ) {
    u_adj <- 
      u_irt / (1 - rho)  # J x N

    within_group <- 
      exp(u_adj) / 
      matrix(
        exp(inclusive_value_i1rt / (1 - rho)),
        nrow = nrow(u_irt),
        ncol = ncol(u_irt),
        byrow = TRUE
      )  # J x N

    across_group <- 
      matrix(
        exp(inclusive_value_i1rt) / exp(inclusive_value_irt),
        nrow = nrow(u_irt),
        ncol = ncol(u_irt),
        byrow = TRUE
      )  # J x N

    s_irt <- within_group * across_group

    return(s_irt)
  }

compute_share_irt_wrapper <- 
  function(
    x_rt,
    price_rt,
    xi_rt,
    sigma_d,
    tau_d_t,
    d_rt,
    num_consumer,
    alpha,
    beta,
    pi_alpha,
    pi_beta,
    rho
  ) {
    delta_rt <- 
      compute_delta_rt(
        x_rt = x_rt,
        price_rt = price_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d_t = tau_d_t,
        alpha = alpha,
        beta = beta
      )

    mu_irt <- 
      compute_mu_irt(
        x_rt = x_rt,
        d_rt = d_rt,
        price_rt = price_rt,
        num_consumer = num_consumer,
        pi_alpha = pi_alpha,
        pi_beta = pi_beta
      )

    u_irt <- 
      compute_u_irt(
        delta_rt = delta_rt,
        mu_irt = mu_irt,
        num_consumer = num_consumer
      )

    if (rho == 1) {
      s_irt <- 
        exp(u_irt) / 
        matrix(
          1 +
          colSums(exp(u_irt)),
          nrow = nrow(u_irt),
          ncol = ncol(u_irt),
          byrow = TRUE
        )
    } else {
      inclusive_i1rt <- 
        compute_inclusive_value_i1rt(
          u_irt = u_irt,
          rho = rho
        )

      inclusive_irt <- 
        compute_inclusive_value_irt(
          inclusive_value_i1rt = inclusive_i1rt
        )

      s_irt <- 
        compute_share_irt(
          u_irt = u_irt,
          inclusive_value_i1rt = inclusive_i1rt,
          inclusive_value_irt = inclusive_irt,
          rho = rho
        )
    }

    return(s_irt)
  }

compute_share_rt <- 
  function(
    share_irt
  ) {
    share_rt <- 
      rowMeans(
        share_irt
      ) %>%
      as.matrix()

    return(share_rt)
  }

compute_share_rt_wrapper <-
  function(
    x_rt,
    price_rt,
    xi_rt,
    sigma_d,
    tau_d_t,
    d_rt,
    num_consumer,
    alpha,
    beta,
    pi_alpha,
    pi_beta,
    rho
  ) {
    s_irt <- 
      compute_share_irt_wrapper(
        x_rt = x_rt,
        d_rt = d_rt,
        price_rt = price_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d_t = tau_d_t,
        num_consumer = num_consumer,
        alpha = alpha,
        beta = beta,
        pi_alpha = pi_alpha,
        pi_beta = pi_beta,
        rho = rho
      )
    share_rt <- 
      compute_share_rt(
        share_irt = s_irt
      )
    return(share_rt)
  }

adjust_owner_rt_with_kappa <- 
  function(
    owner_rt,
    kappa
  ) {
    owner_rt[1, 2] <- kappa
    owner_rt[2, 1] <- kappa
    owner_rt[1, 3] <- kappa
    owner_rt[3, 1] <- kappa
    return(owner_rt)
  }

compute_jacobian_rt <- 
  function(
    s_irt,
    d_rt,
    alpha,
    pi_alpha,
    rho
  ) {
  J <- nrow(s_irt)
  N <- ncol(s_irt)

  alpha_i <- 
    as.vector(
      alpha + pi_alpha * d_rt
    )
  s_i1 <- colSums(s_irt)
  s_i0 <- 1 - s_i1
  s_cond <- 
    s_irt / 
    matrix(
      s_i1, 
      nrow = J, 
      ncol = N, 
      byrow = TRUE
    )

  inv_L <- 1 / (1 - rho)

  jac <- 
    matrix(
    0, 
    J, 
    J
    )

  for (j in 1: J) {
    s_ij <- s_irt[j, ]
    for (k in 1: J) {
      deriv_i <-
        alpha_i * s_ij *
        ( ((j == k) - s_cond[k, ]) * inv_L +
          s_i0 * s_cond[k, ])
      jac[j, k] <- mean(deriv_i)
    }
  }
  jac
}

compute_marginal_cost_rt <- 
  function(
    w_rt,
    sigma_s,
    tau_s_t,
    mu_s_r,
    eta_rt,
    gamma
  ) {
    mc_rt <- 
      w_rt * gamma + 
      sigma_s + 
      tau_s_t + 
      mu_s_r +
      eta_rt

    return(mc_rt)
  }

add_markup_rt <-
  function(
    x_rt,
    w_rt,
    d_rt,
    owner_rt,
    sigma_s,
    tau_s_t,
    mu_s_r,
    eta_rt,
    share_rt,
    jacobian_rt,
    alpha,
    pi_alpha,
    rho,
    gamma
  ) {

    mc_rt <-
      compute_marginal_cost_rt(
        w_rt = w_rt,
        sigma_s = sigma_s,
        tau_s_t = tau_s_t,
        mu_s_r = mu_s_r,
        eta_rt = eta_rt,
        gamma = gamma
      )

    markup <-
      - solve(
        owner_rt * t(jacobian_rt),
        share_rt
      )

    price_rt <- mc_rt + markup
    return(price_rt)
  }

set_endogenous <- 
  function(
    constant
  ) {
    share <-
      foreach (
        t = seq_len(
          constant$num_period_before + 
          constant$num_period_after
        )
      ) %do% {
        share_t <-
          foreach (
            r = seq_len(
              constant$num_market
            )
          ) %do% {
            share_rt <-
              matrix(
                1 / (1 + constant$num_firm),
                nrow = constant$num_firm,
                ncol = 1
              )
            return(share_rt)
          }
        return(share_t)
      }

    price <-
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
            price_rt <-
              matrix(
                0,
                nrow = constant$num_firm,
                ncol = 1
              )
            return(price_rt)
          }
        return(price_t)
      }

    return(
      list(
        share = share,
        price = price
      )
    )
  }

update_price_rt <-
  function(
    x_rt,
    w_rt,
    d_rt,
    owner_rt,
    price_rt,
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
    pi_alpha,
    pi_beta,
    rho,
    gamma
  ) {
    s_irt <- 
      compute_share_irt_wrapper(
        x_rt = x_rt,
        d_rt = d_rt,
        price_rt = price_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d_t = tau_d_t,
        num_consumer = num_consumer,
        alpha = alpha,
        beta = beta,
        pi_alpha = pi_alpha,
        pi_beta = pi_beta,
        rho = rho
      )
    jacobian_rt <-
      compute_jacobian_rt(
        s_irt = s_irt,
        d_rt = d_rt,
        alpha = alpha,
        pi_alpha = pi_alpha,
        rho = rho
      )
    share_rt <-
      compute_share_rt(
        share_irt = s_irt
      )
    price_rt <-
      add_markup_rt(
        x_rt = x_rt,
        w_rt = w_rt,
        d_rt = d_rt,
        owner_rt = owner_rt,
        sigma_s = sigma_s,
        tau_s_t = tau_s_t,
        mu_s_r = mu_s_r,
        eta_rt = eta_rt,
        share_rt = share_rt,
        jacobian_rt = jacobian_rt,
        alpha = alpha,
        pi_alpha = pi_alpha,
        rho = rho,
        gamma = gamma
      )
    return(price_rt)
  }

solve_endogenous_rt <-
  function(
    x_rt,
    w_rt,
    d_rt,
    owner_rt,
    price_rt,
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
    pi_alpha,
    pi_beta,
    rho,
    gamma
  ) {
    fn <-
      function(
        price_rt 
      ) {
        p_rt_new <-
          update_price_rt(
            x_rt = x_rt,
            w_rt = w_rt,
            d_rt = d_rt,
            owner_rt = owner_rt,
            price_rt = price_rt,
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
            pi_alpha = pi_alpha,
            pi_beta = pi_beta,
            rho = rho,
            gamma = gamma
          )
        return(
          price_rt - p_rt_new
        )
      }
    solution <-
      nleqslv::nleqslv(
        x = price_rt,
        fn = fn
      )
    price_rt <- 
      solution$x %>%
      as.matrix()
    share_rt <-
      compute_share_rt_wrapper(
        x_rt = x_rt,
        d_rt = d_rt,
        price_rt = price_rt,
        xi_rt = xi_rt,
        sigma_d = sigma_d,
        tau_d_t = tau_d_t,
        num_consumer = num_consumer,
        alpha = alpha,
        beta = beta,
        pi_alpha = pi_alpha,
        pi_beta = pi_beta,
        rho = rho
      )
    return(
      list(
        share_rt = share_rt,
        price_rt = price_rt
      )
    )
  }

set_equilibrium <- 
  function(
    constant,
    lambda = 1e-6,
    max_iter = 10000
  ) {
  
    # set parameter ---------------------------------------------------------------

    parameter <- set_parameter()

    # draw exogenous variable -----------------------------------------------------

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

    endogenous <- 
      set_endogenous(
        constant = constant
      )

    return(
      list(
        constant = constant,
        parameter = parameter,
        exogenous = exogenous,
        shock = shock,
        endogenous = endogenous
      )
    )
  }

solve_equilibrium <- 
  function(
    equilibrium
  ) {

    endogenous <-
      foreach (
        t = seq_along(
          equilibrium$exogenous$x
        )
      ) %do% {
        endogenous_t <-
          foreach (
            r = seq_along(
              equilibrium$exogenous$x[[t]]
            )
          ) %do% {
            if (
              t > equilibrium$constant$num_period_before
            ) {
              owner_rt <- 
                adjust_owner_rt_with_kappa(
                  owner_rt = equilibrium$exogenous$owner[[t]][[r]],
                  kappa = equilibrium$parameter$conduct$kappa
                )
            } else {
              owner_rt <- equilibrium$exogenous$owner[[t]][[r]]
            }
            endogenous_rt <- 
              solve_endogenous_rt(
                x_rt = equilibrium$exogenous$x[[t]][[r]],
                w_rt = equilibrium$exogenous$w[[t]][[r]],
                d_rt = equilibrium$exogenous$d[[t]][[r]],
                owner_rt = owner_rt,
                price_rt = equilibrium$endogenous$p[[t]][[r]],
                xi_rt = equilibrium$shock$demand$xi[[t]][[r]],
                sigma_d = equilibrium$shock$demand$sigma_d,
                tau_d_t = equilibrium$shock$demand$tau_d[[t]],
                sigma_s = equilibrium$shock$cost$sigma_s,
                tau_s_t = equilibrium$shock$cost$tau_s[[t]],
                mu_s_r = equilibrium$shock$cost$mu_s[[r]],
                eta_rt = equilibrium$shock$cost$eta[[t]][[r]],
                num_consumer = equilibrium$constant$num_consumer,
                alpha = equilibrium$parameter$demand$alpha,
                beta = equilibrium$parameter$demand$beta,
                pi_alpha = equilibrium$parameter$demand$pi_alpha,
                pi_beta = equilibrium$parameter$demand$pi_beta,
                rho = equilibrium$parameter$demand$rho,
                gamma = equilibrium$parameter$cost$gamma
              ) 
            return(endogenous_rt)
          }
        return(endogenous_t)
      }

    for (
      t in seq_along(
        equilibrium$exogenous$x
      )
    ) {
      for (
        r in seq_along(
          equilibrium$exogenous$x[[t]]
        )
      ) {
        equilibrium$endogenous$share[[t]][[r]] <-
          endogenous[[t]][[r]]$share_rt
        equilibrium$endogenous$price[[t]][[r]] <-
          endogenous[[t]][[r]]$price_rt
      }
    }
    return(
      equilibrium
    )
  }

average_firm_metric <- 
  function(
    metric_list,
    period_idx
  ) {
    Reduce(
      "+",
      lapply(
        metric_list[period_idx],
        function(period) {
          Reduce("+", period)
        }
      )
    ) / 
    (
      length(period_idx) * 
      length(metric_list[[1]])
    )
  }
