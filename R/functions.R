set_constant <- function() {
  list(
    T_pre = 10,         # Pre-merger periods
    T_post = 10,        # Post-merger periods
    R = 10,             # Number of regional markets
    J = 4,             # Total number of potential products
    N = 500,            # Number of consumers per market
    Kx = 3,             # Number of demand-side product characteristics
    Kw = 2              # Number of cost-side product characteristics
  )
}

set_parameter <- function(constant) {
  with(constant, {
    alpha <- -0.1

    beta <- 
      rnorm(
        Kx,
        mean = 0,
        sd = 1
      )
    beta[1] <- 1

    Pi <- 
      rnorm(
        Kx + 1,
        mean = 0,
        sd = 0.001
      )
    Pi[1] <- 0.001

    gamma <- 
        abs(rnorm(
          Kw,
          mean = 0,
          sd = 0.3
        ))
    gamma[1] <- -gamma[1]

    kappa <- 0.5

    rho <- 0.6

    list(
      alpha = alpha,
      beta  = beta,
      Pi    = Pi,
      gamma = gamma,
      kappa = kappa,
      rho = rho
    )
  })
}

set_product_characteristics <- function(constant) {
  with(constant, {
    
    X <- 
      matrix(
        2 * rnorm(
          J * (Kx - 1)
        ),
        nrow = J
      )

    X <- 
      cbind(
        rep(1, J),   # intercept
        X
      )

    colnames(X) <- 
      paste(
        "x", 
        1:Kx, 
        sep = "_"
      )

    X <- 
      data.frame(
        j = 1:J,
        X
      ) %>%
      tibble::as_tibble()

    X <- 
      rbind(
        rep(0, dim(X)[2]),
        X
      )

    return(X)
  })
}

set_market_product <- function(constant) {
  with(constant, {

    # Index set
    M <- 
      expand.grid(
        j = 1:J,
        r = 1:R,
        t = 1:(T_pre + T_post)
      ) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(j, r, t)

    # Cost-side covariates (w_1 to w_Kw)
    for (k in 1:Kw) {
        if (k == 1) {
            # For w_1: indicator variable, equals 1 if j is 2 or 3 and t > 10, else 0
            M[[paste0("w_", k)]] <- as.numeric((M$j %in% c(2, 3)) & (M$t > 10))
        } else {
            M[[paste0("w_", k)]] <- 0.01 * abs(rnorm(nrow(M)))
        }
    }

    # Draw fixed effects
    sigma_s <- 
        abs(rnorm(J, mean = 0, sd = 0.1))
    tau_s <- 
        abs(rnorm(T_pre + T_post, mean = 0, sd = 0.1))
    mu_s <- 
        abs(rnorm(R, mean = 0, sd = 0.1))
    sigma_d <- rnorm(J) * 0.1
    tau_d   <- rnorm(T_pre + T_post) * 0.1

    # Assign fixed effects to rows
    M$sigma_s <- sigma_s[M$j]
    M$mu_s    <- mu_s[M$r]
    M$tau_s   <- tau_s[M$t]
    M$sigma_d <- sigma_d[M$j]
    M$tau_d   <- tau_d[M$t]

    # Initialize prices
    M$p <- 0

    # Add outside option (j = 0)
    outside <-
      expand.grid(
        j = 0,
        r = 1:R,
        t = 1:(T_pre + T_post)
      ) %>%
      tibble::as_tibble()

    for (k in 1:Kw) {
      outside[[paste0("w_", k)]] <- 0
    }

    outside <- 
      outside %>%
      dplyr::mutate(
        sigma_s = 0,
        mu_s    = 0,
        tau_s   = 0,
        sigma_d = 0,
        tau_d   = 0,
        p       = 0
      )

    M <- 
      dplyr::bind_rows(
        M,
        outside
      ) %>%
      dplyr::arrange(
        t, r, j
      )

    return(M)
  })
}

set_demographic <- function(
  constant
) {
  with(constant, {

    D <- 
      expand.grid(
        i = 1:N,
        r = 1:R,
        t = 1:(T_pre + T_post)
      ) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(
        r, 
        t, 
        i
      ) %>%
      dplyr::group_by(
        r, t
      ) %>%
      dplyr::mutate(
        D = rnorm(
          n = 1
        )
      ) %>%
      dplyr::ungroup()

    return(D)
  })
}

set_exogenous <- function(
  constant
) {
  X <- set_product_characteristics(
    constant = constant
  )

  M <- set_market_product(
    constant = constant
  )

  D <- set_demographic(
    constant = constant
  )

  return(
    list(
      X = X,
      M = M,
      D = D
    )
  )
}

set_shock <- function(
  constant
) {
  with(constant, {

    # Main shocks (excluding outside option)
    shock <- 
      expand.grid(
        j = 1:J,
        r = 1:R,
        t = 1:(T_pre + T_post)
      ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        xi = 0.2 * rnorm(
          n = nrow(.)
        ),
        eta = 0.2 * abs(
            rnorm(
            n = nrow(.)
            )
        )    
      )

    # Outside option (j = 0) with xi = eta = 0
    outside <- 
      expand.grid(
        j = 0,
        r = 1:R,
        t = 1:(T_pre + T_post)
      ) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        xi = 0,
        eta = 0
      )

    # Combine and return
    dplyr::bind_rows(
      shock,
      outside
    ) %>%
    dplyr::arrange(
      t, 
      r, 
      j
    )
  })
}

attach_X <- function(
    M, 
    X
    ) {
  dplyr::left_join(
    M, 
    X, 
    by = "j"
    )
}

attach_shock <- function(
    M, 
    shock
    ) {
  dplyr::left_join(
    M, 
    shock, 
    by = c("j", "r", "t")
    )
}

compute_delta <- function(
  M, 
  X, 
  shock, 
  parameter
) {

  M_augmented <- 
    M %>%
    attach_X(
      X = X
    ) %>%
    attach_shock(
      shock = shock
    )

  M_augmented %>%
    dplyr::mutate(
      xb = as.numeric(
        as.matrix(
          dplyr::select(
            ., 
            dplyr::starts_with("x_")
          )
        ) %*% parameter$beta
      ),
      delta = xb + 
              parameter$alpha * p + 
              sigma_d + 
              tau_d + 
              xi
    ) %>%
    dplyr::select(
      j, 
      r, 
      t, 
      delta
    )
}

attach_D <- function(
  M, 
  D
) {
  dplyr::left_join(
    D, 
    M, 
    by = c("r", "t")
  ) %>%
  dplyr::select(
    i, 
    j, 
    r, 
    t, 
    everything()
  )
}

compute_mu <- function(
  M, 
  X, 
  D, 
  parameter
) {

  M_augmented <- 
    attach_D(
      M = M,
      D = D
    ) %>%
    attach_X(
      X = X
    )

  M_augmented %>%
    dplyr::mutate(
      mu = as.numeric(
        as.matrix(
          dplyr::bind_cols(
            dplyr::select(
              ., 
              p
            ),
            dplyr::select(
              ., 
              dplyr::starts_with("x_")
            )
          )
        ) %*% parameter$Pi * D
      )
    ) %>%
    dplyr::select(
      i, 
      j, 
      r, 
      t, 
      mu
    ) %>%
    dplyr::arrange(
      r, 
      t, 
      i 
)
}

compute_indirect_utility <- function(
  delta, 
  mu
) {
  dplyr::left_join(
    mu, 
    delta, 
    by = c("j", "r", "t")
  ) %>%
  dplyr::mutate(
    u = delta + mu
  ) %>%
  dplyr::select(
    i, 
    j, 
    r, 
    t, 
    u
  ) %>%
  dplyr::arrange(
    r, 
    t, 
    i,
    j 
  )
}

classify_group <- function(
  u_df
) {
  u_df %>%
    dplyr::mutate(
      g = ifelse(
        j == 0, 
        0, 
        1
      )
    ) %>%
    dplyr::select(
      i, 
      j, 
      r, 
      t, 
      g, 
      u
    ) %>%
    dplyr::arrange(
      r,
      t,
      i,
      j
    )
}

adjust_utility <- function(
  u_df, 
  rho
) {
  u_df %>%
    dplyr::mutate(
      g = ifelse(j == 0, 0, 1),      # group: 0 = outside good, 1 = inside goods
      u_adj = u / (1 - rho)         # adjusted utility
    ) %>%
    dplyr::arrange(
      r,
      t,
      i,
      j
    )
}

compute_group_inclusive_value <- function(
  u_df, 
  rho
) {
  I_g_df <- 
    u_df %>%
    dplyr::filter(g == 1) %>%
    dplyr::group_by(
      i, 
      r, 
      t, 
      g
    ) %>%
    dplyr::summarise(
      I_g = (1 - rho) * log(sum(exp(u_adj))),
      .groups = "drop"
    )

  u_df %>%
    dplyr::left_join(
      I_g_df,
      by = c("i", "r", "t", "g")
    ) %>%
    dplyr::mutate(
      I_g = ifelse(g == 0, 0, I_g)
    ) %>%
    dplyr::arrange(
      r,
      t,
      i,
      j
    )
}

compute_market_inclusive_value <- function(
  u_df
) {
  I_df <- 
    u_df %>%
    dplyr::filter(g == 1) %>%
    dplyr::distinct(i, r, t, .keep_all = TRUE) %>%
    dplyr::mutate(
      I = log(1 + exp(I_g))
    ) %>%
    dplyr::select(i, r, t, I)

  u_df %>%
    dplyr::left_join(
      I_df,
      by = c("i", "r", "t")
    ) %>%
    dplyr::arrange(
      r,
      t,
      i,
      j
    )
}

compute_share_df <- function(
  u_df,
  rho
) {
  u_df %>%
    dplyr::mutate(
      s_i = exp((u - I_g) / (1 - rho) + I_g - I)
    ) %>%
    dplyr::arrange(
        i,
        j,
        r, 
        t
        )
}

compute_market_share_df <- function(u_df) {
  u_df %>%
    dplyr::group_by(
        j, 
        r, 
        t
        ) %>%
    dplyr::summarise(
      s = mean(s_i),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
        r, 
        t, 
        j
        )
}

compute_marginal_cost <- function(
  M, 
  shock, 
  parameter
) {
  M %>%
    dplyr::left_join(shock, by = c("j", "r", "t")) %>%
    dplyr::mutate(
      mc = as.numeric(
        as.matrix(dplyr::select(., dplyr::starts_with("w_"))) %*% parameter$gamma
      ) + sigma_s + tau_s + mu_s + eta
    ) %>%
    dplyr::arrange(
      r,
      t,
      j
    )
}

set_ownership <- function(
  M,
  T_pre,
  kappa
) {
  omega <- list()

  for (r_k in unique(M$r)) {
    omega[[r_k]] <- list()

    for (t_k in unique(M$t)) {
      M_rt <- 
        M %>%
        dplyr::filter(
          r == r_k,
          t == t_k,
          j > 0
        ) %>%
        dplyr::distinct(j, .keep_all = TRUE) %>%
        dplyr::arrange(j)

      firm <- M_rt$j

      if (length(firm) == 0) {
        omega[[r_k]][[t_k]] <- NULL
        next
      }

      if (t_k > T_pre) {
        firm <- ifelse(firm %in% c(2, 3), 23, firm)
      }

      omega_rt <- 
        outer(
          firm,
          firm,
          FUN = function(f1, f2) as.numeric(f1 == f2)
        )

      if (t_k > T_pre) {
        for (i in 1:length(firm)) {
          for (j in 1:length(firm)) {
            if ((firm[i] == 1 && firm[j] == 23) ||
                (firm[i] == 23 && firm[j] == 1)) {
              omega_rt[i, j] <- kappa
            }
          }
        }
      }

      omega[[r_k]][[t_k]] <- omega_rt
    }
  }

  return(omega)
}

compute_jacobian_rt <- function(
  u_df,       # enriched u_df with s_i
  D,          # demographic df with i, r, t, D
  r,
  t,
  parameter
) {
  # Filter to inside goods only (j != 0)
  u_rt <- 
    u_df %>%
    dplyr::filter(r == !!r, t == !!t, j != 0)

  D_rt <- 
    D %>%
    dplyr::filter(r == !!r, t == !!t) %>%
    dplyr::arrange(i)

  alpha_i <- parameter$alpha + parameter$Pi[1] * D_rt$D
  rho     <- parameter$rho
  N       <- length(unique(u_rt$i))
  J       <- length(unique(u_rt$j))  # only inside goods

  # Wide matrix: rows = i, cols = j (j != 0)
  s_mat <- 
    u_rt %>%
    dplyr::distinct(i, j, .keep_all = TRUE) %>%
    dplyr::select(i, j, s_i) %>%
    tidyr::pivot_wider(
      names_from = j,
      values_from = s_i
    ) %>%
    dplyr::select(-i) %>%
    as.matrix()

  # Group membership of products (ordered by j)
  group_vec <- 
    u_rt %>%
    dplyr::distinct(j, g) %>%
    dplyr::arrange(j) %>%
    dplyr::pull(g)

  # Matrix to hold Jacobian
  jacobian <- matrix(0, nrow = J, ncol = J)

  # Efficient loop with vectorization inside
  for (j in 1:J) {
    for (k in 1:J) {
      g_same <- group_vec[j] == group_vec[k]
      s_j <- s_mat[, j]
      s_k <- s_mat[, k]

      if (g_same) {
        if (j == k) {
          d_ijk <- alpha_i * (1 / (1 - rho)) * s_j * (1 - s_j)
        } else {
          d_ijk <- - alpha_i * (1 / (1 - rho)) * s_j * s_k
        }
      } else {
        d_ijk <- -alpha_i * s_j * s_k
      }

      jacobian[j, k] <- mean(d_ijk)
    }
  }

  return(jacobian)
}

compute_jacobian_all <- function(
  u_df,
  D,
  constant,
  parameter
) {
  T <- constant$T_pre + constant$T_post
  foreach(r = 1:constant$R) %do% {
    foreach(t = 1:T) %do% {
      compute_jacobian_rt(
        u_df = u_df,
        D = D,
        r = r,
        t = t,
        parameter = parameter
      )
    }
  }
}

update_price <- function(
  logp,
  M,
  X,
  D,
  parameter,
  constant,
  shock,
  omega
) {
  p <- exp(logp)
  M[M$j>0, "p"] <- p

  cost <- 
    compute_marginal_cost(
        M = M,
        shock = shock,
        parameter = parameter
    )
    
  # 3. Continue with the rest of the update process
  # Compute delta and mu to get indirect utility u
  delta <- compute_delta(
    M = M,
    X = X,
    shock = shock,
    parameter = parameter
  )
  
  mu <- compute_mu(
    M = M,
    X = X,
    D = D,
    parameter = parameter
  )
  
  u_df <- compute_indirect_utility(
    delta = delta,
    mu = mu
  )

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

  
  # 4. Compute individual-level shares (s_i)
  u_df <- compute_share_df(
    u_df = u_df,
    rho = parameter$rho
  )

  s_df <- compute_market_share_df(
    u_df = u_df
  )
  
  # 5. Compute Jacobian for market (r, t)
  jacobian <- compute_jacobian_all(
    u_df = u_df,
    D = D,
    constant = constant,
    parameter = parameter
  )
  
  p_new <-
    foreach(
        rr = 1:constant$R,
        .packages = c(
          "foreach", 
          "magrittr", 
          "Economics"
        ),
        .combine = "rbind"
    ) %dopar% {
        foreach(
            tt = 1:(constant$T_pre + constant$T_post),
            .packages = c(
            "foreach", 
            "magrittr", 
            "Economics"
            ),
            .combine = "rbind"
        ) %dopar% {
            omega_rt <- omega[[rr]][[tt]]
            s_rt_df <- s_df %>%
                dplyr::filter(r == rr, t == tt, j > 0)
            s_rt  <- as.matrix(s_rt_df$s)
            cost_rt <- cost %>% 
                dplyr::filter(r == rr, t == tt, j > 0)
            mc_rt <- as.matrix(cost_rt$mc)
            jacobian_rt <- jacobian[[rr]][[tt]]
            foc_matrix <- omega_rt * t(jacobian_rt)
            markup_rt <- 
                - solve(
                    foc_matrix,
                    s_rt
                )
            p_new_rt <- mc_rt + markup_rt
            return(p_new_rt)
        }
    }
  return(p_new)
}

summarize