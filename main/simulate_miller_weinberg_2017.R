# initialize -------------------------------------------------------------------

rm(list = ls())
devtools::load_all(".")
library(foreach)
library(magrittr)
library(codetools)
set.seed(1)

# set dimension --------------------------------------------------------------- 

constant <- 
  set_constant()


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

# Summary on Equilibrium ---------------------------------------------------
equilibrium$endogenous$price[[1]][[1]]


equilibrium_summary <- list()

T_pre <- 
  seq_len(
    constant$num_period_before
  )

T_post <- 
  seq(
    constant$num_period_before + 1,
    constant$num_period_before + constant$num_period_after
  )

equilibrium_summary$average_pre_merge_price <- 
  average_firm_metric(
    equilibrium$endogenous$price,
    T_pre
  )

equilibrium_summary$average_post_merge_price <- 
  average_firm_metric(
    equilibrium$endogenous$price,
    T_post
  )

equilibrium_summary$average_pre_merge_share <- 
  average_firm_metric(
    equilibrium$endogenous$share,
    T_pre
  )

equilibrium_summary$average_post_merge_share <- 
  average_firm_metric(
    equilibrium$endogenous$share,
    T_post
  )

summary_table <- 
  data.frame(
    Firm = seq_len(constant$num_firm),
    Avg_Price_Pre = equilibrium_summary$average_pre_merge_price,
    Avg_Price_Post = equilibrium_summary$average_post_merge_price,
    Price_Change_Pct = 
      100 * 
      (
        equilibrium_summary$average_post_merge_price - 
        equilibrium_summary$average_pre_merge_price
      ) / equilibrium_summary$average_pre_merge_price,
    Avg_Share_Pre = equilibrium_summary$average_pre_merge_share,
    Avg_Share_Post = equilibrium_summary$average_post_merge_share,
    Share_Change_Pct = 
      100 * 
      (
        equilibrium_summary$average_post_merge_share - 
        equilibrium_summary$average_pre_merge_share
      ) / equilibrium_summary$average_pre_merge_share
  )

summary_table <- 
  round(
    summary_table,
    4
  )

print(summary_table)