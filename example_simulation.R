library(TeenSTD)
library(doParallel)
library(foreach)
library(data.table)
library(ggplot2)

#### Simulation wrapper: simulate the disease outcomes using a set of parameters
sim_wrapper <- function(par_i) {
  ## create a list of parameters used in the disease model
  tmp_ls <- create_parameters(burnin = 2,
                              treat_fact_f = par_i$treat_fact_f,
                              treat_fact_m = par_i$treat_fact_m,
                              screen_fact_m = par_i$screen_fact_m,
                              beta = par_i$beta,
                              factor = par_i$factor,
                              alpha_m = par_i$alpha_m,
                              incubate = par_i$incubate,
                              fr_sex = par_i$fr_sex,
                              p_condom = par_i$p_condom,
                              p_condom_l = par_i$p_condom_l,
                              r_preg = par_i$r_preg,
                              red_condom = par_i$red_condom,
                              oral_h10 = par_i$oral_h10,
                              ratio_oral_mean_over_h10 = par_i$mean_oral / par_i$oral_h10,
                              IUD_h10 = par_i$IUD_h10,
                              ratio_IUD_mean_over_h10 = par_i$mean_IUD / par_i$IUD_h10,
                              IUD_factor = par_i$IUD_factor,
                              I_m = par_i$I_m,
                              I_f = par_i$I_f,
                              T_m = par_i$T_m,
                              T_f = par_i$T_f)

  ## Get the list of parameters
  par <- tmp_ls$par
  ## Get the vector of initial states
  initials <- tmp_ls$initials

  ## Use ODE solver to get the outcomes at each month
  chl <- ode(times = seq(0, par$max_time, 1),
             y = initials,
             parms = par,
             func = transmission_model)

  ## Converting monthly results to annual results
  output <- process_output(chl, par)

  output$sol <- NULL
  output$flow <- NULL

  return(output)
}

#### Obtain 1000 calibrated posterior parameter sets
parmdf <- posterior_set

colnames(parmdf) <- c("treat_fact_f", "treat_fact_m", "screen_fact_m", "beta",
                      "factor", "alpha_m", "incubate", "fr_sex", "p_condom",
                      "p_condom_l", "r_preg", "red_condom", "mean_oral", "oral_h10",
                      "mean_IUD", "IUD_h10", "IUD_factor", "I_m", "I_f", "T_m", "T_f", "ll")

## Get the number of simulations
nsim <- nrow(parmdf)

## Obtain the target data
target <- c(MDH_data$obsM1[1:9], MDH_data$obsF1[1:9], preg_data$obs[1:9])
target_df <- data.table(time = rep(c(2005:2013), 3),
                        y = target,
                        type = rep(c("chlamydia: male", "chlamydia: female", "pregnancy"), each = 9))
setkeyv(target_df, c("type", "time"))


#### Run parallelized simulations across all 1000 parameter sets
a <- Sys.time()

ncore <- parallel::detectCores() - 3
c1 <- makeCluster(ncore)
registerDoParallel(c1)

outlist <- foreach(i = c(1:nsim),
                   .packages = c("TeenSTD", "deSolve", "data.table")) %dopar% {
                     par_i <- parmdf[i, ]
                     res <- sim_wrapper(par_i)

                     sim_v <- c(res$fluxM, res$fluxF, res$fluxP)
                     llk <- calculate_likelihood(sim_v, target)

                     flow_df <- data.table(time = rep(c(2005:2013), 3),
                                           y = sim_v,
                                           type = rep(c("chlamydia: male", "chlamydia: female", "pregnancy"), each = 9),
                                           sim_no = i)
                     setkeyv(flow_df, c("type", "time"))

                     out <- list(flow_df = flow_df, llk = llk)
                   }

stopCluster(c1)

print(Sys.time() - a)

#### Combine all 1000 sets of simulation results
sim_flow <- rbindlist(lapply(outlist, function(x) x$flow_df))

## Calculate posterior mean and 95% credible intervals
sum_flow <- sim_flow[, list(mean_y = mean(y),
                            ub_y = quantile(y, prob = 0.975),
                            lb_y = quantile(y, prob = 0.025)),
                     by = c("type", "time")]

## See results in figures
ggplot(data = sum_flow, aes(x = time, y = mean_y, group = type)) +
  geom_line(data = sum_flow, aes(x = time, y = mean_y),
            colour = "royalblue", size = 1) +
  geom_ribbon(data = sum_flow,
              aes(ymin = lb_y, ymax = ub_y, x = time),
              fill = "lightskyblue", alpha = 0.3) +
  geom_point(data = target_df, aes(x = time, y = y),
             shape = 21, stroke = 1.2, colour = "gray20") +
  facet_wrap(~ type, scales = "free_y") +
  xlab("year") +
  ylab("proportion of population") +
  ggtitle(label = NULL,
          subtitle = "The points are the observed chalmydial incidence by men and women, and the observed pregnancy incidence.\nThe blue lines represent the posterior mean trends generated from the calibrated parameters.\nThe shadded blue areas show the uncertainty around the trends generated from the calibrated parameters.") +
  theme_bw()


