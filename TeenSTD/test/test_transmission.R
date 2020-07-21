library(TeenSTD)
library(data.table)
library(ggplot2)

parmdf <- readRDS("TeenSTD/test/posterior_set.RDS")

colnames(parmdf) <- c("treat_fact_f", "treat_fact_m", "screen_fact_m", "beta",
                      "factor", "alpha_m", "incubate", "fr_sex", "p_condom",
                      "p_condom_l", "r_preg", "red_condom", "mean_oral", "oral_h10",
                      "mean_IUD", "IUD_h10", "IUD_factor", "I_m", "I_f", "T_m", "T_f", "ll")
samp_i <- sample(c(1:nrow(parmdf)), 1)
parmdf <- parmdf[samp_i, ]
print(parmdf)



tmp_ls <- create_parameters(burnin = 2,
                            treat_fact_f = parmdf$treat_fact_f,
                            treat_fact_m = parmdf$treat_fact_m,
                            screen_fact_m = parmdf$screen_fact_m,
                            beta = parmdf$beta,
                            factor = parmdf$factor,
                            alpha_m = parmdf$alpha_m,
                            incubate = parmdf$incubate,
                            fr_sex = parmdf$fr_sex,
                            p_condom = parmdf$p_condom,
                            p_condom_l = parmdf$p_condom_l,
                            r_preg = parmdf$r_preg,
                            red_condom = parmdf$red_condom,
                            oral_h10 = parmdf$oral_h10,
                            ratio_oral_mean_over_h10 = parmdf$mean_oral / parmdf$oral_h10,
                            IUD_h10 = parmdf$IUD_h10,
                            ratio_IUD_mean_over_h10 = parmdf$mean_IUD / parmdf$IUD_h10,
                            IUD_factor = parmdf$IUD_factor,
                            I_m = parmdf$I_m,
                            I_f = parmdf$I_f,
                            T_m = parmdf$T_m,
                            T_f = parmdf$T_f)
par <- tmp_ls$par
initials <- tmp_ls$initials


system.time(chl <- ode(times = seq(0, par$max_time, 1),
                       y = initials,
                       parms = par,
                       func = transmission_model))

system.time(output <- process_output(chl, par))

target <- c(MDH_data$obsM1[1:9], MDH_data$obsF1[1:9], preg_data$obs[1:9])

flow_df <- data.table(time = rep(c(2005:2013), 3),
                      y = c(output$fluxM, output$fluxF, output$fluxP),
                      type = rep(c("chlamydia: male", "chlamydia: female", "pregnancy"), each = 9))
setkeyv(flow_df, c("type", "time"))

target_df <- data.table(time = rep(c(2005:2013), 3),
                        y = target,
                        type = rep(c("chlamydia: male", "chlamydia: female", "pregnancy"), each = 9))
setkeyv(target_df, c("type", "time"))

# flow_df <- rbindlist(list(flow_df, target_df))

ggplot(data = flow_df, aes(x = time, y = y, group = type)) +
  geom_line(data = flow_df, aes(x = time, y = y)) +
  geom_point(data = target_df, aes(x = time, y = y)) +
  facet_wrap(~ type) +
  theme_bw()

