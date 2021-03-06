#' @title Create parameters for the TeenSTD model
#'
#' @param days number of days per month
#' @param month number of months per year
#' @param n_gp total number of groups including two genders (men and women) and two risk groups (high and low)
#' @param years number of years of simulation (2005-2013)
#' @param burnin burnin period in years (not part of the `years` parameter)
#' @param treat_fact_f time to seek treatment among women (in days)
#' @param treat_fact_m time to seek treatment among men (in days)
#' @param screen_fact_m factor adjusting the probability that a man is screened for CT, modeled as a proportion
#'        of the probability of CT screening among women
#' @param beta transmission rate from men to women per sex act
#' @param factor factor adjusting the transmission rate from women to men, ranging from 0 to 1
#' @param alpha_m proportion of asymptomatic CT infection (same for both men and women)
#' @param incubate incubation period in days
#' @param fr_sex frequency of sex per month
#' @param p_condom average probability of condom use per sex act
#' @param p_condom_l probability of condom use per sex act among low risk group
#' @param r_preg rate of pregnancy per sex act
#' @param red_condom reduction of transmission due to condom use
#' @param oral_h10 probability of oral contraception use among women who are at high risk at the end of year 9
#' @param ratio_oral_mean_over_h10 average probability of oral contraception use among women is a
#'        proportion of the probability of oral contraception use among women at high risk
#' @param IUD_h10 probability of IUD use among women at high risk at the end of year 9
#' @param ratio_IUD_mean_over_h10 average probability of IUD use among women is a
#'        proportion of the probability of IUD use among women at high risk
#' @param IUD_factor the probability of IUD use among women at year t is calculated as the probability
#'        of IUD use of year t+1 multiplied by `IUD_factor`
#' @param I_m proportion of men in infected state initially
#' @param I_f proportion of women in infected state initially
#' @param T_m proportion of men in treated state initially
#' @param T_f proportion of women in treated state initially
#' @param sens sensitivity of CT testing
#' @param spec specificity of CT testing
#' @param pos probability of testing positive among the infected
#' @param pInsure proportion of the adolescent who are insured
#' @param fail_condom probability of failed condom use per sex act
#' @param fail_oral probability of failed oral contraception per sex act
#' @param fail_IUD probability of failed IUD use per sex act
#' @param fail_condom_IUD probability of failed condom use and IUD use per sex act
#' @param fail_condom_oral probability of failed condom use and oral contraception per sex act
#' @param model_type whether the model considers both chlamydial infections and pregnancy (`"dual"`) or
#'        only chlamydial infections
#'
#' @export
create_parameters <- function(days = 30,
                              month = 12,
                              n_gp = 4,
                              years = 9,
                              burnin = 0,
                              treat_fact_f = 83.77591, #calibrate
                              treat_fact_m = 67.51961, #calibrate
                              screen_fact_m = 0.142667, #calibrate
                              beta = 0.09999317, #calibrate
                              factor = 1.62593, #calibrate
                              alpha_m = 0.8800686, #calibrate
                              incubate = 18.46923, #calibrate
                              fr_sex = 5.186936, #calibrate
                              p_condom = 0.6825322, #calibrate
                              p_condom_l = 0.9998123, #calibrate
                              r_preg = 0.05181339, #calibrate
                              red_condom = 0.1656216, #calibrate
                              oral_h10 = 0.6561537, #calibrate
                              ratio_oral_mean_over_h10 = 0.2630736 / 0.4051145, #calibrate
                              IUD_h10 = 0.2456997, #calibrate
                              ratio_IUD_mean_over_h10 = 0.09600559 / 0.2456997, #calibrate
                              IUD_factor = 0.8340665, #calibrate
                              I_m = 0.0208744, #calibrate
                              I_f = 0.09887864, #calibrate
                              T_m = 7.260635e-05, #calibrate
                              T_f = 0.06556457, #calibrate
                              sens = 0.93,
                              spec = 0.94,
                              pos = 0.93,
                              pInsure = 0.97,
                              fail_condom = 0.1122,
                              fail_oral = 0.0535,
                              fail_IUD = 0.0002,
                              fail_condom_IUD = 0.00000258,
                              fail_condom_oral = 0.0009,
                              model_type = "dual") {


  if(model_type != "dual") {
    r_preg <- 0
    oral_h10 <- 0
    IUD_h10 <- 0
    IUD_factor <- 0
  }

  #### Number of time steps in the calibration time window
  cycles <- month * years + 1

  #### Time horizon including the burnin period
  max_time <- month * (burnin + years)

  #### transition between risk groups (this is not used in the model)
  tmp_name <- c(paste("mo", c(0, 13, 49, 85), sep = ""))
  risk_tr_ary <- sex_behave_data$risk.tr.ary[c(1:2), , ]
  risk_tr_m <- -log(1 - matrix(rep(risk_tr_ary[, 4, "male"], 4), ncol = 2, byrow = T))
  risk_tr_f <- -log(1 - matrix(rep(risk_tr_ary[, 4, "female"], 4), ncol = 2, byrow = T))

  p_screen_f <- c(rep(0.3508, burnin),
                  c(0.3580, 0.3786, 0.4204, 0.4660, 0.4764, 0.4799, 0.4745, 0.4629, 0.4486, 0.4638))
  names(p_screen_f) <- c((2013 - (years + burnin)) : 2013)
  p_screen_m <- screen_fact_m * p_screen_f

  n_partner <- list(low_m = 0.470, high_m = 3.122, low_f = 0.568, high_f = 3.586)

  avg_dist_risk <- sex_behave_data$avg.dist.risk.ary
  avg_l <- n_partner$low_f * avg_dist_risk[1, 1, "female"]
  avg_h <- n_partner$high_f * avg_dist_risk[2, 1, "female"]
  p_female_l <- avg_l / (avg_l + avg_h)
  p_female_h <- 1 - p_female_l

  avg_l <- n_partner$low_m * avg_dist_risk[1, 1, "male"]
  avg_h <- n_partner$high_m * avg_dist_risk[2, 1, "male"]
  p_male_l <- avg_l / (avg_l + avg_h)
  p_male_h <- 1 - p_male_l
  rr_pop_female <- 1


  pop_low <- (avg_dist_risk[1, 1, "female"] * rr_pop_female +
                    avg_dist_risk[1, 1, "male"] * 1) / (rr_pop_female + 1)
  pop_high <- 1 - pop_low


  #### Contraception methods
  mean_oral <- oral_h10 * ratio_oral_mean_over_h10
  mean_IUD <- IUD_h10 * ratio_IUD_mean_over_h10

  oral_l10 <- (mean_oral - avg_dist_risk[2, 1, "female"] * oral_h10) / avg_dist_risk[1, 1, "female"]
  oral_l10 <- (oral_l10 <= 0) * 0 + (oral_l10 > 0) * oral_l10

  IUD_l10 <- (mean_IUD - avg_dist_risk[2, 1, "female"] * IUD_h10) / avg_dist_risk[1, 1, "female"]
  IUD_l10 <- (IUD_l10 <= 0) * 0 + (IUD_l10 > 0) * IUD_l10

  p_bc <- c(mean = mean_oral + mean_IUD,
            low = oral_l10 + IUD_l10,
            high = oral_h10 + IUD_h10)
  max_condom_use <- 1 - p_bc

  p_condom_h <- (p_condom - pop_low * p_condom_l) / pop_high
  p_condom_h <- (p_condom_h < 0) * 0 + (p_condom_h >= 0) * p_condom_h

  v_IUD_l <- IUD_l10 * (IUD_factor ^ c(years : 0))
  p_IUD_l <- approx(c(rep(v_IUD_l[1], burnin), v_IUD_l),
                        n = (month * burnin + cycles))$y

  v_IUD_h <- IUD_h10 * (IUD_factor ^ c(years : 0))
  p_IUD_h <- approx(c(rep(v_IUD_h[1], burnin), v_IUD_h),
                        n = (month * burnin + cycles))$y
  p_IUD <- cbind(low = p_IUD_l, high = p_IUD_h)

  p_oral_l <- rep(oral_l10, month * burnin + cycles)
  p_oral_h <- rep(oral_h10, month * burnin + cycles)
  p_oral <- cbind(low = p_oral_l, high = p_oral_h)

  alpha_f <- alpha_m

  pTreat_f <- 1 / ((incubate + treat_fact_f) / days)
  pTreat_m <- 1 / ((incubate + treat_fact_m) / days)

  p_condom_l0 <- rep(p_condom_l, (burnin + years) * month + 1)
  p_condom_h0 <- rep(p_condom_h, (burnin + years) * month + 1)
  p_condom <- cbind(low = p_condom_l0, high = p_condom_h0)

  n_partner_mo <- lapply(n_partner, function(x) x / month)

  psi_s_m <- (1/month) * p_screen_m * pInsure
  psi_s_f <- (1/month) * p_screen_f * pInsure
  psi_t_m <- pTreat_m * pInsure
  psi_t_f <- pTreat_f * pInsure

  v_gamma <- c(gamma_n = 1 / month, gamma_t = 1 / (7 / days))

  gr <- 1 / (5 * 12)

  p_dual <- (p_IUD + p_oral) * p_condom
  dual_condom_oral <- p_oral * p_condom
  dual_condom_IUD <- p_IUD * p_condom

  only_condom <- p_condom - p_dual
  only_oral <- p_oral - dual_condom_oral
  only_IUD <- p_IUD - dual_condom_IUD

  preg_lambda <- c(nobc = r_preg,
                   only_condom = r_preg * fail_condom,
                   only_oral = r_preg * fail_oral,
                   only_IUD = r_preg * fail_IUD,
                   condom_oral = r_preg * fail_condom_oral,
                   condom_IUD = r_preg * fail_condom_IUD)

  prob_bc <- p_oral + p_IUD

  if (model_type == "dual") {
    p_condom_given_bc <- p_dual / prob_bc
    p_condom_and_nobc <- p_condom * (1 - prob_bc)
    p_condom_given_nobc <- p_condom_and_nobc / (1 - prob_bc)
  } else {
    p_condom_given_bc <- p_condom
    p_condom_given_bc[, 1] <- 0
    p_condom_given_bc[, 2] <- 0
    p_condom_and_nobc <- p_condom
    p_condom_given_nobc <- p_condom
  }

  rec_prg <- 1/((7.41+7.44+7.47+7.54+7.5+7.64+7.49+7.63+7.71)/9)
  preg_scr <- 1 / 3
  alpha <- c(alpha_m = alpha_m, alpha_f = alpha_f)

  R_m <- -log(1-(1-((1-beta) ^ ((1-p_condom) * fr_sex) *
                      (1-beta*(1-red_condom))^((p_condom)*fr_sex))))
  R_bc_f <- -log(1-(1-((1-(beta * factor))^((1-p_condom_given_bc)*fr_sex)*
                         (1-(beta*factor)*(1-red_condom))^((p_condom_given_bc)*fr_sex))))
  R_nobc_f <- -log(1-(1-((1-(beta*factor))^((1-p_condom_given_nobc)*fr_sex)*
                               (1-(beta*factor)*(1-red_condom))^((p_condom_given_nobc)*fr_sex))))
  R_f <- cbind(low = (prob_bc[, "low"]*R_bc_f[, "low"] +
                            (1-prob_bc[, "low"])*R_nobc_f[, "low"]),
               high = (prob_bc[, "high"]*R_bc_f[, "high"] +
                         (1-prob_bc[, "high"])*R_nobc_f[, "high"]))

  R_m <- cbind(low = n_partner_mo[["low_m"]]*R_m[, "low"],
               high = n_partner_mo[["high_m"]]*R_m[, "high"])
  R_f <- cbind(low = n_partner_mo[["low_f"]]*R_f[, "low"],
               high = n_partner_mo[["high_f"]]*R_f[, "high"])

  p_all_bc <- only_condom + only_IUD + only_oral +
    dual_condom_IUD + dual_condom_oral
  preg_mon_nobc <- 1 - ((1 - preg_lambda["nobc"]) ^ (fr_sex * (1-p_condom_given_nobc))*
                        ((1 - preg_lambda["only_condom"]) ^ (fr_sex * p_condom_given_nobc)))
  preg_mon_bc <- 1 - (((1 - preg_lambda["only_oral"]) ^ (fr_sex * (only_oral/prob_bc))) *
                        ((1-preg_lambda["only_IUD"]) ^ (fr_sex*(only_IUD/prob_bc)))*
                      ((1-preg_lambda["condom_oral"]) ^ (fr_sex * (dual_condom_oral/prob_bc))) *
                        ((1-preg_lambda["condom_IUD"]) ^ (fr_sex*(dual_condom_IUD/prob_bc))))
  r_preg_mon_nobc <- -log(1-preg_mon_nobc)
  r_preg_mon_bc <- -log(1-preg_mon_bc)

  preg_mix1 <- n_partner_mo[["low_f"]] * (prob_bc[, "low"] * r_preg_mon_bc[, "low"] +
                                          (1-prob_bc[, "low"])*r_preg_mon_nobc[, "low"])
  preg_mix2 <- n_partner_mo[["high_f"]] * (prob_bc[, "high"] * r_preg_mon_bc[, "high"] +
                                           (1-prob_bc[, "high"])*r_preg_mon_nobc[, "high"])
  preg_mix_mon <- cbind(low = preg_mix1, high = preg_mix2)

  R_Inf_nopreg <- -log(1-(1-exp(-R_f))*exp(-preg_mix_mon))
  R_noInf_preg<- -log(1-exp(-R_f)*(1-exp(-preg_mix_mon)))
  R_Inf_preg <- -log(1-(1-exp(-R_f))*(1-exp(-preg_mix_mon)))

  dist_nsex_male_age1 <- avg_dist_risk[c(1:2), 1, "male"]
  dist_nsex_female_age1 <- avg_dist_risk[c(1:2), 1, "female"]

  dpop <- c(dist_nsex_male_age1, dist_nsex_female_age1 * rr_pop_female)
  names(dpop) <- c("m:lo", "m:hi", "f:lo", "f:hi")
  n_gp <- n_gp
  sub_gp <- n_gp / 2
  I_male <- I_m * dpop[1:(n_gp/2)]
  I_female <- I_f * dpop[(n_gp/2+1):n_gp]
  T_male <- T_m * dpop[1:(n_gp/2)]
  T_female <- T_f * dpop[(n_gp/2+1):n_gp]

  Ia <- c(I_male * alpha_m, I_female * alpha_f)
  Is <- c(I_male * (1 - alpha_m), I_female * (1 - alpha_f))
  T <- c(T_male, T_female)
  P <- 0 * dist_nsex_female_age1
  Ia_P <- rep(0, sub_gp)
  Is_P <- rep(0, sub_gp)
  Scr_P <- rep(0, sub_gp)
  S <- dpop - (Ia + Is + T + c(0, 0, P))
  cum_chl <- rep(0, n_gp)
  cum_inf <- rep(0, n_gp)
  cum_preg <- rep(0, sub_gp)

  initials <- c(S = S,
                Ia = Ia,
                Is = Is,
                T = T,
                P = P,
                Ia_P = Ia_P,
                Is_P = Is_P,
                Scr_P = Scr_P,
                cum_chl = cum_chl,
                cum_inf = cum_inf,
                cum_preg = cum_preg)

  names(initials) <- c(rep(paste(c("S"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("Ia"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("Is"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("T"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("P"), as.character(c(1:sub_gp)), sep = "")),
                       rep(paste(c("Ia_P"), as.character(c(1:sub_gp)), sep = "")),
                       rep(paste(c("Is_P"), as.character(c(1:sub_gp)), sep = "")),
                       rep(paste(c("Scr_P"), as.character(c(1:sub_gp)), sep = "")),
                       rep(paste(c("cum_chl"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("cum_inf"), as.character(c(1:n_gp)), sep = "")),
                       rep(paste(c("cum_preg"), as.character(c(1:sub_gp)), sep = "")))
  n_state <- length(c("S", "Ia", "Is", "T", "P", "Ia_P", "Is_P", "Scr_P"))

  v_nu <- gr * dpop

  #### Create .apf
  tmp_time <- c(0, month * c(1 : (burnin + years)) + 1)
  psi_t_m_apf <- approxfun(tmp_time, rep(psi_t_m, length(tmp_time)), method = "linear")
  psi_s_m_apf <- approxfun(tmp_time, psi_s_m, method = "linear")
  psi_t_f_apf <- approxfun(tmp_time, rep(psi_t_f, length(tmp_time)), method = "linear")
  psi_s_f_apf <- approxfun(tmp_time, psi_s_f, method = "linear")

  tmp_time <- c(0 : (month * (burnin + years) + 1))
  R_Inf_nopreg3_apf <- approxfun(tmp_time, c(R_Inf_nopreg[, 1], R_Inf_nopreg[nrow(R_Inf_nopreg), 1]))
  R_Inf_nopreg4_apf <- approxfun(tmp_time, c(R_Inf_nopreg[, 2], R_Inf_nopreg[nrow(R_Inf_nopreg), 2]))
  R_noInf_preg3_apf <- approxfun(tmp_time, c(R_noInf_preg[, 1], R_Inf_preg[nrow(R_noInf_preg), 1]))
  R_noInf_preg4_apf <- approxfun(tmp_time, c(R_noInf_preg[, 2], R_Inf_preg[nrow(R_noInf_preg), 2]))
  R_Inf_preg3_apf <- approxfun(tmp_time, c(R_Inf_preg[, 1], R_Inf_preg[nrow(R_Inf_preg), 1]))
  R_Inf_preg4_apf <- approxfun(tmp_time, c(R_Inf_preg[, 2], R_Inf_preg[nrow(R_Inf_preg), 2]))
  preg_mix3_apf <- approxfun(tmp_time, c(preg_mix_mon[, 1], preg_mix_mon[nrow(preg_mix_mon), 1]))
  preg_mix4_apf <- approxfun(tmp_time, c(preg_mix_mon[, 2], preg_mix_mon[nrow(preg_mix_mon), 2]))

  apfs <- list(R_Inf_nopreg3_apf = R_Inf_nopreg3_apf, R_Inf_nopreg4_apf = R_Inf_nopreg4_apf,
               R_noInf_preg3_apf = R_noInf_preg3_apf, R_noInf_preg4_apf = R_noInf_preg4_apf,
               R_Inf_preg3_apf = R_Inf_preg3_apf, R_Inf_preg4_apf = R_Inf_preg4_apf,
               preg_mix3_apf = preg_mix3_apf, preg_mix4_apf = preg_mix4_apf,
               psi_t_m_apf = psi_t_m_apf, psi_t_f_apf = psi_t_f_apf,
               psi_s_m_apf = psi_s_m_apf, psi_s_f_apf = psi_s_f_apf)

  par <- list(preg_mix_mon = preg_mix_mon,
              dpop = dpop,
              R_m1 = R_m[100, "low"], R_m2 = R_m[100, "high"],
              R_f1 = R_f[100, "low"], R_f2 = R_f[100, "high"],
              alpha = alpha,
              pos = pos,
              v_gamma = v_gamma,
              nu = v_nu,
              mu = rep(0, 4),
              gr = gr,
              preg_scr = preg_scr,
              rec_prg = rec_prg,
              n_gp = n_gp,
              sub_gp = sub_gp,
              p_male_l0 = p_male_l,
              p_male_h0 = p_male_h,
              p_female_l0 = p_female_l,
              p_female_h0 = p_female_h,
              days = days,
              month = month,
              years = years,
              burnin = burnin,
              max_time = max_time,
              n_state = n_state)
  par <- c(par, apfs)

  par$v_d_empty <- list(dS = rep(NA, n_gp),
                    dIa = rep(NA, n_gp),
                    dIs = rep(NA, n_gp),
                    dT = rep(NA, n_gp),
                    dP = rep(NA, sub_gp),
                    dIa_P = rep(NA, sub_gp),
                    dIs_P = rep(NA, sub_gp),
                    dScr_P = rep(NA, sub_gp),
                    d_cum_chl = rep(NA, n_gp),
                    d_cum_inf = rep(NA, n_gp),
                    d_cum_preg = rep(NA, sub_gp))

  return(list(par = par, initials = initials))
}


##Created dataframe with parameter ranges for calibrated model inputs
# parm_name<-(c("treat_fact_f", "treat_fact_m", "beta","factor",
#               "alpha_m", "incubate", "screen_fact_m", "fr_sex",
#               "p_condom", "p_condom_l", "r_preg", "red_condrom",
#               "oral_h10","ratio_oral_mean_over_h10",
#               "IUD_h10","ratio_IUD_mean_over_h10", "oral_fact",
#               "IUD_factor", "I_m", "T_m", "I_f", "T_f"))
# parm_lb<-(c(1,1,0.03,1,0.65,1,0,0,0.329,0.329,0,0.15,0.143,NA,0.056,NA,1,0.5,0.0176,NA,0.0176,NA))
# parm_ub<-(c(120,120,0.8,2,0.96,21,0.2,50,0.728,1,0.094,0.9,0.266,NA,0.104,NA,1,1,0.1464,NA,0.1464,NA))
# parm_df<-data.frame(parm_name, parm_lb, parm_ub)

