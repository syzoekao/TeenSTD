#' @title Transmission model at a time step
#'
#' @param t time
#' @param initials epidemic states at each time step (population distribution in each of S, I, T, and P states)
#' @param par a list of parameter values
#'
#' @export
transmission_model <- function(t, initials, par) {
  list2env(par, envir = .GlobalEnv)
  list2env(v_d_empty, envir = .GlobalEnv)

  S  <- initials[paste('S', c(1:n_gp), sep = '')]
  Ia <- initials[paste('Ia', c(1:n_gp), sep = '')]
  Is <- initials[paste('Is', c(1:n_gp), sep = '')]
  T  <- initials[paste('T', c(1:n_gp), sep = '')]
  P     <- initials[paste('P', c(1:sub_gp), sep = '')]
  Ia_P  <- initials[paste('Ia_P', c(1:sub_gp), sep = '')]
  Is_P  <- initials[paste('Is_P', c(1:sub_gp), sep = '')]
  Scr_P <- initials[paste('Scr_P', c(1:sub_gp), sep = '')]

  alpha_vec <- c(rep(alpha['alpha_m'], sub_gp), rep(alpha['alpha_f'], sub_gp))

  preg_mix3 <- preg_mix3_apf(t)
  preg_mix4 <- preg_mix4_apf(t)

  R_Inf_nopreg3 <- R_Inf_nopreg3_apf(t)
  R_Inf_nopreg4 <- R_Inf_nopreg4_apf(t)
  R_noInf_preg3 <- R_noInf_preg3_apf(t)
  R_noInf_preg4 <- R_noInf_preg4_apf(t)
  R_Inf_preg3 <- R_Inf_preg3_apf(t)
  R_Inf_preg4 <- R_Inf_preg4_apf(t)

  nu1 = nu[1]
  nu2 = nu[2]
  nu3 = nu[3]
  nu4 = nu[4]

  rktr1 = 0
  rktr2 = 0
  rktr3 = 0
  rktr4 = 0

  psi_t_m <- psi_t_m_apf(t)
  psi_t_f <- psi_t_f_apf(t)
  psi_s_m <- psi_s_m_apf(t)
  psi_s_f <- psi_s_f_apf(t)

  p_male_l <- p_male_l0/sum(c(S[1], Ia[1], Is[1], T[1]))
  p_male_h <- p_male_h0/sum(c(S[2], Ia[2], Is[2], T[2]))
  p_female_l <- p_female_l0/sum(c(S[3], Ia[3], Is[3], T[3], P[1], Ia_P[1], Is_P[1], Scr_P[1]))
  p_female_h <- p_female_h0/sum(c(S[4], Ia[4], Is[4], T[4], P[2], Ia_P[2], Is_P[2], Scr_P[2]))

  # Male: Low risk
  Inft1 <- (R_m1 * p_female_l * (Ia[3] + Is[3] + Ia_P[1] + Is_P[1]) +
              R_m1 * p_female_h * (Ia[4] + Is[4] + Ia_P[2] + Is_P[2])) * S[1]
  a1 <- psi_s_m * pos * Ia[1] + (psi_t_m + psi_s_m) * pos * Is[1]
  dS[1]  <- (nu[1] + rktr2*S[2]) + (v_gamma[['gamma_t']]*T[1]) +
    (v_gamma[['gamma_n']]*Ia[1]) + (v_gamma[['gamma_n']]*Is[1]) -
    Inft1 - (gr + rktr1)*S[1] # dS
  dIa[1] <- rktr2*Ia[2] + alpha_vec[1]*Inft1 -
    psi_s_m*pos*Ia[1] - v_gamma[['gamma_n']]*Ia[1] -
    (gr + rktr1)*Ia[1] # dIa
  dIs[1] <- rktr2*Is[2] + (1 - alpha_vec[1])*Inft1 -
    (psi_t_m + psi_s_m)*pos*Is[1] - v_gamma[['gamma_n']]*Is[1] -
    (gr + rktr1)*Is[1] # dIs
  dT[1]  <- rktr2*T[2] + a1 -
    (v_gamma[['gamma_t']] + gr + rktr1)*T[1] # dT
  d_cum_chl[1] <- a1
  d_cum_inf[1] <- Inft1

  # Male: High risk
  Inft2 <- (R_m2*p_female_l*(Ia[3] + Is[3] + Ia_P[1] + Is_P[1]) +
              R_m2*p_female_h*(Ia[4] + Is[4] + Ia_P[2] + Is_P[2]))*S[2]
  a2 <- psi_s_m*pos*Ia[2] + (psi_t_m + psi_s_m)*pos*Is[2]
  dS[2]  <- (nu[2] + rktr1*S[1]) + (v_gamma[['gamma_t']]*T[2]) +
    (v_gamma[['gamma_n']]*Ia[2]) + (v_gamma[['gamma_n']]*Is[2]) -
    Inft2 - (gr + rktr2)*S[2] # dS
  dIa[2] <- rktr1*Ia[1] + alpha_vec[2]*Inft2 -
    psi_s_m*pos*Ia[2] - v_gamma[['gamma_n']]*Ia[2] -
    (gr + rktr2)*Ia[2] # dIa
  dIs[2] <- rktr1*Is[1] + (1 - alpha_vec[2])*Inft2 -
    (psi_t_m + psi_s_m)*pos*Is[2] - v_gamma[['gamma_n']]*Is[2] -
    (gr + rktr2)*Is[2] # dIs
  dT[2]  <- rktr1*T[1] + a2 -
    v_gamma[['gamma_t']]*T[2] - (gr + rktr2)*T[2] # dT
  d_cum_chl[2] <- a2
  d_cum_inf[2] <- Inft2

  # Female: Low risk
  Inft3_wo_preg <- (R_Inf_nopreg3*p_male_l*(Ia[1] + Is[1]) +
                      R_Inf_nopreg3*p_male_h*(Ia[2] + Is[2]))*S[3]
  Inft3_preg <- (R_Inf_preg3*p_male_l*(Ia[1] + Is[1]) +
                   R_Inf_preg3*p_male_h*(Ia[2] + Is[2]))*S[3]
  Inft3_in_preg <- (R_f1*p_male_l*(Ia[1] + Is[1]) + R_f1*p_male_h*(Ia[2] + Is[2]))*P[1]

  a3_1 <- (psi_s_f*pos*Ia[3] + (psi_t_f + psi_s_f)*pos*Is[3])
  a3_2 <- preg_scr*pos*(Ia_P[1] + Is_P[1])
  b3_1 <- preg_mix3*p_male_l*S[1]*S[3]+preg_mix3*p_male_h*S[2]*S[3] # For assymmetry, we consider the mixing with low and high risk as well.
  b3_2 <- R_noInf_preg3*p_male_l*(Ia[1] + Is[1])*S[3] + R_noInf_preg3*p_male_h*(Ia[2] + Is[2])*S[3]
  b3_3 <- preg_mix3*(sum(c(S[1:2], Ia[1:2], Is[1:2])))*Ia[3]  # assuming that the compartment T is really small and negligible_
  b3_4 <- preg_mix3*(sum(c(S[1:2], Ia[1:2], Is[1:2])))*Is[3]  # assuming that the compartment T is really small and negligible_
  dS[3]  <- (nu[3] + rec_prg*P[1] + rec_prg*Scr_P[1] + rktr4*S[4]) +
    (v_gamma[['gamma_t']]*T[3]) + (v_gamma[['gamma_n']]*Ia[3]) + (v_gamma[['gamma_n']]*Is[3]) -
    (Inft3_wo_preg + Inft3_preg) - b3_1 - b3_2 -
    (gr + rktr3)*S[3] # dS
  dIa[3] <- rktr4*Ia[4] + rec_prg*Ia_P[1] + alpha_vec[3]*Inft3_wo_preg -
    psi_s_f*pos*Ia[3] - v_gamma[['gamma_n']]*Ia[3] -
    b3_3 - (gr + rktr3)*Ia[3] # dIa
  dIs[3] <- rktr4*Is[4] + rec_prg*Is_P[1] + (1 - alpha_vec[3])*Inft3_wo_preg -
    (psi_t_f + psi_s_f)*pos*Is[3] - v_gamma[['gamma_n']]*Is[3] -
    b3_4 - (gr + rktr3)*Is[3] # dIs
  dT[3]  <- rktr4*T[4] + a3_1 - v_gamma[['gamma_t']]*T[3] - (gr + rktr3)*T[3] # dT
  dP[1]  <- b3_1 + v_gamma[['gamma_n']]*(Ia_P[1] + Is_P[1]) + b3_2 + v_gamma[['gamma_t']]*Scr_P[1] -
    (rec_prg + gr)*P[1] - Inft3_in_preg
  dIa_P[1] <- b3_3 + alpha_vec[3]*Inft3_preg + alpha_vec[3]*Inft3_in_preg -
    (rec_prg + preg_scr*pos + v_gamma[['gamma_n']] + gr)*Ia_P[1]
  dIs_P[1] <- b3_4 + (1 - alpha_vec[3])*Inft3_preg + (1 - alpha_vec[3])*Inft3_in_preg -
    (rec_prg + preg_scr*pos + v_gamma[['gamma_n']] + gr)*Is_P[1]
  dScr_P[1] <- a3_2 - (rec_prg + v_gamma[['gamma_t']] + gr)*Scr_P[1]
  d_cum_chl[3] <- a3_1 + a3_2
  d_cum_inf[3] <- Inft3_wo_preg + Inft3_preg + Inft3_in_preg
  d_cum_preg[1] <- b3_1 + b3_2 + b3_3 + b3_4 + Inft3_preg

  # Female: High risk
  Inft4_wo_preg <- (R_Inf_nopreg4*p_male_l*(Ia[1] + Is[1]) +
                      R_Inf_nopreg4*p_male_h*(Ia[2] + Is[2]))*S[4]
  Inft4_preg <- (R_Inf_preg4*p_male_l*(Ia[1] + Is[1]) +
                   R_Inf_preg4*p_male_h*(Ia[2] + Is[2]))*S[4]
  Inft4_in_preg <- (R_f2*p_male_l*(Ia[1] + Is[1]) + R_f2*p_male_h*(Ia[2] + Is[2]))*P[2]

  a4_1 <- (psi_s_f*pos*Ia[4] + (psi_t_f + psi_s_f)*pos*Is[4])
  a4_2 <- preg_scr*pos*(Ia_P[2] + Is_P[2])
  b4_1 <- preg_mix4*p_male_l*S[1]*S[4] + preg_mix4*p_male_h*S[2]*S[4]
  b4_2 <- R_noInf_preg4*p_male_l*(Ia[1] + Is[1])*S[4] + R_noInf_preg4*p_male_h*(Ia[2] + Is[2])*S[4]
  b4_3 <- preg_mix4*(sum(c(S[1:2], Ia[1:2], Is[1:2])))*Ia[4]
  b4_4 <- preg_mix4*(sum(c(S[1:2], Ia[1:2], Is[1:2])))*Is[4]
  dS[4]  <- (nu[4] + rec_prg*P[2] + rec_prg*Scr_P[2] + rktr3*S[3]) +
    (v_gamma[['gamma_t']]*T[4]) + (v_gamma[['gamma_n']]*Ia[4]) + (v_gamma[['gamma_n']]*Is[4]) -
    (Inft4_wo_preg + Inft4_preg) - b4_1 - b4_2 -
    (gr + rktr4)*S[4] # dS
  dIa[4] <- rktr3*Ia[3] + rec_prg*Ia_P[2] + alpha_vec[4]*Inft4_wo_preg -
    psi_s_f*pos*Ia[4] - v_gamma[['gamma_n']]*Ia[4] -
    b4_3 - (gr + rktr4)*Ia[4] # dIa
  dIs[4] <- rktr3*Is[3] + rec_prg*Is_P[2] + (1 - alpha_vec[4])*Inft4_wo_preg -
    (psi_t_f + psi_s_f)*pos*Is[4] - v_gamma[['gamma_n']]*Is[4] -
    b4_4 - (gr + rktr4)*Is[4] # dIs
  dT[4]  <- rktr3*T[3] + a4_1 - v_gamma[['gamma_t']]*T[4] - (gr + rktr4)*T[4] # dT
  dP[2]  <- b4_1 + v_gamma[['gamma_n']]*(Ia_P[2] + Is_P[2]) + b4_2 + v_gamma[['gamma_t']]*Scr_P[2] -
    (rec_prg + gr)*P[2] - Inft4_in_preg
  dIa_P[2] <- b4_3 + alpha_vec[4]*Inft4_preg + alpha_vec[4]*Inft4_in_preg -
    (rec_prg + preg_scr*pos + v_gamma[['gamma_n']] + gr)*Ia_P[2]
  dIs_P[2] <- b4_4 + (1 - alpha_vec[4])*Inft4_preg + (1 - alpha_vec[4])*Inft4_in_preg -
    (rec_prg + preg_scr*pos + v_gamma[['gamma_n']] + gr)*Is_P[2]
  dScr_P[2] <- a4_2 - (rec_prg + v_gamma[['gamma_t']] + gr)*Scr_P[2]
  d_cum_chl[4] <- a4_1 + a4_2
  d_cum_inf[4] <- Inft4_wo_preg + Inft4_preg + Inft4_in_preg
  d_cum_preg[2] <- b4_1 + b4_2 + b4_3 + b4_4 + Inft4_preg

  list(c(dS, dIa, dIs, dT, dP, dIa_P, dIs_P, dScr_P, d_cum_chl, d_cum_inf, d_cum_preg))
}

#' @export
process_output <- function(sim_out, par) {
  max_time <- par$max_time
  n_gp <- par$n_gp

  tot_male <- rowSums(sim_out[, c(rep(paste(c("S"), as.character(c(1:2)), sep = "")),
                              rep(paste(c("Ia"), as.character(c(1:2)), sep = "")),
                              rep(paste(c("Is"), as.character(c(1:2)), sep = "")),
                              rep(paste(c("T"), as.character(c(1:2)), sep = "")))])
  # print("======1")

  tot_female <- rowSums(sim_out[, c(rep(paste(c("S"), as.character(c(3:4)), sep = "")),
                                rep(paste(c("Ia"), as.character(c(3:4)), sep = "")),
                                rep(paste(c("Is"), as.character(c(3:4)), sep = "")),
                                rep(paste(c("T"), as.character(c(3:4)), sep = "")),
                                rep(paste(c("P"), as.character(c(1:2)), sep = "")),
                                rep(paste(c("Ia_P"), as.character(c(1:2)), sep = "")),
                                rep(paste(c("Is_P"), as.character(c(1:2)), sep = "")),
                                rep(paste(c("Scr_P"), as.character(c(1:2)), sep = "")) )])


  flux <- matrix(NA, nrow = nrow(sim_out), ncol = (2*n_gp + 5), byrow = TRUE)
  colnames(flux) <- c("time", rep(paste("fluxM", as.factor(c(1:2)), sep = "")),
                      rep(paste("fluxF", as.factor(c(1:2)), sep = "")),
                      "tot_fluxM", "tot_fluxF",
                      rep(paste0("infM", c(1:2))),
                      rep(paste0("infF", c(1:2))),
                      "tot_infM", "tot_infF")
  flux[, "time"] <- sim_out[, "time"]
  flux[1, ] <- 0

  flux[, "fluxM1"] <- c(0, diff(sim_out[, "cum_chl1"]) / tot_male[c(2:length(tot_male))])
  flux[, "fluxM2"] <- c(0, diff(sim_out[, "cum_chl2"]) / tot_male[c(2:length(tot_male))])
  flux[, "fluxF1"] <- c(0, diff(sim_out[, "cum_chl3"]) / tot_female[c(2:length(tot_female))])
  flux[, "fluxF2"] <- c(0, diff(sim_out[, "cum_chl4"]) / tot_female[c(2:length(tot_female))])
  flux[, "infM1"] <- c(0, diff(sim_out[, "cum_inf1"]) / tot_male[c(2:length(tot_male))])
  flux[, "infM2"] <- c(0, diff(sim_out[, "cum_inf2"]) / tot_male[c(2:length(tot_male))])
  flux[, "infF1"] <- c(0, diff(sim_out[, "cum_inf3"]) / tot_female[c(2:length(tot_female))])
  flux[, "infF2"] <- c(0, diff(sim_out[, "cum_inf4"]) / tot_female[c(2:length(tot_female))])

  flux[, "tot_fluxM"] <- rowSums(flux[, c("fluxM1", "fluxM2")])
  flux[, "tot_fluxF"] <- rowSums(flux[, c("fluxF1", "fluxF2")])
  flux[, "tot_infM"] <- rowSums(flux[, c("infM1", "infM2")])
  flux[, "tot_infF"] <- rowSums(flux[, c("infF1", "infF2")])

  preg_mat <- matrix(NA, nrow = nrow(sim_out), ncol = 4, byrow = TRUE)
  colnames(preg_mat) <- c("time", "preg_low", "preg_high", "preg_tot")
  preg_mat[, "time"] <- sim_out[, "time"]
  preg_mat[1, ] <- 0
  preg_mat[, "preg_low"] <- c(0, diff(sim_out[, "cum_preg1"]) / tot_female[c(2:length(tot_female))])
  preg_mat[, "preg_high"] <- c(0, diff(sim_out[, "cum_preg2"]) / tot_female[c(2:length(tot_female))])
  preg_mat[, "preg_tot"] <- rowSums(preg_mat[, c("preg_low", "preg_high")])

  IaM <- sim_out[, "Ia1"]*rowSums(sim_out[, paste(c("S", "Ia", "Is", "T"), 1, sep = "")]) +
    sim_out[, "Ia2"]*rowSums(sim_out[, paste(c("S", "Ia", "Is", "T"), 2, sep = "")])
  IsM <- sim_out[, "Is1"]*rowSums(sim_out[, paste(c("S", "Ia", "Is", "T"), 1, sep = "")]) +
    sim_out[, "Is2"]*rowSums(sim_out[, paste(c("S", "Ia", "Is", "T"), 2, sep = "")])
  IaF <- rowSums(sim_out[, c("Ia3", "Ia_P1")])*rowSums(sim_out[, c("S3", "Ia3", "Is3", "T3", "P1", "Ia_P1", "Is_P1", "Scr_P1")]) +
    rowSums(sim_out[, c("Ia4", "Ia_P2")])*rowSums(sim_out[, c("S4", "Ia4", "Is4", "T4", "P2", "Ia_P2", "Is_P2", "Scr_P2")])
  IsF <- rowSums(sim_out[, c("Is3", "Is_P1")])*rowSums(sim_out[, c("S3", "Ia3", "Is3", "T3", "P1", "Ia_P1", "Is_P1", "Scr_P1")]) +
    rowSums(sim_out[, c("Is4", "Is_P2")])*rowSums(sim_out[, c("S4", "Ia4", "Is4", "T4", "P2", "Ia_P2", "Is_P2", "Scr_P2")])

  sol <- cbind(sim_out[, "time"],
               sim_out[, paste(rep(c("S", "Ia", "Is", "T"), each = 2), c(1:2), sep = "")]/tot_male,
               sim_out[, c(paste(rep(c("S", "Ia", "Is", "T"), each = 2), c(3:4), sep = ""),
                       paste(rep(c("P", "Ia_P", "Is_P", "Scr_P"), each = 2), c(1:2), sep = ""))]/tot_female,
               IaM, IsM, IaF, IsF,
               flux[, c("tot_fluxM", "tot_fluxF", "tot_infM", "tot_infF")],
               preg_mat[, "preg_tot"])

  f_year <- 2013
  burnin <- par$burnin
  years <- par$years
  # print("======3")
  flux_df <- data.frame(cbind(flux[-1, ], group = rep(c((f_year - (burnin + years) + 1):f_year), each = 12)))
  yr_fluxM <- aggregate(flux_df$tot_fluxM, by=list(Category = flux_df$group), FUN=sum)
  yr_fluxF <- aggregate(flux_df$tot_fluxF, by=list(Category = flux_df$group), FUN=sum)
  yr_infM <- aggregate(flux_df$tot_infM, by=list(Category = flux_df$group), FUN=sum)
  yr_infF <- aggregate(flux_df$tot_infF, by=list(Category = flux_df$group), FUN=sum)
  flux_pred <- cbind(yr_fluxM$x, yr_fluxF$x, yr_infM$x, yr_infF$x)

  preg_df <- data.frame(cbind(preg_mat[-1, ], group = rep(c((f_year - (burnin + years) + 1):f_year), each = 12)))
  yr_preg <- aggregate(preg_df$preg_tot, by=list(Category = preg_df$group), FUN=sum)
  preg_pred <- yr_preg$x

  return(list(fluxM = yr_fluxM$x[c((burnin + 1) : (max_time / 12))],
              fluxF = yr_fluxF$x[c((burnin + 1) : (max_time / 12))],
              fluxP = yr_preg$x[c((burnin + 1) : (max_time / 12))],
              flow = flux_pred, preg = preg_pred, sol = sol))
}




