% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\name{create_parameters}
\alias{create_parameters}
\title{Create parameters for the TeenSTD model}
\usage{
create_parameters(
  days = 30,
  month = 12,
  n_gp = 4,
  years = 9,
  burnin = 0,
  treat_fact_f = 83.77591,
  treat_fact_m = 67.51961,
  screen_fact_m = 0.142667,
  beta = 0.09999317,
  factor = 1.62593,
  alpha_m = 0.8800686,
  incubate = 18.46923,
  fr_sex = 5.186936,
  p_condom = 0.6825322,
  p_condom_l = 0.9998123,
  r_preg = 0.05181339,
  red_condom = 0.1656216,
  oral_h10 = 0.6561537,
  ratio_oral_mean_over_h10 = 0.2630736/0.4051145,
  IUD_h10 = 0.2456997,
  ratio_IUD_mean_over_h10 = 0.09600559/0.2456997,
  IUD_factor = 0.8340665,
  I_m = 0.0208744,
  I_f = 0.09887864,
  T_m = 7.260635e-05,
  T_f = 0.06556457,
  sens = 0.93,
  spec = 0.94,
  pos = 0.93,
  pInsure = 0.97,
  fail_condom = 0.1122,
  fail_oral = 0.0535,
  fail_IUD = 2e-04,
  fail_condom_IUD = 2.58e-06,
  fail_condom_oral = 9e-04,
  model_type = "dual"
)
}
\arguments{
\item{days}{number of days per month}

\item{month}{number of months per year}

\item{n_gp}{total number of groups including two genders (men and women) and two risk groups (high and low)}

\item{years}{number of years of simulation (2005-2013)}

\item{burnin}{burnin period in years (not part of the `years` parameter)}

\item{treat_fact_f}{time to seek treatment among women (in days)}

\item{treat_fact_m}{time to seek treatment among men (in days)}

\item{screen_fact_m}{factor adjusting the probability that a man is screened for CT, modeled as a proportion
of the probability of CT screening among women}

\item{beta}{transmission rate from men to women per sex act}

\item{factor}{factor adjusting the transmission rate from women to men, ranging from 0 to 1}

\item{alpha_m}{proportion of asymptomatic CT infection (same for both men and women)}

\item{incubate}{incubation period in days}

\item{fr_sex}{frequency of sex per month}

\item{p_condom}{average probability of condom use per sex act}

\item{p_condom_l}{probability of condom use per sex act among low risk group}

\item{r_preg}{rate of pregnancy per sex act}

\item{red_condom}{reduction of transmission due to condom use}

\item{oral_h10}{probability of oral contraception use among women who are at high risk at the end of year 9}

\item{ratio_oral_mean_over_h10}{average probability of oral contraception use among women is a
proportion of the probability of oral contraception use among women at high risk}

\item{IUD_h10}{probability of IUD use among women at high risk at the end of year 9}

\item{ratio_IUD_mean_over_h10}{average probability of IUD use among women is a
proportion of the probability of IUD use among women at high risk}

\item{IUD_factor}{the probability of IUD use among women at year t is calculated as the probability
of IUD use of year t+1 multiplied by `IUD_factor`}

\item{I_m}{proportion of men in infected state initially}

\item{I_f}{proportion of women in infected state initially}

\item{T_m}{proportion of men in treated state initially}

\item{T_f}{proportion of women in treated state initially}

\item{sens}{sensitivity of CT testing}

\item{spec}{specificity of CT testing}

\item{pos}{probability of testing positive among the infected}

\item{pInsure}{proportion of the adolescent who are insured}

\item{fail_condom}{probability of failed condom use per sex act}

\item{fail_oral}{probability of failed oral contraception per sex act}

\item{fail_IUD}{probability of failed IUD use per sex act}

\item{fail_condom_IUD}{probability of failed condom use and IUD use per sex act}

\item{fail_condom_oral}{probability of failed condom use and oral contraception per sex act}

\item{model_type}{whether the model considers both chlamydial infections and pregnancy (`"dual"`) or
only chlamydial infections}
}
\description{
Create parameters for the TeenSTD model
}
