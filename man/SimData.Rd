% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimData.R
\name{SimData}
\alias{SimData}
\title{Data generation for simulation}
\usage{
SimData(
  sim,
  n,
  tL,
  tU,
  nbreak = 8,
  p01 = 1,
  p02 = 1,
  q_eta = 1,
  t_max = 19,
  k_max = 2,
  distr = "Gaussian",
  m_age = 0,
  std_age = 2,
  sigmau_inv,
  tbtheta,
  tbeta0,
  tbeta1,
  tsigma,
  tsigmad,
  ttheta,
  teta,
  tgamma,
  lambda0 = 0.08
)
}
\arguments{
\item{sim}{number of simulated datasets to be generated.}

\item{n}{sample size for each dataset.}

\item{tL}{lower limit of time point for spline basis.}

\item{tU}{upper limit of time point for spline basis.}

\item{nbreak}{number of knots for spline basis.}

\item{p01}{number of covariates in the first biomarker.}

\item{p02}{number of covariates in the second biomarker.}

\item{q_eta}{number of survival covariates.}

\item{t_max}{maximum number of repeated measurements allowed for each subject.}

\item{k_max}{dimension of random effects.}

\item{distr}{distributional assumption of random effects.}

\item{m_age}{true mean value of age variable in the longitudinal sub-model.}

\item{std_age}{true standard deviation of age variable in the longitudinal sub-model.}

\item{tbtheta}{true btheta matrix.}

\item{tbeta0}{true value of beta0 matrix.}

\item{tbeta1}{true value of beta1 vector.}

\item{tsigma}{true value of sigma square vector.}

\item{tsigmad}{true value of variance of random effects.}

\item{ttheta}{true value of theta vector.}

\item{teta}{true value of survival fixed effects.}

\item{tgamma}{true value of latent association parameter.}

\item{lambda0}{true baseline hazard.}

\item{sigma_inv}{sigma_inv.}
}
\description{
Data generation for simulation
}
