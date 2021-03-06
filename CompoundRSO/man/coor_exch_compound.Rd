% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{coor_exch_compound}
\alias{coor_exch_compound}
\title{function for constructing the Adaptive-RSO designs using coordinate exchange}
\usage{
coor_exch_compound(
  x,
  x_matrix,
  C,
  model_order,
  C_dtype,
  W,
  freeze_rows,
  x_matrix_Dopt,
  x_matrix_Iopt,
  w,
  K,
  me_index,
  qe_index,
  two_fi_index,
  me_index_daug,
  qe_index_daug,
  two_fi_index_daug,
  D_opt_type = "Bayesian",
  telescoping = FALSE
)
}
\arguments{
\item{x}{settings matrix of the design}

\item{x_matrix}{design matrix of the specified order including the intercept}

\item{C}{list of candidate points for all factors in the design}

\item{model_order}{string indicating linear,quadratic (Main Effects+Quadratic Effects), 2FI (Main Effects + TwoFactor Interactions),full, non_standard for a user defined order}

\item{C_dtype}{data type of the factors in the design either "cont"(continuous) or "cat"(categorical)}

\item{W}{Moment matrix computed for a given lambda}

\item{freeze_rows}{the rows of the matrix to freeze while augmenting, 0 when type is "non_augment"}

\item{x_matrix_Dopt}{the D-optimal design matrix}

\item{x_matrix_Iopt}{the I-optimal design matrix}

\item{w}{a numeric value of weight that is used to compute the compound optimal design}

\item{K}{a diagonal matrix specifying the prior variance in Bayesian D-optimal design set to 0.001 for potential terms and 0 for primary terms}

\item{me_index}{vector of index of Main Effects to be included in the I-optimal design for model_order "non_standard"}

\item{qe_index}{vector index of Quadratic Effects to be included in the I-optimal design for model_order "non_standard"}

\item{two_fi_index}{list of index of Two Factor Interactions to be included in the I-optimal design for model_order "non_standard"}

\item{me_index_daug}{vector of index of Main Effects to be included in the D-optimal design for model_order "non_standard"}

\item{qe_index_daug}{vector index of Quadratic Effects to be included in the D-optimal design for model_order "non_standard"}

\item{two_fi_index_daug}{list of index of Two Factor Interactions to be included in the D-optimal design for model_order "non_standard"}

\item{D_opt_type}{Type of D-optimal design Bayesian or regular D-optimal design}

\item{telescoping}{a boolean value taking TRUE or FALSE. Telescoping reduces the candidate factors the local region of the optimum design
\itemize{
      \item x_matrix- the design matrix for a given model order
      \item opt_val- the compound optimal objective value
      \item x_best- matrix of factor levels for the optimal design
      \item deff_best- the d-efficiency of the optimal design
      \item ieff_best- the uniform efficiency of the optimal design
      }}

\item{type}{"non_augment" or "augment" (augment an existing deisgn)}
}
\description{
function for constructing the Adaptove-RSO designs using coordinate exchange. The Adaptive-RSO designs is the compound optimal design between Bayesian D-optimal design and I-optimal design
}
\examples{
\dontrun{
coor_exch_compound(x,x_matrix,C,model_order,me_index = c(),qe_index = c(),two_fi_index = c(),C_dtype, W, freeze_rows=0,x_matrix_Dopt,x_matrix_Iopt,w,K=0,D_opt_type="Bayesian",me_index_daug, qe_index_daug,two_fi_index_daug,telescoping=FALSE)
}
}
