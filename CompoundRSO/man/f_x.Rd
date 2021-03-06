% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_X_matrix_v2.0.R
\name{f_x}
\alias{f_x}
\title{function to create a row of the design matrix of a specified model order}
\usage{
f_x(
  x,
  model_order,
  C_dtype,
  me_index = c(),
  qe_index = c(),
  two_fi_index = c(),
  cubic_index = c(),
  quatric_index = c()
)
}
\arguments{
\item{x}{a vector or row matrix of factor levels to generate the design matrix for}

\item{model_order}{a string indicating the order of the regression model. "linear","quadratic" (Main Effects+Quadratic Effects), "2FI" (Main Effects + TwoFactor Interactions),"full", "non_standard"- for a user defined model order order}

\item{C_dtype}{a string vector indicating the data type of the experimental factors in the design either "cont"(continuous) or "cat"(categorical)}

\item{me_index}{vector of index of Main Effects to be included in the design for model_order "non_standard"}

\item{qe_index}{vector index of Quadratic Effects to be included in the design for model_order "non_standard"}

\item{two_fi_index}{list of index of Two Factor Interactions to be included in the design for model_order "non_standard"}

\item{cubic_index}{vector of index of Cubic Effects to be included in the design for model_order "non_standard" (although a cubic design can be sampled, current design algorithms don't support cubic models)}

\item{quatric_index}{vector of index of Quatric Effects to be included in the design for model_order "non_standard" (although a quatric design can be sampled, current design algorithms don't support cubic models)}
}
\value{
a row of the design matrix for the specified factor levels and  model_order
 @examples
  f_x(x,model_order="full",C_dtype, me_index=c(),qe_index=c(),two_fi_index=c(),cubic_index=c(),quatric_index=c())
}
\description{
function to create a row of the design matrix of a specified model order
}
