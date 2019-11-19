\name{Goodness of Fit : Mean Relative Approximation Error}

\alias{gofMRAE}

\title{
Mean Relative Approximation Error
}

\description{
Calculates and returns mean relative approximation error (MRAE).
}

\usage{
gofMRAE(Obs, Prd, dgt = 3)
}

\arguments{
  \item{Obs}{
Observed values or target vector.
}
  \item{Prd}{
  Predicted values. Values produced by approximation or regression.
}
  \item{dgt}{
Number of digits in decimal places. Default is 3.
}
}

\value{
 \item{MeanRelativeApproximationError}{Goodness of fit - mean relative approximation error (MRAE)}
}

\references{
The Connection Dependent Threshold Model for Finite Sources -A Generalization of the Engset Multirate Loss Model - Ioannis D. Moscholios and Michael D. Logothetis.

Competitive adsorption equilibrium modeling of volatile organic compound (VOC) and water vapor onto activated carbon - Imranul I. Laskara, Zaher Hashishoa,⁎, John H. Phillipsb, James E. Andersonc, Mark Nichols.

A new decision tree based algorithm for prediction of hydrogen sulfide solubility in various ionic liquids - Reza Soleimani, Amir Hossein Saeedi Dehaghani, Alireza Bahadori.
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe
}

\examples{
# dummy inputs, independent variable
# integers from 0 to 19
inputs <- 0:19

# dummy targets/observed values, dependent variable
# a product of 2 times inputs minus 5 with some normal noise
targets <- -5 + inputs*1.2 + rnorm(20)

# linear regression model
model<-lm(targets~inputs)

# About the model
summary(model)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit : mean relative approximation error (MRAE)
gofMRAE(targets, predicted)
}
