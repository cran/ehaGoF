\name{Goodness of Fit : Coefficient of Determination (R-Squared)}

\alias{gofCoD}

\title{
Coefficient of Determination (R-Squared)
}

\description{
Calculates and returns coefficient of determination (R-squared).
}

\usage{
gofCoD(Obs, Prd, dgt = 3)
}

\arguments{
  \item{Obs}{
  Observed or measured values or target vector.
}
  \item{Prd}{
Predicted or fitted values by the model. Values produced by approximation or regression.
}
  \item{dgt}{
Number of digits in decimal places. Default is 3.
}
}

\value{
\item{CoefficientofDetermination}{Goodness of fit - coefficient of determination (R-squared)}
}

\references{
Comparison of Different Data Mining Algorithms for Prediction of Body Weight From Several Morphological Measurements in Dogs - S Celik, O Yilmaz.

A new decision tree based algorithm for prediction of hydrogen sulfide solubility in various ionic liquids - Reza Soleimani, Amir Hossein Saeedi Dehaghani, Alireza Bahadori.
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe
}

\examples{
# dummy inputs, independent variable
# integers from 0 to 99
inputs <- 0:99

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*1.2 + rnorm(100)

# linear regression model
model<-lm(targets~inputs)

# About the model
summary(model)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit : coefficient of determination (R-squared)
gofCoD(targets, predicted)
}

