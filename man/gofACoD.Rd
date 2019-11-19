\name{Goodness of Fit : Adjusted Coefficient of Determination (Adjusted R-Squared)}

\alias{gofACoD}

\title{
Adjusted Coefficient of Determination (Adjusted R-Squared)
}

\description{
Calculates and returns adjusted coefficient of determination (adjusted R-squared).
}

\usage{
gofACoD(Obs, Prd, nTermInAppr = 2, dgt = 3)
}

\arguments{
  \item{Obs}{
  Observed or measured values or target vector.
}
  \item{Prd}{
Predicted or fitted values by the model. Values produced by approximation or regression.
}
  \item{nTermInAppr}{
Number of terms in approximation or regression models formula, interception included. For simple linear regression with one independent variable is simply 2. Default is 2.
}
  \item{dgt}{
Number of digits in decimal places. Default is 3.
}
}

\value{
\item{AdjustedCoefficientofDetermination}{Goodness of fit - adjusted coefficient of determination (adjusted R-squared)}
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
# integers from 0 to 19
inputs <- 0:19

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*1.2 + rnorm(20)

# linear regression model
model<-lm(targets~inputs)

# About the model
summary(model)

# Number of Terms
n = length(model$coefficients)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit : adjusted coefficient of determination (adjusted R-squared)
gofACoD(targets, predicted, dgt=4,nTermInAppr=n)
}
