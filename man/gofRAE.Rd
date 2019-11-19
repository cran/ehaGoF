\name{Goodness of Fit : Global Relative Approximation Error}

\alias{gofRAE}

\title{
Global Relative Approximation Error
}

\description{
Calculates and returns global relative approximation error (RAE).
}

\usage{
gofRAE(Obs, Prd, dgt = 3)
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
\item{RelativeApproximationError}{Global relative approximation error (RAE)}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan - Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas, Wilhelm Grzesiak.

Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019.

Comparison of Different Data Mining Algorithms for Prediction of Body Weight From Several Morphological Measurements in Dogs - S Celik, O Yilmaz.

The Connection Dependent Threshold Model for Finite Sources -A Generalization of the Engset Multirate Loss Model - Ioannis D. Moscholios and Michael D. Logothetis.
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

# Goodness of fit : global relative approximation error (RAE)
gofRAE(targets, predicted)
}

