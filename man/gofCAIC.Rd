\name{Goodness of Fit : Corrected Akaike's Information Criterion}

\alias{gofCAIC}

\title{
Corrected Akaike's Information Criterion
}

\description{
Calculates and returns corrected Akaike's information criterion.
}

\usage{
gofCAIC(Obs, Prd, nTermInAppr = 2, dgt = 3)
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
  \item{CorrectedAkaikesInformationCriterion}{Goodness of fit - corrected Akaike's information criterion (cAIC)}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan - Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas, Wilhelm Grzesiak.

Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019.
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe
}

\note{
When \code{n/k}  is greater than 40, where \code{n} is the umber of observations and \code{k} is the number of terms in approximation, Akaike's Information Criterion (\code{gofAIC}) is used.
}

\examples{
# dummy inputs, independent variable
# integers fron 0 to 79
inputs <- 0:79

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*1.2 + rnorm(80)

# linear regression model
model<-lm(targets~inputs)

# About the model
summary(model)

# Number of Terms
n = length(model$coefficients)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for Goodness of Fit function (GoF)
library(ehaGoF)

# Goodness of Fit : Corrected Akaike's Information Criterion
gofCAIC(targets, predicted, dgt=4, nTermInAppr=n)
}
