\name{Goodness of Fit : Akaike's Information Criterion}

\alias{gofAIC}

\title{
Akaike's Information Criterion
}

\description{
Calculates and returns Akaike's information criterion (AIC).
}

\usage{
gofAIC(Obs, Prd, nTermInAppr=2, dgt = 3)
}

\arguments{
  \item{Obs}{
  Observed or measured values or target vector.
}
  \item{Prd}{
Predicted or fitted values by the model. Values produced by approximation or regression.
}
  \item{nTermInAppr}{
  Number of terms in approximation or regression models formula, including interception. For simple linear regression with one independent variable is simply 2. Default is 2.
}
  \item{dgt}{
Number of digits in decimal places. Default is 3.
}
}

\value{
 \item{AkaikesInformationCriterion }{Akaike's information criterion (AIC)}
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
When \code{n/k}  is not greater than 40, where \code{n} is the number of observations and \code{k} is the number of terms in approximation, Corrected Akaike's Information Criterion (\code{gofCAIC}) is used.
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

# Number of Terms
n = length(model$coefficients)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit : Akaike's information criterion (AIC)
gofAIC(targets, predicted, dgt=4, nTermInAppr=n)
}
