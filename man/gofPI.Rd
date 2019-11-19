\name{Goodness of Fit : Performance Index}

\alias{gofPI}

\title{
Performance Index
}

\description{
Calculates and returns performance index (PI).
}

\usage{
gofPI(Obs, Prd, dgt = 3)
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
\item{PerformanceIndex}{Goodness of fit: performance index (PI)}
}

\references{
Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019
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

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit - performance index (PI)
gofPI(targets, predicted)
}

