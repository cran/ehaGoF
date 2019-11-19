\name{Goodness of Fit - Coefficient of Variation}

\alias{gofCV}

\title{
Coefficient of Variation.
}

\description{
Calculates and returns goodness of fit - coefficient of variation (CV).
}

\usage{
gofCV(Obs, Prd, dgt=3)
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
\item{CoeficientOfVariation}{Goodness of fit - coefficient of variation (CV).}
}

\references{
Comparison of Different Data Mining Algorithms for Prediction of Body Weight From Several Morphological Measurements in Dogs - S Celik, O Yilmaz
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

# Goodness of fit - coefficient of variation (CV)
gofCV(targets, predicted)
}

