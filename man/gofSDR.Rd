\name{Goodness of Fit - Standard Deviation Ratio}

\alias{gofSDR}

\title{
Standard Deviation Ratio
}

\description{
Calculates and returns standard deviation ratio (SDR).
}

\usage{
gofSDR(Obs, Prd, dgt=3)
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
\item{StandardDeviationRatio}{Standard deviation ratio (SDR) of given set.}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan, Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas and Wilhelm Grzesiak, Pakistan J. Zool., vol. 49(1), pp 257-265, 2017.

Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019

Comparison of Different Data Mining Algorithms for Prediction of Body Weight From Several Morphological Measurements in Dogs - S Celik, O Yilmaz
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gülbe
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

# Input values, independent variable
input <- 0:4

# Target vector, observed values
target <- c(1.9, 4.1, 5.89, 7.9, 10.01)

# Simple linear regression target across input like target = a * input + b,
# where a and b are coefficients.
model <- lm(target~input)

# Information about the model
summary(model)

# Values predicted by the model
predicted <- predict(model)

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit - standard deviation ratio (SDR)
gofSDR(target, predicted)
}


