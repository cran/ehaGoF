\name{Goodness of fit - Relative Root Mean Square Error}

\alias{gofRRMSE}

\title{
Relative Root Mean Square Error
}

\description{
Calculates and returns relative root mean square error (RRMSE) of the model. The ratio of the mean of square root of residuals squared to the mean of observed values.
}

\usage{
gofRRMSE(Obs, Prd, dgt = 3)
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

\details{
RRMSE is calculated by dividing RMSE by the mean of observed values.
}

\value{
\item{RelativeRootMeanSquareError}{Relative root mean square error (RRMSE) of given set.}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan, Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas and Wilhelm Grzesiak, Pakistan J. Zool., vol. 49(1), pp 257-265, 2017.

Indirect Estimation of Structural Parameters in South African Forests Using MISR-HR and LiDAR Remote Sensing Data - Precious Nokuthula Wistebaar Mahlangu, Renaud Mathieu, Konrad Wessels, Laven Naidoo, Michel M Verstraete, Gregory P Asner, Russell Main
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe
}

\note{
For more information:
Despotovic, M., Nedic, V., Despotovic, D., & Cvetanovic, S., Evaluation of empirical models for predicting monthly mean horizontal diffuse solar radiation, Renewable and Sustainable Energy Reviews, Volume 56, April 2016, Pages 246-260, ISSN 1364-0321, http://dx.doi.org/10.1016/j.rser.2015.11.058.
}

\examples{
# Input values, independent variable
input <- 0:4

# Target vector, observed values, dependent variable
target <- c(1.9, 4.1, 5.89, 7.9, 10.01)

# Simple linear regression, target across input like: target = a * input + b,
# where a and b are coefficients.
model <- lm(target~input)

# Information about the model
summary(model)

# Values predicted by the model
predicted <- predict(model)

# using library ehaGoF for goodness of fit
library(ehaGoF)

# Goodness of fit - relative root mean square error (RRMSE)
gofRRMSE(target, predicted)
}

