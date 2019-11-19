\name{Goodness of Fit - Root Mean Square Error}

\alias{gofRMSE}

\title{
Root Mean Square Error
}

\description{
Calculates and returns root mean square error (RMSE).
}

\usage{
gofRMSE(Obs, Prd, dgt = 3)
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
\item{RootMeanSquareError}{Root mean square error (RMSE)}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan - Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas, Wilhelm Grzesiak

Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019

Indirect Estimation of Structural Parameters in South African Forests Using MISR-HR and LiDAR Remote Sensing Data - Precious Nokuthula Wistebaar Mahlangu, Renaud Mathieu, Konrad Wessels, Laven Naidoo, Michel M Verstraete, Gregory P Asner, Russell Main, Remote Sens. 2018, 10, 1537 ; doi:10.3390/rs10101537.
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

# Goodness of fit - root mean square error (RMSE)
gofRMSE(targets, predicted)
}

