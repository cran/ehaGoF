\name{Goodness of Fit - Mean Absolute Percentage Error}

\alias{gofMAPE}

\title{
Mean Absolute Percentage Error
}

\description{
Calculates and returns the goodness of fit criterion: mean absolute percentage error (MAPE), a.k.a. mean absolute percentage deviation or MAPD.
}

\usage{
gofMAPE(Obs, Prd, dgt=3)
}

\arguments{
  \item{Obs}{
  Observed values or targets.
}
\item{Prd}{
Predicted or expre-ected values produced by the model.
}
\item{dgt}{
Number of digits in decimal places. Default is 3.
}
}

\details{
Mean absolute percentage error (MAPE) is a measure of prediction accuracy of a forecasting method in statistics. It is commonly used as a loss function for regression problems and in model evaluation, for its very intuitive interpretation in terms of relative error. It usually expresses accuracy as a percentage.
}

\value{
\item{MAPE}{Mean absolute percentage error (MAPE) of given set.}
}

\references{
Comparison of the Predictive Capabilities of Several Data Mining Algorithms and Multiple Linear Regression in the Prediction of Body Weight by Means of Body Measurements in the Indigenous Beetal Goat of Pakistan, Ecevit Eyduran, Daniel Zaborski, Abdul Waheed, Senol Celik, Koksal Karadas and Wilhelm Grzesiak, Pakistan J. Zool., vol. 49(1), pp 257-265, 2017.

Prediction of Selected Reproductive Traits of Indigenous Harnai Sheep under the Farm Management System via various Data Mining Algorithms - Daniel Zaborski, Muhammad Ali, Ecevit Eyduran, Wilhelm Grzesiak, Mohammad Masood Tariq, Ferhat Abbas, Abdul Waheed, Cem Tirink - Pakistan journal of zoology, 2019

Comparison of Different Data Mining Algorithms for Prediction of Body Weight From Several Morphological Measurements in Dogs - S Celik, O Yilmaz
}

\author{
Prof. Dr. Ecevit Eyduran,
TA. Alper Gulbe
}

\note{
For more information look at these papers:

Rob J. Hyndman, Anne B. Koehler,
Another look at measures of forecast accuracy,
International Journal of Forecasting,
Volume 22, Issue 4,
2006,
Pages 679-688,
ISSN 0169-2070,


Arnaud de Myttenaere, Boris Golden, Bénédicte Le Grand, Fabrice Rossi,
Mean Absolute Percentage Error for regression models,
Neurocomputing,
Volume 192,
2016,
Pages 38-48,
ISSN 0925-2312,
}

\examples{
# dummy inputs, independent variable
# integers from 0 to 9
inputs <- 0:9

# dummy targets/observed values, dependent variable
# a product of 2*times inputs minus 5 with some normal noise
targets <- -5 + inputs*2 + rnorm(10)

# linear regression model
model<-lm(targets~inputs)

# model's predicted values against targets
predicted<-model$fitted.values

# using library ehaGoF for goodness of fit.
library(ehaGoF)

# Goodness of fit - mean absolute percentage error statistics
gofMAPE(targets, predicted)
}

