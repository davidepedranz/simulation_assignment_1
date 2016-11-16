source('routines.R')

# load the data
dataset <- read.csv('./pedranz.csv')$x

# remove the mean
estimated.mean <- EstimateMean(dataset)
dataset.zeromean <- dataset - estimated.mean

# compute the autocorrelation
autocorrelation <- ComputeAutocorrelation(dataset.zeromean)

# plot
source('plots.R')
PlotAutocorrelation(autocorrelation[1:400])
