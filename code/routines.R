# for kurtosis estimation
library('e1071')

ComputeAutocorrelation <- function(dataset) {
  # Computes the autocorrelation function for a given dataset.
  # We condider 3/4 of the result with a lag ranging
  # from 1 to the whole dataset.
  # This function uses the native R autocorrelation function,
  # but apply a correction to denormalize it.
  # 
  # Args:
  #   dataset: The dataset of samples.
  #
  # Returns:
  #   The autocorrelation of the dataset.
  
  # compute length of the dataset
  n <- length(dataset)
  
  # compute the "R autocorrelation"
  autocorrelation <- acf(dataset, lag.max=n, type="covariance", plot = FALSE)
  
  # apply the correction
  # take the first 75% of the result
  last.value <- floor(n * 0.75)
  return((autocorrelation$acf[2:(n+1)] / n:1 * n)[1:last.value])
}

EstimateMean <- function(dataset) {
  # Computes the estimated mean for a dataset of samples taken
  # from a population with the following model:
  #
  #    x(t) = A*sin(2*pi*f*t) + Y + Z
  #    Y ~ Logistic(mu, s), Z ~ Gaussian(0, sigma)
  #
  # Args:
  #   dataset: The dataset of samples.
  #
  # Returns:
  #   The estimated mean of the population.

  return(mean(dataset))
}

EstimateFrequency <- function(dataset, lag.range, sampling.frequency) {
  # Computes the estimated frequency for a dataset
  # of samples taken from a population with the following model:
  #
  #    x(t) = A*sin(2*pi*f*t) + Y + Z
  #    Y ~ Logistic(0, s), Z ~ Gaussian(0, sigma)
  #
  # Args:
  #   dataset:            The dataset of samples.
  #   lag.range:          The range of possible lags.
  #   sampling.frequency: The frequency of the sampling.
  #
  # Returns:
  #   The estimated frequency of the population.

  # autocorrelation of the dataset
  # this should be a cosine function with the same frequency of
  # the original sine and some residual noise due to the variances of Y and Z
  autocorrelation <- ComputeAutocorrelation(dataset)

  # possible lags to test
  lag.possibiles <- seq(lag.range[1], lag.range[2], 1)

  # compute the best lag from the given ones
  peaks.means <- sapply(lag.possibiles, function(lag) {
    peaks.indexes <- seq(1, length(autocorrelation), lag * 2)
    return(mean(abs(autocorrelation[peaks.indexes])))
  })

  # extract the best lag
  # the '-1' is because R is 1-indexed
  lag.best <- lag.possibiles[1] + which.max(peaks.means) - 1

  # we consider for lag the number of measures
  # between a peak and the successive zero of the signal
  # thus, each period contains 4 lags
  return(sampling.frequency / (lag.best * 4))
}

EstimateAmplitude <- function(dataset, estimated.frequency, sampling.frequency) {
  # Computes the estimated amplitude for a dataset
  # of samples taken from a population with the following model:
  #
  #    x(t) = A*sin(2*pi*f*t) + Y + Z
  #    Y ~ Logistic(0, s), Z ~ Gaussian(0, sigma)
  #
  # Args:
  #   dataset:            The dataset of samples.
  #   lag.range:          The estimated frequency of the population.
  #   sampling.frequency: The frequency of the sampling.
  #
  # Returns:
  #   The estimated amplitude of the population.

  # autocorrelation of the dataset
  # this should be a cosine function with the same frequency of
  # the original sine and some residual noise due to the variances of Y and Z
  autocorrelation <- ComputeAutocorrelation(dataset)

  # compute the number of samples between 2 following peaks
  # (positive and negative) of the autocorrelation function
  # this is equal to half of the period of the sinusoide
  half.period <- sampling.frequency / estimated.frequency / 2

  # computes the mean over all peaks in the autocorrelation sinusoide
  # NB: we consider the absolute value of each peak
  #     because we count also the negative peaks
  peaks.indexes <- seq(1, length(autocorrelation), half.period)
  autocorrelation.amplitude <- mean(abs(autocorrelation[peaks.indexes]))

  # the cosine is correlated to the original sine in this way:
  #    x(t) = A * sin(2 * pi * t * f)
  #    R(x(t)) = 1/2 * A^2 * cos(2 * pi * t * f)
  # use this relation to extract the estimated amplitude
  return(sqrt(2 * autocorrelation.amplitude))
}

EstimateVariances <- function(dataset) {
  # Computes the variances for a dataset of samples
  # taken from a population with the following model:
  #
  #    X = Y + Z
  #    Y ~ Logistic(0, s), Z ~ Gaussian(0, sigma)
  #
  # Args:
  #   dataset:            The dataset of samples.
  #
  # Returns:
  #   A list of the estimated variances of Y and Z.

  # notation:
  # var.y = variance of Y
  # var.z = variance of Z
  # var.x = variance of the dataset X = Y + Z
  # kur.x = excess kurtosis of X

  # we can build the following system of equations, where Var[X], Kur[X] - 3
  # are computed from the dataset

  # Y and Z are independent, so Var[X] = Var[Y] + Var[Z]
  # Y and Z are independent, so Kur[X] - 3 = (1.2 * sd.y^4) / (sd.x^2 + sd.y^2)

  # compute the variance Var[X] of the dataset
  var.x <- var(dataset)

  # compute the excess kurtosis (Kur[X] - 3) of the dataset
  kur.x <- kurtosis(dataset)

  # now we solve the system of equations
  # the system admits 2 solutions, but we discard the negative
  # value for the variances, since it does not make any sense
  var.y <- sqrt((kur.x * var.x^2) / 1.2)
  var.z <- var.x - var.y

  # return the list of variances
  return(list(y = var.y, z = var.z))
}

MeanConfidenceIntervalBatchMeans <- function(dataset, batch.size, confidence=0.99) {
  # Computes the confidence interval of the mean mu for a dataset
  # of samples taken from a population with the following model:
  #
  #    x(t) = A*sin(2*pi*f*t) + Y + Z
  #    Y ~ Logistic(mu, s), Z ~ Gaussian(0, sigma)
  #
  # Args:
  #   dataset:    The dataset of samples.
  #   batch.size: The number of samples to use for each batch.
  #   confidence: The required confidence.
  #
  # Returns:
  #   Epsilon of the simmetric confidence interval obtained
  #   using the batch means method.

  # compute the batches
  batches <- split(dataset, ceiling(seq_along(dataset) / batch.size))

  # take the mean of each batch
  # now we should have i.i.d. values and
  # we can treat them as samples from a Gaussian population
  batches.means <- sapply(batches, mean)

  # we take a simmetric confidence interval
  alpha <- qnorm(1 - ((1 - confidence) / 2))

  # compute the needed statistics on batches.means
  batches.sd <- sd(batches.means)
  batches.length <- length(batches.means)

  # compute the epsilon
  return(alpha * batches.sd / sqrt(batches.length))
}

SolveAssignment1 <- function(sampling.time, sampling.frequency, dataset.raw) {
  
  # compute some useful info about the experiment
  sampling.points <- sampling.time * sampling.frequency                      # points
  sampling.times <- (1:sampling.points) * (sampling.time / sampling.points)  # s
  
  #####################################
  # estimate the parameters
  #####################################
  
  # 1st parameter = mean
  estimated.mean <- EstimateMean(dataset.raw)
  
  # make the population zero-mean
  # this centers the sinusioide on the y-axis
  # and makes easier to estimate its frequency
  dataset.zeromean <- dataset.raw - estimated.mean
  
  # 2nd parameter = frequency
  estimated.frequency <- EstimateFrequency(dataset.zeromean, c(20, 30), sampling.frequency)
  
  # 3rd parameter = amplitude
  estimated.amplitude <- EstimateAmplitude(dataset.zeromean, estimated.frequency, sampling.frequency)
  
  # remove the sinusoide from the dataset
  dataset.nosinusoide <- dataset.zeromean - (estimated.amplitude * sin(2 * pi * estimated.frequency * sampling.times))
  
  # 4th, 5th parameters = s and sigma
  estimated.variances <- EstimateVariances(dataset.nosinusoide)
  estimated.s <- sqrt(estimated.variances$y * 3) / pi
  estimated.sigma <- sqrt(estimated.variances$z)
  
  
  #####################################
  # confidence interval for mean
  #####################################
  
  # using the batch means, taking 2 wavelenghts as batch size
  batch.size <- sampling.frequency / estimated.frequency * 2
  raw.confidence.95 <- MeanConfidenceIntervalBatchMeans(dataset.raw, batch.size, 0.95)
  raw.confidence.99 <- MeanConfidenceIntervalBatchMeans(dataset.raw, batch.size, 0.99)
  
  # working on the data without the sinusoide
  # NB: it makes no difference for the computation to
  #     remove or not the estimated mean
  # NB: we test also batches of 1 unit,
  #     which means treating the original population as Gaussian
  cleaned.confidences.95 <- sapply(1:200, function(size) MeanConfidenceIntervalBatchMeans(dataset.nosinusoide, size, 0.95))
  cleaned.confidence.95 <- min(cleaned.confidences.95)
  cleaned.confidences.99 <- sapply(1:200, function(size) MeanConfidenceIntervalBatchMeans(dataset.nosinusoide, size, 0.99))
  cleaned.confidence.99 <- min(cleaned.confidences.99)
  
  
  #####################################
  # results
  #####################################
  return(list(
    mean=estimated.mean,
    frequency=estimated.frequency,
    amplitude=estimated.amplitude,
    s=estimated.s,
    sigma=estimated.sigma,
    confidence.95=raw.confidence.95,
    confidence.99=raw.confidence.99,
    confidence2.95=cleaned.confidence.95,
    confidence2.99=cleaned.confidence.99
  ))
}

