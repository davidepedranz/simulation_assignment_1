# simulate assignment
SimulateData <- function(seed, sampling.time=3, sampling.frequency=20000, 
                         mean=3, frequency=200, amplitude=1.25, s=4.5, sigma=1) {
  
  # set the seed 
  set.seed(seed)
  
  # compute some useful info about the experiment
  sampling.points <- sampling.time * sampling.frequency                      # points
  sampling.times <- (1:sampling.points) * (sampling.time / sampling.points)  # s
  
  # compute the sinusoide
  sinusoide <- amplitude * sin(2 * pi * frequency * sampling.times)
  
  # generate the logistic
  y <- rlogis(sampling.points, location=mean, scale=s)
  
  # generate the gaussian
  z <- rnorm(sampling.points, mean=0, sd=sigma)
  
  # compute dataset
  return(sinusoide + y + z)
}

# load the needed routines
source('routines.R')

# set some experiment parameters
sampling.time <- 3           # s
sampling.frequency <- 20000  # Hz

# generate the data
dataset <- SimulateData(1012, sampling.time, sampling.frequency)

# solve the simulation
solution <- SolveAssignment1(sampling.time, sampling.frequency, dataset)
print(solution)
