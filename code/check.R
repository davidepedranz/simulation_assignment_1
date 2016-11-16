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

# real data
real.mean <- 3
real.frequency <- 200
real.amplitude <- 1.25
real.s <- 4.5
real.sigma <- 1

# generate the data
dataset <- SimulateData(2500, sampling.time, sampling.frequency,
                        real.mean, real.frequency, real.amplitude, real.s, real.sigma)

# solve the simulation
solution <- SolveAssignment1(sampling.time, sampling.frequency, dataset)
print(solution)

# print errors
errors <- list(
  mean =  abs(real.mean - solution$mean) / real.mean * 100, 
  frequency = abs(real.frequency - solution$frequency) / real.frequency * 100,
  amplitude = abs(real.amplitude - solution$amplitude) / real.amplitude * 100,
  s = abs(real.s - solution$s) / real.s * 100,
  sigma = abs(real.sigma - solution$sigma) / real.sigma * 100
)
print(errors)
