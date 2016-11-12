# load the needed routines
source('routines.R')

# set some experiment parameters
sampling.time <- 3           # s
sampling.frequency <- 20000  # Hz

# load the data
dataset <- read.csv('pedranz.csv')$x

# solve the assignment
solution <- SolveAssignment1(sampling.time, sampling.frequency, dataset)
print(solution)
