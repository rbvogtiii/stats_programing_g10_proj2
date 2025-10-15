## Work was evenly distributed between all group members
## First, we each created our own solution to the project
## Then we compiled our versions and kept the best/most efficient code
## We worked together on the compiled version to make edits, add comments, and optimise

## Code to simulate the SEIR model of disease spread:
## SEIR stands for Susceptible / Exposed / Infectious / Recovered.
## We assume that a given part of the population is infectious,
## and the rest of the population is susceptible, on the first day.
## We then simulate the movements of the population between the four categories.
## People may be infected by people in their household, 
## who they see the most, with a given daily probability.
## There is another set daily probability that they may be infected
## by a member of their regular network or contacts (which is also simulated).
## Finally, they may be infected by a random person with a daily probability
## which is dependent on their "sociability."
## As well as the above probabilities of moving from state S to E,
## we model the movement from state E to I with a daily probability of gamma,
## and the movement from state I to R with a daily probability of delta.


get_h <- function(n, hmax = 5) {
  ## given a population size of n, the population is sorted into households
  ## hmax is the cap on household sizes, with size uniformly distributed from 1 to hmax
  sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])
}


net_helper <- function(idx, probs) {
  ## given a person (idx) and an adjacency matrix (probs),
  ## return a list of that person's connections
  connections <- which(probs[idx, ] == 1 | probs[, idx] == 1) ## find connections
  if (length(connections) == 0) { ## no connections
    NA
  } else {
    connections ## return vector of connections
  }
}

get.net <- function(beta, h, nc = 15) {
  ## Given a vector containing a sociability parameter for each member of
  ## the population (beta), a vector containing each person's household (h),
  ## and the average number of connections each person should have, return a
  ## list where the element corresponding to a given person contains their
  ## network. The probablility that two people are connected is proportional
  ## to their sociability parameter.
  n <- length(beta)
  link_const <- nc / ((mean(beta)**2) * (n - 1))
  probs <- outer(beta, beta) * link_const ## initialize probabilities for all pairings
  for (i in 1:n) probs[i, i:n] <- 0 ## get rid of double counting
  probs[outer(h, h, FUN = "==")] <- 0 ## remove connections within households
  rands <- runif(n * n)
  probs[probs >= rands] <- 1
  lapply(seq_along(beta), net_helper, probs) ## create list containing network
}


nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
  # Simulate the movement of people between susceptible, exposed,
  # infected, and recovered groups, over a specified
  # number of days. Infected people recover with probability
  # delta and exposed people become infected with probability gamma.
  # Susceptible people can be exposed through their household with probability,
  # alpha[1], through their social network with probability alpha[2], or randomly
  # with a probability proportional to the product of their sociability and the
  # sociability of each currently infected person.
  
  
  # Initialization of data 
  n <- length(beta)
  pop <- 1:n
  
  # Randomly infect a small proportion of the population to start the simulation.
  I <- sample(pop, n * pinf) 
  S <- pop[!(pop %in% I)] # put rest of population in S
  E <- c()
  R <- c()
  
  # Identify and compute a scaling constant for random mixing infections. 
  infection_const <- (alpha[3] * nc) / ((mean(beta)**2) * (n - 1))
  
  # nNitialize storgae vectors for counts
  S_out <- c(length(S))
  E_out <- c(length(E))
  I_out <- c(length(I))
  R_out <- c(length(R))
  sum_out <- c(length(S) + length(E) + length(I) + length(R))
  t <- c(0)
  
  for (day in 1:nt) {
    # Moving from E to I Transition
    # Each exposed person becomes infectious with daily probability delta.
    I_prob <- gamma - runif(length(E))
    I <- c(I, E[I_prob >= 0])
    E <- E[I_prob < 0]
    
    # Moving from I to R Transition 
    # Each infectious person recovers with daily probability gamma.
    R_prob <- delta - runif(length(I))
    R <- c(R, I[R_prob >= 0])
    I <- I[R_prob < 0]
    
    # If there are no more infectious individuals i.e the epidemic ends, record results and skip to the following day.
    if (length(I) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }
    
    # Household exposures
    # We calculate the number of infected people in each
    # susceptible person's household, using that to find the probability
    # that they are not exposed by any of them, then use that to
    # determine the probability they are exposed by at least one
    E_prob <- 1 - (((1 - alpha[1])**tabulate(h[pop %in% I], nbins = max(h)))[h[S]]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]
    
    # If the epidemic end out after this step, record and continue to the next day.
    if (length(I) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }
    
    # Random exposures
    # We calculate the probability that each infected person will not
    # expose each member of the susceptible group, then use that to
    # determine the probability they are exposed to at least one
    E_prob <- 1 - ((t(matrix(beta[S], nrow = length(S), ncol = length(I))) * beta[I]) * infection_const)
    E_prob <- 1 - apply(E_prob, 2, prod) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]
    
    # If the epidemic ends after this step, record and continue to the next day.
    if (length(I) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }
    
    if (length(unlist(alink[I])) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }
    
    # Nettwork exposures
    # We calculate the number of infected people in each
    # susceptible person's network and use that to find the probability
    # that they are not exposed by any of them. Then, we determine 
    # the probability they are exposed to at least one
    E_prob <- 1 - ((1 - alpha[2])**tabulate(unlist(alink[I]), nbins = n)[S]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]
    
    # Record the daily count/compartment sizes.
    S_out <- c(S_out, length(S))
    E_out <- c(E_out, length(E))
    I_out <- c(I_out, length(I))
    R_out <- c(R_out, length(R))
    sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
    t <- c(t, day)
  }
  out <- list(S_out, E_out, I_out, R_out, sum_out, t)
  names(out) <- c("S", "E", "I", "R", "Sum", "t")
  out
}

plot_nseir <- function(result, n) {
  ## function to plot solutions, 
  ## where "result" is an output from the nseir function
  ## so it contains vectors S,E,I,R,t
  ## we plot each of S,E,I,R against time
  ## n is used in scaling the graph axes and grid

  ## setting up the empty plot with appropriate axis lengths  
  plot(result$t, result$S, type = "n", xlim = range(result$t), ylim = c(0, n),
       xlab = "Day", ylab = "Number of Individuals")
  
  ## adding a grid for ease of reading results
  abline(v = seq(0, length(result$t)-1, 10), lty = 1, col = "gray90")
  abline(h = seq(0, n, n/10), lty = 1, col = "gray90")
  
  ## plotting the number of people in each category over time
  lines(result$t, result$S, col = "lightseagreen", lwd = 3)
  lines(result$t, result$E, col = "goldenrod2", lwd = 3)
  lines(result$t, result$I, col = "orangered3", lwd = 3)
  lines(result$t, result$R, col = "darkgreen", lwd = 3)
  
  legend("left", 
         legend = c("Susceptible", "Exposed", "Infectious", "Recovered"),
         col = c("lightseagreen", "goldenrod2", "orangered3", "darkgreen"),
         lwd = 3, bty = "n", cex=0.6)
}



## run the simulation with select parameters:
n <- 10000 ## population size
nc <- 15 ## average number of contacts each person has
h <- get_h(n)
beta <- runif(n)
alink <- get.net(beta, h, nc)

## scenario 1: full model with default parameters
s1 <- nseir(beta, h, alink)
## scenario 2: removing the household and network structure (random mixing only)
s2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))
## scenario 3: removing the variability of beta
const_beta <- rep(mean(beta), length(beta))
alink <- get.net(const_beta, h, nc)
s3 <- nseir(const_beta, h, alink)
## scenario 4: both random mixing only and constant beta
s4 <- nseir(const_beta, h, alink, alpha = c(0, 0, 0.04))

# plotting scenarios:
par(mfrow = c(2, 2))

plot_nseir(s1,n)
title(main = "Full Model")

plot_nseir(s2,n)
title(main = "Random Mixing Only")

plot_nseir(s3,n)
title(main = "Constant Beta")

plot_nseir(s4,n)
title(main = "Random Mixing & Constant Beta")
