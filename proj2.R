## part 1
get_h <- function(n, hmax = 5) {
  sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])
}

## part 2
net_helper <- function(idx, probs) { # creates a list of connections for the given person
  connections <- which(probs[idx, ] == 1 | probs[, idx] == 1) # find connections
  if (length(connections) == 0) { # no connections
    NA
  } else {
    connections # return vector of connections
  }
}

get.net <- function(beta, h, nc = 15) {
  n <- length(beta)
  link_const <- nc / ((mean(beta)**2) * (n - 1))
  probs <- outer(beta, beta) * link_const # initialize probabilities for all pairings
  for (i in 1:n) probs[i, i:n] <- 0 # get rid of double counting
  probs[outer(h, h, FUN = "==")] <- 0 # remove connections within households
  # rands <- matrix(runif(n * n), n, n)
  rands <- runif(n * n)
  probs[probs >= rands] <- 1
  lapply(seq_along(beta), net_helper, probs) # create list containing network
}

## part 3
nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
  n <- length(beta)
  pop <- 1:n
  I <- sample(pop, n * pinf) # randomly choose pinf% of population to start infected
  S <- pop[!(pop %in% I)] # put rest of population in S
  E <- c()
  R <- c()
  infection_const <- (alpha[3] * nc) / ((mean(beta)**2) * (n - 1))

  S_out <- c(length(S))
  E_out <- c(length(E))
  I_out <- c(length(I))
  R_out <- c(length(R))
  sum_out <- c(length(S) + length(E) + length(I) + length(R))
  t <- c(0)

  for (day in 1:nt) {
    # move from E to I with prob gamma
    I_prob <- gamma - runif(length(E))
    I <- c(I, E[I_prob >= 0])
    E <- E[I_prob < 0]

    # move from I to R with prob delta
    R_prob <- delta - runif(length(I))
    R <- c(R, I[R_prob >= 0])
    I <- I[R_prob < 0]

    if (length(I) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }

    # household exposures
    E_prob <- 1 - (((1 - alpha[1])**tabulate(h[pop %in% I], nbins = max(h)))[h[S]]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]

    if (length(I) == 0) {
      S_out <- c(S_out, length(S))
      E_out <- c(E_out, length(E))
      I_out <- c(I_out, length(I))
      R_out <- c(R_out, length(R))
      sum_out <- c(sum_out, (length(S) + length(E) + length(I) + length(R)))
      t <- c(t, day)
      next
    }

    # random exposures
    E_prob <- 1 - ((t(matrix(beta[S], nrow = length(S), ncol = length(I))) * beta[I]) * infection_const)
    E_prob <- 1 - apply(E_prob, 2, prod) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]

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

    # network exposures
    E_prob <- 1 - ((1 - alpha[2])**tabulate(unlist(alink[I]), nbins = n)[S]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]

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

# DESCRIPTION: Implements an SEIR model incorporating household structure, regular contact networks
# and random mixing. This function simulates the spread of a disease through a population over a specified
# number of days (t).

# input beta An n-vector of sociability weights (βi) per person.
# input h An n-vector of Household IDs.
# input alink A contact list defining the regular (non-household) contacts for each person (from get.net).
# input alpha A vector of 3 infection probabilities: c(α_h) for household, (α_c) for contact, and (α_r) for random mixing.
# input delta The daily transition rate of an Exposed person becoming Infectious (E -> I).
# input gamma The dailytransition rate of an Infectious person Recovering (I -> R).
# input nc The average number of (random mixing) contacts per person per day.
# input nt The number of days to simulate.
# input pinf The fraction/seed of the population to start in the Infectious state.

# output A list containing vedctors S, E, I, R (total daily population counts) and the day number t(i).

nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
  
  # #Initializing conditions of the epidemic. 
  n <- length(beta)
  alpha_h <- alpha[1]; alpha_c <- alpha[2]; alpha_r <- alpha[3] # Involves infecting a small subset at the start.
  
  infection_const <- ((alpha_r * nc) / (mean(beta)^2 * (n - 1)))  # Constant used to scale random-mixing transmission probabilities. 

  status <- rep(1, n) # Define population states using integers for efficiency.

  # Randomly assign the initial infected individuals and ensure at least one person is infected to start the simulation.
  num_initial_infected <- max(1, round(n * pinf))
  initial_infected_idx <- sample(n, num_initial_infected)
  status[initial_infected_idx] <- 3 # Start in state I

  # Create storage vectors to record the number of people in each stage on day 1. 
  S <- E <- I <- R <- numeric(nt)
  counts <- tabulate(status, nbins = 4)

  S[1] <- counts[1]; E[1] <- counts[2]; I[1] <- counts[3]; R[1] <- counts[4]

  # Main Simulation Loop identifying indices of people in each stage
   for (t in 2:nt) {
   # Identify indices of individuals in each state at the start of the day.
   s_idx <- which(status == 1)
   e_idx <- which(status == 2)
   i_idx <- which(status == 3)
   
   #If no person is infected or exposed, the epidemic will have ended.
   if (length(e_idx) == 0 && length(i_idx) == 0) {
     # Fill the remaining days with the last recorded values and break the loop.
     S[t:nt] <- S[t-1]; E[t:nt] <- 0; I[t:nt] <- 0; R[t:nt] <- R[t-1]
     break
   }
   
   # State Transitions: depend on the state at the start of the day.
   
   # Each infectious person recovers (I -> R) with daily probability delta.
   if (length(i_idx) > 0) {
     status[i_idx[runif(length(i_idx)) < delta]] <- 4
   }
   
   # [cite_start]Each exposed person becomes infectious (E -> I) with daily probability gamma. 
   if (length(e_idx) > 0) {
     status[e_idx[runif(length(e_idx)) < gamma]] <- 3
   }
   
   # State S -> E Transition
   # New infections are based on who was infectious at the starting point of the day (i_idx).
   if (length(s_idx) > 0 && length(i_idx) > 0) {
     
     # Initialize probability of avoiding infection for each susceptible person (from all sources)
     p_avoid <- rep(1, length(s_idx))
     
     # a) Household Transmission
     i_per_hh <- tabulate(h[i_idx], nbins = max(h)) # Count the number of infectious individuals in each household.
     k <- i_per_hh[h[s_idx]] # For each susceptible person, find the number of infectious housemates (k).
     p_avoid <- p_avoid * ((1 - alpha_h)^k) # Update avoidance probability.
     
     # b) Regular Contact Network Transmission 
     is_infectious_bool <- status == 3 
     m <- sapply(alink[s_idx], function(contacts) sum(is_infectious_bool[contacts]))  # For each susceptible person, count their infected contacts (m).
     p_avoid <- p_avoid * ((1 - alpha_c)^m)  # Update avoidance probability:
     
     # c) Random Mixing Transmission : this calculates the total random mixing pressure from all infected persons
     total_random_pressure <- infection_const * sum(beta[i_idx]) * beta[s_idx]  
     p_avoid <- p_avoid * exp(-total_random_pressure) # Probability of avoiding infection from random mixing (using a Poisson approximation).
     
     p_infect <- 1 - p_avoid  # The total probability of getting infected is 1 minus the total probability of avoidance.
     
     newly_exposed_idx <- s_idx[runif(length(s_idx)) < p_infect] # Find and infect newly exposed people using a Bernoulli Trial (to see if they get infected) 
     status[newly_exposed_idx] <- 2 # Move newly infected people to Exposed state.
   }
   
   # Record the day's results
   counts <- tabulate(status, nbins = 4)
   S[t] <- counts[1]; E[t] <- counts[2]; I[t] <- counts[3]; R[t] <- counts[4]
 }
return(list(S = S, E = E, I = I, R = R, t = 1:nt))
}


## STEP 4
plot_nseir <- function(result) {
  
  y_max <- max(result$S, result$E, result$I, result$R)

  plot(result$t, result$S, type = "n", 
       xlim = range(result$t), ylim = c(0, y_max),
       xlab = "Day", ylab = "Number of Individuals")

  # grid(nx=10, ny = 5, col = "gray90", lty = 1)
  abline(v = seq(0, length(result$t)-1, 10), lty = 1, col = "gray90")
  abline(h = seq(0, y_max, 100), lty = 1, col = "gray90")
  
  ## https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
  lines(result$t, result$S, col = "lightseagreen", lwd = 3)
  lines(result$t, result$E, col = "goldenrod2", lwd = 3)
  lines(result$t, result$I, col = "orangered3", lwd = 3)
  lines(result$t, result$R, col = "darkgreen", lwd = 3)

  legend("left", 
         legend = c("Susceptible", "Exposed", "Infectious", "Recovered"),
         col = c("lightseagreen", "goldenrod2", "orangered3", "darkgreen"),
         lwd = 3, bty = "n", cex=0.6)
}

# part 5
# Rprof()

n <- 10000
nc <- 15
h <- get_h(n)
beta <- runif(n)
alink <- get.net(beta, h, nc)

s1 <- nseir(beta, h, alink)
s2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))
const_beta <- rep(mean(beta), length(beta))
alink <- get.net(const_beta, h, nc)
s3 <- nseir(const_beta, h, alink)
s4 <- nseir(const_beta, h, alink, alpha = c(0, 0, 0.04))

# Rprof(NULL)
# print(summaryRprof())

# plotting scenarios:
par(mfrow = c(2, 2))

plot_nseir(s1)
title(main = "Full Model")

plot_nseir(s2)
title(main = "Random Mixing Only")

plot_nseir(s3)
title(main = "Constant Beta")

plot_nseir(s4)
title(main = "Random Mixing & Constant Beta")
