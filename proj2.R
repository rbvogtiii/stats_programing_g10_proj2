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
  return(alink)
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

## function to update state of each person in network daily
nseir_2 <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005) {
    beta_avg <- mean(beta)
    mixing_const <- (alpha[3] * nc) / (beta_avg^2 * (n - 1))
    
    # states: 1=S, 2=E, 3=I, 4=R
    state <- rep(1, n)
    
    # randomly select initial infected individuals
    start_I <- sample(1:n, round(n*pinf))
    state[start_I] <- 3
    
    # create vectors for each state
    S <- c(nt + 1)
    E <- c(nt + 1)
    I <- c(nt + 1)
    R <- c(nt + 1)
    t <- 0:nt
    
    # initial states
    S[1] <- sum(state == 1)
    E[1] <- 0
    I[1] <- sum(state == 3)
    R[1] <- 0
    
    for (day in 1:nt) {
      new_states <- states
      
      # transitioning from I to R
      current_I <- which(states == 3)
      if (length(current_I) > 0) {
        new_R <- runif(length(current_I)) < delta
        new_states[current_I[new_R]] <- 4
      }
      
      # transitioning from E to I
      current_E <- which(states == 2)
      if (length(current_E) > 0) {
        new_I <- runif(length(current_E)) < gamma
        new_states[current_E[new_I]] <- 3
      }
      
      # transitioning from S to E (getting exposed/infected)
      current_I <- which(states == 3)
      current_S <- which(states == 1)
      
      if (length(current_I) > 0 && length(current_S) > 0) {
        
        new_E <- rep(FALSE, n)
        
        for (i in current_I) {
          # household
          household_S <- which(h == h[i] & state == 1)
          if (length(household_S) > 0) {
            infections <- runif(length(household_S)) < alpha[1]
            new_E[household_S[infections]] <- TRUE
          }
          
          # contacts
          contacts <- alink[[i]]
          if (length(contacts) > 0) {
            contacts_S <- contacts[state[contacts] == 1]
            if (length(contacts_S) > 0) {
              infections <- runif(length(contacts_S)) < alpha[2]
              new_E[contacts_S[infections]] <- TRUE
            }
          }
          
          # random mixing
          if (length(current_S) > 0) {
            prob_random_mixing <- mixing_const * beta[i] * beta[current_S]
            infections <- runif(length(current_S)) < prob_random_mixing
            new_E[current_S[infections]] <- TRUE
          }
        }
        new_states[new_E] <- 2
      }
      states <- new_states
      
      S[day + 1] <- sum(state == 1)
      E[day + 1] <- sum(state == 2)
      I[day + 1] <- sum(state == 3)
      R[day + 1] <- sum(state == 4)
    }
    list(S = S, E = E, I = I, R = R, t = t)
  }
## "Note that while looping through individuals in the I state is probably inevitable, 
## the code can still be made to run in a few seconds for n = 10000, 
## especially with careful use of expressions like x[ind1][ind2] <- y in places
## (assignment to a subvector of a subvector)." - go back and edit to use this!!


nseir_3<- function(beta, h, alink, alpha = c(.1, .01, .01), delta = 0.2, gamma = 0.4, nc = 15, nt = 100, pinf = .005) {
  n <- length(beta)
  alpha_h <- alpha[1]; alpha_c <- alpha[2]; alpha_r <- alpha[3]
  
  infection_const <- (alpha_r * nc) / (mean(beta)^2 * (n - 1))
  status <- rep(1, n); status[sample(n, max(1, round(n * pinf)))] <- 3
  S <- E <- I <- R <- numeric(nt)
  counts <- tabulate(status, nbins = 4)
  S[1] <- counts[1]; E[1] <- counts[2]; I[1] <- counts[3]; R[1] <- counts[4]
  random_const <- (alpha_r * nc) / (mean(beta)^2 * (n - 1))
  
  # Main Simulation Loop
  for (t in 2:nt) {
    s_idx <- which(status == 1)
    e_idx <- which(status == 2)
    i_idx <- which(status == 3)
    
    if (length(e_idx) == 0 && length(i_idx) == 0) {
      S[t:nt] <- S[t-1]; E[t:nt] <- 0; I[t:nt] <- 0; R[t:nt] <- R[t-1]
      break 
    }
    # State Transitions
    # I -> R and E -> I 
    if (length(i_idx) > 0) status[i_idx[runif(length(i_idx)) < delta]] <- 4
    if (length(e_idx) > 0) status[e_idx[runif(length(e_idx)) < gamma]] <- 3
    
    # S -> E Transition
    if (length(s_idx) > 0 && length(i_idx) > 0) {
      
      # Initialise probability of avoiding infection for each susceptible person
      p_avoid <- rep(1, length(s_idx))
      
      # a) Household Infection 
      i_by_hh <- table(h[i_idx]) # Count infected in each household
      s_hh_id <- h[s_idx] # Get household ID for each susceptible
      k <- i_by_hh[as.character(s_hh_id)] # Lookup count of infected for each susceptible 
      k[is.na(k)] <- 0 # If household has no infected, count is 0
      p_avoid <- p_avoid * ((1 - alpha_h)^k) # Probability of avoiding infection from k infectious housemates
      
      # b) Regular Contact Infection
      is_inf <- status == 3 # Create a logical vector for quick lookups
      m <- unlist(lapply(alink[s_idx], function(contacts) sum(is_inf[contacts])))
      p_avoid <- p_avoid * ((1 - alpha_c)^m)
      
      # c) Random Mixing Infection 
      # This calculates the total random mixing pressure from all infected at once
      total_random_pressure <- infection_const * sum(beta[i_idx]) * beta[s_idx]
      p_avoid <- p_avoid * exp(-total_random_pressure)
      
      # Final Bernoulli trial to find newly exposed individuals
      p_infect <- 1 - p_avoid
      newly_exposed_idx <- s_idx[runif(length(s_idx)) < p_infect]
      status[newly_exposed_idx] <- 2
    }
    
    # Record Results
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
