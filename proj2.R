# distributing n people into households
n = 1000 # population size
hmax = 5 # maximum number of people per household

h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])

# connection network
get.net <- function(beta,h,nc=15) {
  beta_avg <- mean(beta)
  alink <- vector("list",n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (h[i] != h[j]) {
        prob_ij = (nc*beta[i]*beta[j])/(beta_avg^2 *(n - 1))

        if (runif(1) < prob_ij) {
          alink[[i]] <- c(alink[[i]], j)
          alink[[j]] <- c(alink[[j]], i)
        }
      }
    }
  }
  return(alink)
}

## function to update state of each person in network daily
nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005) {
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