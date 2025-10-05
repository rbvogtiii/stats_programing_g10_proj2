n <- 6
hmax <- 2

## part 1
h <- rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n]

## part 2
net_helper <- function(idx, probs) { # creates a list of connections for the given person
  connections <- which(probs[idx, ] == 1 | probs[, idx] == 1) # find connections
  if (length(connections) == 0) { # no connections
    NA
  } else {
    connections
  }
}

get.net <- function(beta, nc = 15) {
  link_const <- nc / ((mean(beta)**2) * (n - 1))
  probs <- matrix(beta %*% t(beta), n, n) * link_const # initialize probabilities for all pairings
  probs[h[col(probs)] == h[row(probs)] | col(probs) >= row(probs)] <- 0 # get rid of double counting
  rands <- matrix(runif(n * n), n, n)
  probs[(probs - rands) <= 0] <- 0
  probs[probs != 0] <- 1
  lapply(seq_along(beta), net_helper, probs)
}

## part 3
nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
  pop <- 1:n
  I <- sample(pop, n * pinf) # randomly choose pinf% of population to start infected
  S <- pop[!(pop %in% I)] # put rest of population in S
  E <- c()
  R <- c()
  infection_const <- (alpha[3] * nc) / ((mean(beta)**2) * (n - 1))

  for (day in 1:nt) {
    # move from E to I with prob gamma
    I_prob <- gamma - runif(length(E))
    I <- c(I, E[I_prob >= 0])
    E <- E[I_prob < 0]

    # move from I to R with prob delta
    R_prob <- delta - runif(length(I))
    R <- c(R, I[R_prob >= 0])
    I <- I[I_prob < 0]

    # for (i in I) {
      # each member of household exposed with prob alpha[1]

      # each member of network exposed with prob alpha[2]

      # random ppl exposed w prob beta[i]*beta[j]*infection_const
    # }
  }
}

beta <- runif(n)
alink <- get.net(beta, 1)

nseir(beta, h, alink, pinf = 0.50)
