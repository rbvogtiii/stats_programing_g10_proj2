## part 1
get_h <- function(n, hmax = 5) {
  rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n]
}

## part 2
net_helper <- function(idx, probs) { # creates a list of connections for the given person
  connections <- which(probs[idx, ] == 1 | probs[, idx] == 1) # find connections
  if (length(connections) == 0) { # no connections
    NA
  } else {
    connections
  }
}

get.net <- function(beta, h, nc = 15) {
  n <- length(beta)
  link_const <- nc / ((mean(beta)**2) * (n - 1))
  probs <- outer(beta, beta) * link_const # initialize probabilities for all pairings
  for (i in 1:n) probs[i, i:n] <- 0 # get rid of double counting
  probs[outer(h, h, FUN = "==")] <- 0 # remove connections within households
  rands <- matrix(runif(n * n), n, n)
  probs[probs >= rands] <- 1
  lapply(seq_along(beta), net_helper, probs)
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
      next
    }

    # alternative for household exposures
    E_prob <- 1 - (((1 - alpha[1])**tabulate(h[pop %in% I], nbins = max(h)))[h[S]]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]

    if (length(I) == 0) {
      next
    }

    if (length(unlist(alink[I])) == 0) {
      next
    }

    # alternative for network exposures
    E_prob <- 1 - ((1 - alpha[2])**tabulate(unlist(alink[I]), nbins = n)[S]) - runif(length(S))
    E <- c(E, S[E_prob >= 0])
    S <- S[E_prob < 0]

    if (length(I) == 0) {
      next
    }

    # alternative for random exposures
    # update to use outer
    E_prob <- 1 - ((t(matrix(beta[S], nrow = length(S), ncol = length(I))) * beta[I]) * infection_const)
    E_prob <- 1 - apply(E_prob, 2, prod) - runif(length(S))
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

nseir_2 <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, nc = 15, nt = 100, pinf = .005) {
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

    for (i in I) {
      # each member of household exposed with prob alpha[1]
      household <- S[h[S] == h[i]]
      E_prob <- alpha[1] - runif(length(household))
      E <- c(E, household[E_prob >= 0])

      # each member of network exposed with prob alpha[2]
      network <- S[S %in% alink[i]]
      E_prob <- alpha[2] - runif(length(network))
      E <- c(E, network[E_prob >= 0])

      # random ppl exposed w prob beta[i]*beta[j]*infection_const
      E_prob <- (beta[i] * beta[pop %in% S] * infection_const) - runif(length(S))
      E <- c(E, S[E_prob >= 0])
    }

    E <- unique(E)
    S <- S[!(S %in% E)]

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

# part 4
n <- 10000
nc <- 15
h <- get_h(n)
beta <- runif(n)
alink <- get.net(beta, h, nc)

# print(system.time(s1 <- nseir(beta, h, alink, pinf = 0.5)))
# print(system.time(s2 <- nseir_2(beta, h, alink, pinf = 0.5)))
# plot(s1$E, type='l', ylim = c(0, max(s1$S, s2$S)))
# lines(s1$I)
# lines(s1$S)
# lines(s2$E, col = 2)
# lines(s2$I, col = 2)
# lines(s2$S, col = 2)

s1 <- nseir(beta, h, alink)
s2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))
const_beta <- rep(mean(beta), length(beta))
alink <- get.net(const_beta, h, nc)
s3 <- nseir(const_beta, h, alink)
s4 <- nseir(const_beta, h, alink, alpha = c(0, 0, 0.04))

plot(s1$t, s1$S, type = "l", xlab = "day", ylim = c(0, n))
lines(s1$t, s1$I)
# lines(s1$t, s1$E)

lines(s2$t, s2$S, col = 2)
lines(s2$t, s2$I, col = 2)
# lines(s2$t, s2$E, col = 2)

lines(s3$t, s3$S, col = 3)
lines(s3$t, s3$I, col = 3)
# lines(s3$t, s3$E, col = 3)

lines(s4$t, s4$S, col = 4)
lines(s4$t, s4$I, col = 4)
# lines(s4$t, s4$E, col = 4)
