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
  # probs
}

beta <- c(1, 2, 3, 3, 1, 1)
print(h)
tmp <- get.net(beta, 1)
