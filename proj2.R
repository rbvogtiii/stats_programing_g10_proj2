# distributing n people into households
n = 1000 # population size
hmax = 5 # maximum number of people per household

h <- sample(rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n])


# connection network
get.net <- function(beta,h,nc=15) {
  beta_avg <- mean(beta)
  links <- vector("list",n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (h[i] != h[j]) {
        prob_ij = (nc*beta[i]*beta[j])/(beta_avg^2 *(n - 1))

        if (runif(1) < prob) {
          links[[i]] <- c(links[[i]], j)
          links[[j]] <- c(links[[j]], i)
        }
      }
    }
    return(links)
  }
}
