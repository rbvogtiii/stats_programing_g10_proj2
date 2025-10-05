n <- 1000
hmax <- 5

## part 1
h <- rep(1:n, times = sample(1:hmax, n, replace = TRUE))[1:n]