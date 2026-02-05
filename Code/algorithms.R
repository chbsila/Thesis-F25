library(BiGER)

rank_aggregation <- function(r) {
  
  return(BiGER::BiGER(r,
                      n_r = rep(nrow(r), ncol(r)),
                      n_u = rep(0, ncol(r)))
  )
}

sim <- function(n_items, n_processes, rho=0.5){
  # ```
  # This function simulates ranked lists from the latent variable model proposed
  # by Wang et al. (2025). For now, let's focus on fully ranked lists first.
  
  # (int) num_items: Total number of items being ranked.
  # (int) num_processes: Total number of processes.
  # (vec) rho: The correlation between a gene's local importance and global importance.
  
  # Truth
  sigma_s2=rho^(-2)-1 #study variance
  mu_i <- rnorm(n_items)
  true_rank=rank(-mu_i)
  
  ## simulate data for each process
  w  <- matrix(NA,n_items,n_processes)
  r  <- matrix(NA,n_items,n_processes)
  
  for (p in 1:n_processes) {
    w[,p] <- mu_i + rnorm(n_items, mean=0, sd=sqrt(sigma_s2))
    r[,p] <- rank(-w[,p])
  }
  
  return(list(r=r,true_rank=true_rank))
}

sim_ranking <- function(n_good, n_bad, n_items) {
  # ```
  # This function simulates a batch of good ranking and a batch of
  # byzantine ranking. All the byzantine rankings might collude with each other.
  
  # (int) n_good: Total number of good processes
  # (int) n_bad: Total number of bad processes.
  # (int) n_items: Total number of items being ranked.
  
  good <- sim(n_items, n_good)
  bad <- sim(n_items, n_bad) 
  return(list("good" = good, "bad" = bad))
}

intersection_model_A <- function(n_good, n_bad){
  # ```
  # This function constructs intersection model A.
  # Each good process hears from itself, all t Byzantine processes and (n - 2t - 1) other good processes 
  # sampled uniformly at random
  
  # (int) n_good: Total number of good processes.
  # (int) n_bad: Total number of Byzantine processes.
  
  # network[i, j] = 1 means good receiver i receives a message from sender j.
  
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  byzantine_ids <- (n_good+1):n
  network <- matrix(0, nrow = n_good, ncol = n)
  
  for (i in good_ids){
    network[i, i] <- 1
    network[i, sample(setdiff(good_ids, i), n-2*t-1)] <- 1   # n-2t-1 good processes at random
    for (j in byzantine_ids){
      network[i,j] <- 1 # Byzantine inputs are always present 
    }
  }
  return(network)
}

intersection_model_B <- function(n_good, n_bad){
  # ```
  # This function constructs intersection model B.
  # TO DO: EXPLAIN
  
  # (int) n_good: Total number of good processes.
  # (int) n_bad: Total number of Byzantine processes.

  # network[i, j] = 1 means good receiver i receives a message from sender j.
  
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  network <- matrix(0, nrow = n_good, ncol = n)
  # Each process P_i samples n - t total including itself
  sample_size <- (n - t - 1) # other than itself
  for (i in good_ids) {
    network[i, i] <- 1  # M = {i}
    candidates <- setdiff(1:n, i) 
    S_i <- sample(candidates, sample_size) # uniformly at random 
    for (s in S_i) {
      network[i, s] <- 1  # M = M \union {P_s}
      good_candidates <- setdiff(good_ids, i)
      network[i, sample(good_ids, t)] <- 1   # M = M \union {t good processes at random}
    }
  }
  return(network)
}

intersection_model_C <- function(n_good, n_bad) {
  # ```
  # This function constructs intersection model C.
  # We build nested super sets A_1, A_2, ..., A_{t+1}, where:
  # |A_1| = n-t and each A_{k+1} adds one new process from the remainder
  # Then, each good process selects one A_k uniformly at random and hears from
  # all processes in that super set.
  
  # (int) n_good: Total number of good processes.
  # (int) n_bad: Total number of bad (Byzantine) processes.
  
  # network[i, j] = 1 means good receiver i receives a message from sender j.
  
  n <- n_good + n_bad
  t <- n_bad
  good_ids <- 1:n_good
  bad_ids  <- (n_good + 1):n
  all_ids <- 1:n
  network <- matrix(0, nrow = n_good, ncol = n)
  A_sets <- vector("list", t + 1)
  A_sets[[1]] <- sample(all_ids, n-t)
  remainder <- setdiff(all_ids, A_sets[[1]])
  
  remainder <- sample(remainder, t)  # randomize order
  
  for (k in 2:(t + 1)) {
    A_sets[[k]] <- c(A_sets[[k - 1]], remainder[k - 1])
  }
  
  # Assign each good process a random superset and fill the network 
  membership <- sample(1:(t + 1), size = n_good, replace = TRUE)
  for (i in 1:n_good) {
    heard_from <- A_sets[[membership[i]]]
    network[i, heard_from] <- 1
  }
  
  return(network)
}

sim_round <- function(r, network) {
  # ```
  # This function simulates a round of communication given a network configuration.
  # (matrix) r_all: Matrix of rankings (good + bad) used in this round
  # (matrix) network: Pre-computed communication network/intersection
  # (integer) n_good: number of good processes
  # (integer) n_bad: number of bad processes

    n_good <- ncol(r$good$r)
    n_bad  <- ncol(r$bad$r)
    r_all  <- cbind(r$good$r, r$bad$r)
    
    results <- list()
    trust_matrix <- matrix(0, nrow = n_good, ncol = n_good + n_bad)
    
    # Good processes
    for (i in 1:n_good) { 
      ids_used <- which(network[i, ] == 1)
      r_heard <- r_all[, ids_used, drop = FALSE]
      
      results[[i]] <- list()
      results[[i]][["id_good"]] <- i
      results[[i]][["ids_used"]] <- ids_used
      results[[i]][["ra"]] <- rank_aggregation(r_heard)
      results[[i]][["posterior_ranking"]] <- rank(-results[[i]]$ra$mu)
      
      sigmas <- results[[i]]$ra$sigma2
      sigmas <- (sigmas - min(sigmas)) / (max(sigmas) - min(sigmas))
      trust_matrix[i, ids_used] <- 1 - sigmas
    }
    
    # Bad processes
    for (j in 1:n_bad) {
      idx <- n_good + j
      results[[idx]] <- list()
      results[[idx]][["posterior_ranking"]] <- sample(1:nrow(r_all))
    }
    
    return(list(results = results, trust_matrix = trust_matrix))
  }
  
