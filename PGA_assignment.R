#### Population genetic analysis A1

## Question 1

#Build a wright fisher simulator
simulate_wf <- function(N, seed = 1, s = 0) {
  #Function to make reproducible Wright-Fisher simulator
  set.seed(seed)
  
  #Set up some storage
  c_nt <- data.frame(n = NA, t = 1:500)
  m_max <- 1
  n <- 1
  t <- 1
  
  #Run simulator until fixation or loss
  while (n %in% seq(1, N-1)) {
    
    #Store previous n
    if (t < 500) {c_nt$n[t] <- n}
    
    #Set probability (considering selection)
    p <- n/N
    binom_prob <- p*(1/(1*p + (1 - s)*(1 - p)))
    
    #Do binomial sim
    n_prime <- rbinom(N, 1, binom_prob)
    n <- sum(n_prime)
    
    #Update m_max
    if (n > m_max) {m_max <- n}
    
    #Update t
    t <- t + 1
    
  }
  
  #Ensure rest of ns in c_nt are the same as last value
  if (t < 500) {c_nt$n[t:500] <- n}
  
  return(c(n, t, m_max, c_nt$n))
}


haploid_wf <- sapply(1:1000, simulate_wf, N = 100, s = 0)

#Plot distribution of T where variant fixes and is lost - 
#Consider scaling by proportions
plot_fixed_not <- function(df) {
  fixed <- which(df[1,] == 100)
  
  if (length(fixed) > 0) {
    hist(df[2,fixed], main = 'Distribution of t where variant fixed', xlab = 'T', breaks = 100)
    hist(df[2,-fixed], main = 'Distribution of t where variant lost', xlab = 'T', breaks = 100)
  } else {
    hist(df[2,], main = 'Distribution of t where variant lost', xlab = 'T', breaks = 100)
  }
  
}
par(mfrow = c(1,2))
plot_fixed_not(haploid_wf)


#Plot distribution of m_max
hist(haploid_wf[3,])

#Heatmap plot
c_nt_map <- matrix(0, nrow = 100, ncol = 500)
haploid_wf_n <- haploid_wf[4:dim(haploid_wf)[1],]


make_heatmatrix <- function(df) {
  #Function to convert each column to matrix
  
  #Set up matrix
  c_nt_map <- matrix(0, nrow = 101, ncol = 500)
  
  #Loop over it to set squares
  for (index in 1:1000) {
    for (i in 1:500) {
      c_nt_map[df[i, index] + 1, i] <- c_nt_map[df[i, index] + 1, i] + 1
  
      }
  }
  
  c_nt_map
}

heat_matrix <- make_heatmatrix(haploid_wf_n)
heatmap(sqrt(heat_matrix), Colv = NA, Rowv = NA, scale = 'row', col = rainbow(256))


?heatmap
fig_tb <- as_tibble(heat_matrix)
colnames(fig_tb) <- c(1:500)
fig_tb  <- cbind(tibble(n = 0:100), fig_tb)
fig_tb  <- fig_tb[2:101,] %>% gather('t', 'count', -n)

ggplot(fig_tb, aes(y = n, x = t, fill=count)) + 
  geom_tile() + theme_minimal() + scale_fill_gradientn(colours = rainbow(5))
haploid_wf_df <- data.frame(n = as.vector(unlist(haploid_wf_n)), 
                           t = rep(1:500, times = 1000),
                           run = rep(1:1000, each = 500), 
                           fixed = rep(ifelse(haploid_wf[1,] > 0, T, F), each = 500))







plot_paths <- function(input, title = '') {
  #Function to plot paths of the drift
  
  #Take only rows with paths in (remove top few)
  input_n <- input[4:dim(input)[1],]
  
  #Set rownames
  t <- rownames(input_n) <- 1:500
  
  #Make plot
  plot(x=t, y=input_n[,1], type="l", ylim=c(1,100), col = ifelse(input[1,] > 0, 'green', 'red')[1], 
       xlab = 'T', ylab = 'n', main =  title)
  
  #Add the other lines
  for (i in 2:length(input_n[1,])) {
    lines(x=t ,y=input_n[,i], col = ifelse(input[1,] > 0, 'green', 'red')[i])
  }
  
}

par(mfrow = c(2,2))
plot_paths(diploid_wf_0.02_0)
plot_paths(diploid_wf_0.02_1)
plot_paths(diploid_wf_0.005_0)
plot_paths(diploid_wf_0.005_1)
plot_paths(haploid_wf, 'Haploid no selection')


## Make the colsum things
plot_sums <- function(input, title = '') {
  #Function to plot the distribution of the sums
  
  #Retreive bits we are interested in
  input_n <- input[4:dim(input)[1],]
  
  #Make matrix
  mat <- make_heatmatrix(input_n)
  
  #Get sums - could colour by which ones fixed? - or do a cumulative density plot?
  #plot(apply(mat[2:101,], 2, sum), type = 'l', log = 'y', xlab = 't', ylab = 'count', col = 'green')
  #lines(y = mat[1,], x = 1:length(mat[1,]), col = 'red')
  
  plot(x = 0:100, y = apply(mat, 1, sum), log = 'y')
  
  
}


#### Now repeat with selection
haploid_wf_neg0.02 <- sapply(1:1000, simulate_wf, N = 100, s = -0.02)
haploid_wf_neg0.005 <- sapply(1:1000, simulate_wf, N = 100, s = -0.005)
haploid_wf_0.02 <- sapply(1:1000, simulate_wf, N = 100, s = 0.02)
haploid_wf_0.005 <- sapply(1:1000, simulate_wf, N = 100, s = 0.005)

par(mfrow = c(1,2))
plot_fixed_not(haploid_wf_neg0.02)
plot_fixed_not(haploid_wf_neg0.005)
plot_fixed_not(haploid_wf_0.02)
plot_fixed_not(haploid_wf_0.005)

par(mfrow = c(2,2))
plot_paths(haploid_wf_neg0.02)
plot_paths(haploid_wf_neg0.005)
plot_paths(haploid_wf_0.02)
plot_paths(haploid_wf_0.005)



#### Now do diploid simulation
#Build a wright fisher simulator
simulate_diploid_wf <- function(N, seed = 1, s = 0, h = 0) {
  #Function to make reproducible Wright-Fisher simulator
  set.seed(seed)
  
  #Set up some storage
  c_nt <- data.frame(n = NA, t = 1:500)
  m_max <- 1
  n <- 1
  t <- 1
  
  #Run simulator until fixation or loss
  while (n %in% seq(1, N-1)) {
    
    #Store previous n
    if (t < 500) {c_nt$n[t] <- n}
    
    #Set probability (considering selection)
    p <- n/(2*N)
    
    numerator <- p^2 + p*(1 - p)*(1 - h*s)
    denominator <- p^2 + 2*p*(1 - p)*(1 - h*s) + (1 - p)^2 * (1 - s)
    binom_prob <- numerator/ denominator
    
    #Do binomial sim
    n_prime <- rbinom(2*N, 1, binom_prob)
    n <- sum(n_prime)
    
    #Update m_max
    if (n > m_max) {m_max <- n}
    
    #Update t
    t <- t + 1
    
  }
  
  #Ensure rest of ns in c_nt are the same as last value
  if (t < 500) {c_nt$n[t:500] <- n}
  
  return(c(n, t, m_max, c_nt$n))
}

#Make average path
weighted_heatmap <- heat_matrix
for (i in 0:100) {
  
  
}


### Question 3

#Set up some paramters
total_adults = 30923
HbAA = 25374
HbAS = 5482
HbSS = 67

#Calculate allele frequencies
freq_A = (HbAA*2 + HbAS)/(total_adults*2)
freq_S = (HbSS*2 + HbAS)/(total_adults*2)

#Calculate expected values if HWE = true
exp_AA = (freq_A^2)*total_adults
exp_AS = (2*freq_A*freq_S)*total_adults
exp_SS = (freq_S^2)*total_adults

#Make contingency table
chi_table <- as.table(rbind(c(HbAA, HbAS, HbSS), c(exp_AA, exp_AS, exp_SS)))

dimnames(chi_table) <- list(real = c("Real", "Expected"),
                    genotype = c("HbAA","HbAS", "HbSS"))

(Xsq <- chisq.test(chi_table))  # Prints test summary
Xsq$observed   # observed counts (same as M)
Xsq$expected   # expected counts under the null
Xsq$residuals  # Pearson residuals
Xsq$stdres     # standardized residuals

#3b
N = 40000

#Get birth numbers under HWE
born_AA = (freq_A^2)*N
born_AS = (2*freq_A*freq_S)*N
born_SS = (freq_S^2)*N

#Get genotype fitnesses
W_SS = (HbSS/born_SS)
W_AA = (HbAA/born_AA)/W_SS
W_AS = (HbAS/born_AS)/W_SS

W_SS = 1

#Calculate S and H
s = 1 - W_AA
h = (1 - W_AS)/s



#3b
N = 32000

#Get birth numbers under HWE
born_AA = (freq_A^2)*N
born_AS = (2*freq_A*freq_S)*N
born_SS = (freq_S^2)*N

#Get genotype fitnesses
W_SS = (HbSS/born_SS)
W_AA = (HbAA/born_AA)/W_SS
W_AS = (HbAS/born_AS)/W_SS

W_SS = 1

#Calculate S and H
s = 1 - W_AA
h = (1 - W_AS)/s


#### Question 4
calc_fixation <- function(f_ai,  N = 100) {
  #Function to calculate fixation probability
  
  numerator = 1 - exp(-2 * f_ai)
  denominator = 1 - exp(-2 * f_ai * N)
  
  numerator/denominator
}

range(calc_fixation(rexp(1000, 1/0.01)))



#Do gibbs sampling update and then track number of genes in state 1
simulate_population <- function(cycles, N = 100, genes = 1000, seed = 1, rate_change = 0) {
  #Function to do simulation
  
  #Set up a vector with the genes
  genotype <- rep(0, times = 1000)
  n1 <- rep(NA, times = cycles)
  
  #Make it reproducible
  set.seed(seed)
  
  #Set selection coefficients
  f_ais <- rexp(genes, 1/0.01)
  
  #Get fixation probabilities
  pi_p <- calc_fixation(f_ais)
  pi_p_df <- data.frame(lost = 1 - pi_p, fixed = pi_p)
  
  #Get reverse fixation probabilities
  pi_p_neg <- calc_fixation(-f_ais)
  pi_p_neg_df <- data.frame(lost = 1 - pi_p_neg, fixed = pi_p_neg)
  
  #Work through the simulation i times
  for (i in 1:cycles) {

    #Get outcomes
    ones <- which(genotype == 1)
    
    #Make new genotype
    new_genotype <- genotype
    
    #Allow random reversing of direction
    if (rpois(1, rate_change) > 0) {
      pi_p_neg_df_temp <- pi_p_neg_df
      pi_p_neg_df <- pi_p_df
      pi_p_df <- pi_p_neg_df_temp
      
    }
    
    #Sample according to fixation probabilities (different for pos and neg)
    if (length(ones) > 0) {
      
      #If some in both states sample
      new_genotype[-ones] <- apply(pi_p_df[-ones,], 1, sample, x = c(0, 1), size  =  1, replace = T)
      new_genotype[ones] <- apply(pi_p_neg_df[ones,], 1, sample, x = c(1, 0), size  =  1, replace = T)
    } else if (length(ones) == length(new_genotype)) {
      
      #If all ones
      new_genotype <- apply(pi_p_neg_df, 1, sample, x = c(1, 0), size  =  1, replace = T)
    } else {
      
      #If all zeros
      new_genotype <- apply(pi_p_df, 1, sample, x = c(0, 1), size  =  1, replace = T)
    }
    
    
    #Track statistics
    n1[i] <- sum(new_genotype)
    
    #Update the genotype
    genotype <- new_genotype
  }
  
  return(list(n1, genotype, f_ais))
  
}

n1_time <- simulate_population(100000, rate_change = 1e-4)[[1]]
plot(x = 1:length(n1_time), y = n1_time, type = 'l')

#Could plot a log odds ratio type curve (sigmoid curve?) for which ones end up as 1 and which as zero? 
#Also do this as proportions over time so track how often each is a 1 and which is a zero at equilibrium

#Should record when a switch occurs in terms of t,  so  that can plot that on the graphX









