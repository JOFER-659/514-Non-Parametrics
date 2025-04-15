##########################################################################
# STAT 514 Project
####################
# BY: Jacob Hofer
##########################################################################

# This project is intended to show how kruskal wallis and One-Way ANOVA 
# power differs when effect size and strength of unbalance in the design
# of the experiment vary 

# The first part of the code will show how power differs when we have a balanced
# case for both methods and same effect size, with sample size increasing

# The second part of the code will show how power differs when total sample size
# is the same and effect size is the same, but variability of sample sizes at
# each level increases
# (need to write code for finding 3rd obs when effect size is constant)

# The third part of the code will explore how power differs when total sample
# size and variance is the same but effect (true difference in group means) 
# increases (need to write code for finding 3rd obs when variance is constant).

# NOTE: variance of normal distributions for each level is always equal. 

# Necessary Libraries
library(tidyverse)
library(reshape2)
set.seed(11)
############################# Part 1 ####################################

# For one sample 

# Set parameters

# Group Means
g1mean <- 10
g2mean <- 4
g3mean <- 5
var <- 4
Nseq <- seq(from = 9, to = 51, by = 3)


B <- 1000 #Number of Sims per Nseq
alpha <- 0.05 # Type 1 error rate
app <- numeric(length = length(Nseq))
kpp <- numeric(length = length(Nseq))
for (j in 1:length(Nseq)){
  N <- Nseq[j]
  n1 = n2 = n3 = (N/3)
kp <- numeric(length = B)
ap <- numeric(length = B)
for (i in 1:B){
p1df <- data.frame(group = c(rep('g1', times = n1), rep('g2', times = n2),
                             rep('g3', times = n3)),
                   val = c(rnorm(n1, g1mean, var), rnorm(n2, g2mean, var),
                           rnorm(n3, g3mean, var)))

p1df$group <- factor(p1df$group)

test <- kruskal.test(x = p1df$val, g = p1df$group)
kp[i] <-test$p.value

ap[i] <- anova(lm(val ~ group, data = p1df))$'Pr(>F)'[1]
}
app[j] <- mean(ap < alpha)
kpp[j] <- mean(kp < alpha)
}

#Put simulation results in data frame
df1 <- data.frame(Sample_Size = Nseq, Kruskal = kpp, Anova = app)

#Plot the data
df_long1 <- melt(df1, id.vars = "Sample_Size", variable.name = "Test", value.name = "Power")

ggplot(df_long1, aes(x = Sample_Size, y = Power, color = Test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(title = "Power vs. Sample Size",
       x = "Total Sample Size",
       y = "Power",
       color = "Test Type") +
  theme_minimal()

############################# Part 2 ####################################

#### Holding n3 (middle most mean)
simulate_imbalance <- function(group_means, var, B, alpha, imbalance_vals, fixed_n, fix = "n3") {
  results <- data.frame(Imbalance = numeric(), Kruskal = numeric(), ANOVA = numeric())
  
  g1mean <- group_means[1]
  g2mean <- group_means[2]
  g3mean <- group_means[3]
  
  for (k in imbalance_vals) {
    if (fix == "n3") {
      n1 <- fixed_n - k
      n2 <- fixed_n + k
      n3 <- fixed_n
    } else if (fix == "n2") {
      n1 <- fixed_n - k
      n2 <- fixed_n
      n3 <- fixed_n + k
    } else if (fix == "n1") {
      n1 <- fixed_n
      n2 <- fixed_n - k
      n3 <- fixed_n + k
    }
    
    # Skip invalid combinations
    if (n1 <= 1 || n2 <= 1 || n3 <= 1) next
    
    kp <- numeric(B)
    ap <- numeric(B)
    
    for (i in 1:B) {
      df <- data.frame(
        group = factor(c(rep("g1", n1), rep("g2", n2), rep("g3", n3))),
        val = c(rnorm(n1, g1mean, var), rnorm(n2, g2mean, var), rnorm(n3, g3mean, var))
      )
      
      kp[i] <- kruskal.test(val ~ group, data = df)$p.value
      ap[i] <- anova(lm(val ~ group, data = df))$'Pr(>F)'[1]
    }
    
    results <- rbind(results, data.frame(
      Imbalance = k,
      Kruskal = mean(kp < alpha),
      ANOVA = mean(ap < alpha)
    ))
  }
  
  return(results)
}


####### Running simulations #########

n3result <- simulate_imbalance(c(10, 4, 5), var = 4, B = 10000, alpha = 0.05, 
                               imbalance_vals = -9:9, fixed_n = 12, fix = "n3")
n2result <- simulate_imbalance(c(10, 4, 5), var = 4, B = 10000, alpha = 0.05, 
                               imbalance_vals = -9:9, fixed_n = 12, fix = "n2")
n1result <- simulate_imbalance(c(10, 4, 5), var = 4, B = 10000, alpha = 0.05, 
                               imbalance_vals = -9:9, fixed_n = 12, fix = "n1")

n3df_long <- melt(n3result, id.vars = "Imbalance", variable.name = "Test", value.name = "Power")
ggplot(n3df_long, aes(x = Imbalance, y = Power, color = Test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Power vs. Sample Size Imbalance (Fixing n3 = 12)",
    x = "Imbalance (n2 - n1)",
    y = "Power",
    color = "Test Type"
  ) +
  theme_minimal()


n2df_long <- melt(n2result, id.vars = "Imbalance", variable.name = "Test", value.name = "Power")
ggplot(n2df_long, aes(x = Imbalance, y = Power, color = Test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Power vs. Sample Size Imbalance (Fixing n2 = 12)",
    x = "Imbalance (n3 - n1)",
    y = "Power",
    color = "Test Type"
  ) +
  theme_minimal()


n1df_long <- melt(n1result, id.vars = "Imbalance", variable.name = "Test", value.name = "Power")
ggplot(n1df_long, aes(x = Imbalance, y = Power, color = Test)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Power vs. Sample Size Imbalance (Fixing n1 = 12)",
    x = "Imbalance (n3 - n2)",
    y = "Power",
    color = "Test Type"
  ) +
  theme_minimal()


################ part 3 #####################################################