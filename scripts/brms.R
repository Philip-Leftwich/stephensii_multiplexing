# Load the brms package
library(brms)

# Define the Bayesian binomial mixed-effects model with interactions
model <- 
  brm(cbind(win, loss)| trials(win+loss)   ~ cas9_parent * gRNA_type * pre_cut + (1 | id),
  family = binomial(),
  data = homing_data
)

# Summarize the model
summary(model)

# Plot the results
plot(model)
