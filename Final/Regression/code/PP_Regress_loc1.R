#!/usr/bin/env Rscript
require(dplyr)

MyDF <- read.csv("../data/EcolArchives-E089-51-D1.csv")

MyDF <- na.omit(MyDF[ , c("Predator.mass", "Prey.mass",
                          "Predator.lifestage",
                          "Type.of.feeding.interaction",
                          "Location")])

locations <- unique(MyDF$Location) #only separate by location
#create result data.frame (class)
regress_results <- data.frame(
  Location = character(),
  Type.of.feeding.interaction = character(),
  Predator.lifestage = character(),
  Location = character(),
  Slope = numeric(),
  Intercept = numeric(),
  R.squared = numeric(),
  F.statistic = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

#the loop
for (loc in locations) {
  
  subset_data <- subset(MyDF, Location == loc)
  
  if (nrow(subset_data) > 2) {  # Only regress if more than 2 points, 1point is not enough
    
    # Linear regression, no log transformation
    lm_result <- lm(Predator.mass ~ Prey.mass, data = subset_data)
    lm_summary <- summary(lm_result)
    
    slope     <- lm_summary$coefficients[2, 1]
    intercept <- lm_summary$coefficients[1, 1]
    r_squared <- lm_summary$r.squared
    f_stat    <- lm_summary$fstatistic[1]
    p_val     <- pf(f_stat,
                    lm_summary$fstatistic[2],
                    lm_summary$fstatistic[3],
                    lower.tail = FALSE)
    
    # Append results
    regress_results <- rbind(
      regress_results,
      data.frame(
        Location = loc,
        Type.of.feeding.interaction = as.character(subset_data$Type.of.feeding.interaction[1]),
        Predator.lifestage = as.character(subset_data$Predator.lifestage[1]),
        Location = as.character(subset_data$Location[1]),
        Slope = slope,
        Intercept = intercept,
        R.squared = r_squared,
        F.statistic = f_stat,
        p.value = p_val
      )
    )
  }
}

#save as csv
write.csv(regress_results,
          "../results/PP_Regress_Location_Results1.csv",
          row.names = FALSE)