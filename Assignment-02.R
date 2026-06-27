#
# Missing Data Overview
#
data <- read.table("data2.csv", header = TRUE, sep = ",")
str(data)
n <- nrow(data)
#
# Task A 
#Calculate and present appropriate descriptive statistics for all the variables in the dataset. 
#Also, count the missing values for all the variables.
#
### Sum of missing values for each variable
colSums(is.na(data)) #True for missing and False for observed values 
#BMI 4740(47.4%), alt 5287 (52.9%), ggt20 4418(44.2%), crp 5370 (53.7%)


#turning BMI,gender,ggt20 into factors
data$gender <- factor(data$gender, levels = c("male", "female"))
data$bmi <- factor(data$bmi, levels = c("healthy weight", "underweight", "overweight", "obesity"))
data$ggt20 <- factor(data$ggt20, levels = c("0", "1")) #1 higher than >20U/L, 0 else

### Descriptive statistics
#
# Continuous Variables
#
#mean(SD)  na.rm() is needed because of NAs in alt
round(colMeans(data[, c("age", "alb", "alt")], na.rm = TRUE), 2) #Mean
round(sapply(data[ , c("age", "alb", "alt")], sd, na.rm = TRUE) ,2) #SD

mean_alt_cc <- mean(data$alt, na.rm = TRUE) #complete case mean value of alt
mean_alt_cc
#median(IQR)
round(median(data$crp, na.rm = TRUE), 2) #Median
round(IQR(data$crp, na.rm = TRUE), 2) #IQR

#
#Categorical Variables
#
# Frequency tables
table(data$gender)
table(data$bmi)
table(data$ggt20)
# Proportions / Percentages
prop.table(table(data$gender)) * 100
prop.table(table(data$bmi)) * 100
prop.table(table(data$ggt20)) * 100

#### TASK B ####
#Consider the estimation of the mean Alanine levels in our sample. 

#Fit the simple linear regression model of Alanine to Age, and manually fill in the missing values on Alanine 
#based on the predicted values of this model plus a random error term.

datac <- data #copy of dataset

#Missingness Indicator
datac$miss <- 1*is.na(datac$alt) # missingness indicator for Alanine (1=Missing, 0=Not Missing)

#Missingness Mechanism
fitMiss <- glm(miss ~ age,data = datac,family = binomial()) #logistic regression model to see if Age predicts missing indicator
summary(fitMiss) #p-value<0.05 and positive coef => MAR & as age increases , the probability of `Alanine` being missing increases.


# Standard Regression Imputation
lm <- lm(alt ~ age, data = datac)
datac$pred_values <- predict(lm, newdata = datac)

#Manual multiple imputation (m=20)
sd_errors <- sigma(lm)
M <- 20
imp <- list()
alt_mean_estimates <- numeric(20)
#seed
set.seed(1) #rnorm will generate different errors in each iteration

for(i in 1:M) {
  
  # Temporary storing of predicted values
  current_alt_pred <- datac$pred_values 
  
  # Replacing values with observed values (miss == 0)
  current_alt_pred[datac$miss == 0] <- datac$alt[datac$miss == 0] 
  
  # Adding stochastic error to missing values (miss == 1)
  current_alt_pred[datac$miss == 1] <- current_alt_pred[datac$miss == 1] + rnorm(sum(datac$miss == 1), 0, sd_errors)
  
  # Create a copy of the dataset for this iteration
  imputed_data <- datac
  imputed_data$alt <- current_alt_pred  # Replacing the 'alt' column with the imputed column
  alt_mean_estimates[i] <- mean(current_alt_pred)
  #Saving the imputed dataset in the list we created
  imp[[i]] <- imputed_data
}

# Using Rubin’s rules, manually combining the estimates from the 20 vectors.
complete <- function(imp_object, action) {
  return(imp_object[[action]])
}
means_imp <- sapply(1:M, function(i)  mean(complete(imp,i)$alt))
mean_mi <- mean(means_imp) # rubin's estimator for mean of alt

# Calculating the 2.5th and 97.5th percentiles of the 20 estimators as a measure of variability
variability_percentiles <- quantile(means_imp, probs = c(0.025, 0.975))
print(variability_percentiles)

#Rubin’s variance estimator as a second measure of variability
term1 <- mean(sapply(1:M, function(i)  var(complete(imp,i)$alt)/nrow(datac)))

term2 <- (M+1)*sum((means_imp - mean_mi)^2)/(M*(M-1))

Τ_lm <- term1+term2 #Rubin’s variance estimator

#
##Repeating the above using a flexible regression model with a cubic spline (at age 50)
#
#flexible regression model with a cubic spline
fitImpSpl <- lm(alt ~ ns(age,knots = 50),data = datac[datac$miss==0,])
summary(fitImpSpl)

datac$pred_values_spline <- predict(fitImpSpl, newdata = datac)
sd_errors_spline <- sigma(fitImpSpl)

M <- 20
imp_spline <- list()
alt_mean_estimates_spline <- numeric(20)
set.seed(1)

for(i in 1:M) {
  
  current_alt_pred <- datac$pred_values_spline 
  
  current_alt_pred[datac$miss == 0] <- datac$alt[datac$miss == 0] 
  
  current_alt_pred[datac$miss == 1] <- current_alt_pred[datac$miss == 1] + 
    rnorm(sum(datac$miss == 1), 0, sd_errors_spline)
  
  imputed_data <- datac
  imputed_data$alt <- current_alt_pred 
  alt_mean_estimates_spline[i] <- mean(current_alt_pred)
  
  imp_spline[[i]] <- imputed_data
}

#Rubins mean estimator of alt

means_imp_spline <- sapply(1:M, function(i)  mean(complete(imp_spline,i)$alt))
mean_mi_spline <- mean(means_imp_spline) 

# Calculating the 2.5th and 97.5th percentiles
variability_percentiles_spline <- quantile(means_imp_spline, probs = c(0.025, 0.975))
print(variability_percentiles_spline)

#Rubin’s variance estimator
term1_spline <- mean(sapply(1:M, function(i)  var(complete(imp_spline,i)$alt)/nrow(datac)))

term2_spline <- (M+1)*sum((means_imp_spline - mean_mi_spline)^2)/(M*(M-1))

Τ_spline <- term1_spline+term2_spline #Rubin’s variance estimator for spline model


#
##Inverse Probability Weighting (IPW) estimator for the mean Alanine levels
#
#Alt and age have a positive linear relationship
with(datac, summary(lm(alt ~ age))) 
#and the missingness depends on age, so we can get a corrected estimation of mean through IPW

# Probability of missingness
datac$probMiss <- predict(fitMiss,type = "response", newdata = datac)
# IPW weights
datac$ipw <- 1/(1-datac$probMiss)
datac$Cobs <- 1 - datac$miss # essentially needed for the denominator in the IPW estimator
# IPW estimate of the mean
mean_IPW <- round(sum(datac$alt*datac$ipw*datac$Cobs,na.rm = TRUE)
                  /sum(datac$ipw*datac$Cobs),4)


#
##Repeating IPW but improving the missingness model by using a cubic spline for Age (at 46 and 54)
#

fitMiss_spline <- glm(miss ~ bs(age,knots = c(46,54)),data = datac,family = binomial())

# Probability of missingness
datac$probMiss_spline <- predict(fitMiss_spline,type = "response", newdata = datac)
# IPW weights
datac$ipw_spline <- 1/(1-datac$probMiss_spline)
datac$Cobs <- 1 - datac$miss # essentially needed for the denominator in the IPW estimator
# IPW estimate of the mean
mean_IPW_spline <- round(sum(datac$alt*datac$ipw_spline*datac$Cobs,na.rm = TRUE)
                         /sum(datac$ipw_spline*datac$Cobs),4)


##Examining weights distribution
# IPW Linear Model 
summary(datac$ipw)
hist(datac$ipw, main="Distribution of IPW (Linear Model)", xlab="Weights", col="lightblue") #This tells us that while 99% of your observations have small, reasonable weights, a tiny handful of observations have massive, runaway weights in the thousands.

# IPW Spline Model
summary(datac$ipw_spline)
hist(datac$ipw_spline, main="Distribution of IPW (Spline Model)", xlab="Weights", col="lightgreen") #More reasonable tail ~80 and not 15000, but still not the best

#A maximum weight of ~80 means the rarest observed individuals are counting for 80 people, which is statistically manageable.
#A weight of 15000 means a single person was counting for 15000 people(unrealistic)
#which completely destroys the stability of your estimate.


####COMPARISON####
#~General Overview
# Four out of the five methods (Complete-case, Stochastic, Natural Spline, and IPW Spline) 
#yield highly consistent mean values clustered tightly between 14.77 and 15.32. 
#The exception is the Standard IPW (9.85630),
#which severely underestimates the mean compared to every other method.

#~Imputation Models vs. Complete-Case

##Stochastic Imputation (15.31932): 
###Pulls the mean slightly higher than the complete-case baseline.
###Because stochastic imputation adds random residual noise to preserve natural variance
###, it slightly shifts the central tendency here.

##Natural Spline Imputation (14.97612): 
###By allowing a non-linear relationship (1 knot) in the imputation framework,
###this model yields a mean that is almost identical to the complete-case analysis (a tiny difference of just +0.04).

#~IPW vs IPW splines
##Standard IPW (9.85630)is an extreme outlier. 
### massive outlier weights up to 15,000, and dragged the mean down to an artificial 9.85630.

##IPW with Splines (14.77320) completely fixes this issue.
###By introducing 2 knots to allow for non-linear relationships in the weight-modeling stage, 
###the mean shoots right back up into the expected 14.7-15.3 range.

table_means <- data.frame(
  Model = c("complete-case", "stochastic", "natural spline", "IPW", "IPW spline"),
  Mean_Value = c(mean_alt_cc, mean_mi, mean_mi_spline, mean_IPW, mean_IPW_spline)
)
table_means


#### TASK C ####
#Consider as the scientific model the regression of Albumin to Alt, GGT and CRP. 
#Implementation of mice and jomo

#
#Mice with m=10, maxit=15 and method=pmm for all variables
#
datam <- data #Copy for mice
imp2 <- mice(datam, m = 10, maxit = 15, seed = 10, method = "pmm") # changing methods
#imp2 <- readRDS("MICE_pmm.rds") 
fit_mice <- with(imp2, lm(alb ~ alt+ggt20+crp)) #crp not statistically important. 
summary(pool(fit_mice))

#
#Mice with but changing imputation methods, choosing an appropriate imputation method (of choice, instead of pmm)
#
#We can also use: meth <-make.method(datam) and let R define the default method for imputation
imp3 <- mice(datam, m = 10, maxit = 15, seed = 10, method = c("", "", "polyreg", "", "norm", "logreg", "norm")) # changing methods
#imp3 <- readRDS("MICE_parametric.rds")
imp3$method #checking that it used the correct methods
fit_mice_2 <- with(imp3, lm(alb ~ alt+ggt20+crp))
summary(pool(fit_mice_2))

#
#Mice iwth an interaction termAge*Alt for imputation of GGT and CRP
#
ggplot(data, aes(x = age, y = alt)) +
  geom_point(color = "darkred", size = 3) + # Adds the data points
  labs(title = "Patient Age vs. ALT Levels",
       x = "Age (Years)",
       y = "Alanine transaminase / ALT (U/L)") +
  theme_minimal()                           

## Interactions in the imputation models
data_ext <- cbind(datam,  alt.age = NA) # adding a variable that will represent the interaction of the two variables
meth <- c("", "", "polyreg", "", "norm", "logreg", "norm", "~I(alt*age)")
meth

predM <- make.predictorMatrix(data_ext)
predM[c("alt", "age"), "alt.age"] <- 0 # setting to zero, interaction term will not be used for imputing the main effects of it
predM

imp4 <- mice(data_ext,m = 10, maxit = 15, method = meth, predictorMatrix = predM, seed = 10)
#imp4 <- readRDS("MICE_interaction")
fit_mice_3 <- with(imp4, lm(alb ~ alt+ggt20+crp))
summary(pool(fit_mice_3))

#Diagnostics plot to investigate the imputation procedure.
plot(imp4)
densityplot(imp4)

#
#JOMO imputation nimp=10 nburn=5000 nbetween=2000 using all variables
#
Y <- datam[, c("bmi", "alt","ggt20", "crp")]
datam$cons <- 1 # creating a column of 1s
datam$gender<- as.numeric(datam$gender) -1 # declaring gender as numeric for jomo to work. Argument "-1" is used to make reference group = 0
X <- datam[, c( "cons","age", "gender", "alb")] 

set.seed(2)
imp5 <- jomo(Y = Y, X = X, nimp = 10, nburn = 5000, nbetween = 2000)
#imp5 <- readRDS("jomo_imputed_results.rds")
# apply with < model < pool workflow to fit the same model in the imputed dataset
imp.list <- split(imp5, imp5$Imputation)[-1] #split the dataset by the imputation variable and remove the first component of the list created, which corresponds to the initial dataset
imp.list <- imputationList(imp.list) # creation of the imputationList object
lm_jomo <- with(data = imp.list, lm(alb~alt+ggt20+crp))
print(summary(pool(lm_jomo), conf.int = TRUE), digits = 3)

#Complete Case analysis
lm_cc <- lm(alb ~ alt + ggt20 + crp , data = data)
summary(lm_cc)

#The complete case analysis must be firmly rejected due to the massive loss of data (87%) and the risk of selection bias. 
#Method 3 (MICE with Interaction) is favored as the most appropriate approach.
#It incorporates the complex relationships among the variables, avoiding the underestimation of coefficients and ensuring high statistical power.


# #Saving RDS
# #MICE_pmm.rds
# saveRDS(imp2, file = "MICE_pmm.rds")
# #MICE_parametric.rds
# saveRDS(imp3, file = "MICE_parametric.rds")
# #MICE_interaction
# saveRDS(imp4, file = "MICE_interaction")
# #jomo_imputed_results.rds
# saveRDS(imp5, file = "jomo_imputed_results.rds")
