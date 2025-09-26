%% #model build
#aim for this is to use vol as response variable
#use NA, HT, DA, ACh as predictor variables

#this makes the analysis reproducible
set.seed(1)

#scale = TRUE tells R each variable needs to be scaled to have mean 0 and sd 1
#makes sure no predictor variable has more influence in the model if measured in diff units
#validation = CV says to use k-fod cross-validation to evaluate performance of model
# uses k=10 by default 
#could change this here to LOOCV for leave one out cross valdiation
model <- plsr(els_low ~ ., data=exp_low, scale = TRUE, validation = "CV")


#this tells us the root mean squared error (RMSE) - used to evaluate performance of model
#compares observed values against predicted values
#lower = better, more close
summary(model)

#visualise RMSE with the MSE and R2 based on number of PLS components - helps determine best number of components
validationplot(model)
validationplot(model, val.type = "MSEP")
validationplot(model, val.type = "R2")

%%%%%%%%%%%%%%%% #test and train
#now we can use this model to make predictions on new observations
#define training and testing sets
  #need to edit this section once I have my dataset fully loaded and ready to go
train <- exp_low[1:NUMBERchange, c("els_low", "NA", "disp", "drat", "wt", "qsec")]
y_test <- mtcars[numberchange:nrow(mtcars), c("els_low")]
test <- mtcars[numberchange:nrow(mtcars), c("mpg","disp", "drat", "wt", "qsec")]

#use model to make predictions on test set
model <-plsr(hp~mpg+disp+drat+wt+qsec, data=train, scale=TRUE, validation = "CV")
pcr_pred <- predict(model, test, ncomp=2)

#calculate RMSE
#this needs to be interpreted in terms of the scale of your data
#if RMSE is much smaller than typical values, model makes good predictions
#if RMSE is much larger than typical values, model makes bad predictions
sqrt(mean((pcr_pred - y_test)^2))

cor.test(x=els_low, y = ADRA2A_map, method = 'pearson')
cor.test(x=els_low, y = HT1D_map, method = 'pearson')
cor.test(x=els_low, y = HT2C_map, method = 'pearson')
cor.test(x=els_low, y = HT3, method = 'pearson')
cor.test(x=els_low, y = CHRNA2_map, method = 'pearson')
cor.test(x=els_low, y = CHRNA3_map, method = 'pearson')

ggplot(exp_low, aes(x = els_low, y = HT2C_map)) +
  geom_point(color = "blue") +  # Scatterplot points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add regression line (no confidence interval shading)
  ggtitle("Scatterplot of els_low and HT2c_map") +  # Plot title
  xlab("els_low") +  # x-axis label
  ylab("ADRA2A_map") +  # y-axis label
  theme_minimal()  
%%
