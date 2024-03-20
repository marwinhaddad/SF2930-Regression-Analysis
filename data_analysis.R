setwd("C:/Users/marwi/Documents/RProjects/Project1")

library(caret)
library(ggplot2)
library(MASS)

set.seed(0)

convert2SI <- function(data){
  data$density <- data$density * 1000
  data$weight <- data$weight * 0.4535924
  data$height <- data$height * 0.0254
  data$neck <- data$neck * 0.01
  data$chest <- data$chest * 0.01
  data$abdomen <- data$abdomen * 0.01
  data$hip <- data$hip * 0.01
  data$thigh <- data$thigh * 0.01
  data$knee <- data$knee * 0.01
  data$ankle <- data$ankle * 0.01
  data$biceps <- data$biceps * 0.01
  data$forearm <- data$forearm * 0.01
  data$wrist <- data$wrist * 0.01
  return(data)
}


# read data
data <- read.csv("bodyfatmen.csv")
data <- convert2SI(data)
data.model <- lm(density ~., data=data)

# data summary
summary(data.model)
anova(data.model)

# external studentized residuals of FULL MODEL
rstud <- rstudent(data.model)

# norm plots of r student
qqnorm(rstud)
qqline(rstud)

# res vs fitted
plot(data.model$fitted.values, rstud,
     main="R-student vs. Fitted Values", 
     xlab="Fitted values", ylab="R-sudent"
     )
abline(h = 0, lty = "dashed", col = "gray")

regressors <- names(data)[names(data) != "density"]
par(mfrow = c(4, 4))  # Adjust the layout according to the number of regressors
for (regressor in regressors) {
  plot(data[[regressor]], data.model$residuals,
       main=paste("R-student vs.", regressor), 
       xlab=regressor, ylab="R-student")
  abline(h = 0, lty = "dashed", col = "gray")
}

plot(data.model, which = 1:5)
# 39 seems very influential

library(car)  # avPlots
avPlots(data.model)

res_df <- data.frame(fit = data.model$fit, residuals = rstud)

ggplot(res_df, aes(x = fit, residuals)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0, color = "gray", linetype = "dashed") +
  labs(title = "Studentized Residuals vs. Fitted Values",
       y = "Studentized Residuals",
       x = "Fitted Values") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.background=element_rect(colour="black"))


ggplot(qq_data, aes(x = Theoretical, y = Residuals)) +
  geom_point() +  # Plot all points
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +  # Add a reference line (y = x)
  labs(title = "QQ Plot of Studentized Residuals",
       x = "Theoretical Quantiles", y = "Studentized Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# ------------- question 2 --------------- #

infl <- influence.measures(data.model)
summary(infl)

infl.index <- c(5, 31, 36, 39, 41, 52, 83, 102, 155, 171, 200, 202, 203, 217, 220)
rstud.infl <- rstud[infl.index]

regressors <- names(data)[names(data) != "density"] # Exclude the response variable

data$Label <- ifelse(seq_len(nrow(data)) %in% infl.index, as.character(seq_len(nrow(data))), NA)

# Loop through each regressor to create a plot
for(regressor in regressors) {
  p <- ggplot(data, aes_string(x = regressor, y = "density")) +
    geom_point() +  # Plot all points
    geom_point(data = data[infl.index, ], aes_string(x = regressor, y = "density"), color = "red") +  # Highlight influential points in red
    geom_text(aes_string(label = "Label"), vjust = -0.5, hjust = -0.5, color = "red", size = 3) +  # Add labels to influential points
    geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Add a linear regression line
    labs(title = paste("Density vs", regressor),
         x = regressor, y = "Density") +
    theme_minimal() +  # Use a minimal theme, which includes a grid
    theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
          panel.grid.minor = element_line(color = "grey90", size = 0.25)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)  # Display the plot
}

# without outlier labels
for(regressor in regressors) {
  p <- ggplot(data, aes_string(x = regressor, y = "density")) +
    geom_point() +  # Plot all points
    geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Add a linear regression line
    labs(title = paste("Density vs", regressor),
         x = regressor, y = "Density") +
    theme_minimal() +  # Use a minimal theme, which includes a grid
    theme(panel.grid.major = element_line(color = "grey80", size = 0.5),
          panel.grid.minor = element_line(color = "grey90", size = 0.25)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)  # Display the plot
}

theoretical_quantiles <- qqnorm(rstud, plot.it = FALSE)$x
qq_data <- data.frame(Theoretical = theoretical_quantiles, Residuals = rstud)

# Subset data for influential points to highlight
infl_qq_data <- qq_data[infl.index, ]

# QQ plot using ggplot2
ggplot(qq_data, aes(x = Theoretical, y = Residuals)) +
  geom_point() +  # Plot all points
  geom_point(data = infl_qq_data, aes(x = Theoretical, y = Residuals), color = "red") +  # Highlight influential points
  geom_text(data = infl_qq_data, aes(x = Theoretical, y = Residuals, label = rownames(infl_qq_data)), vjust = -1, hjust = 1, color = "red", size = 4) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +  # Add a reference line (y = x)
  labs(title = "QQ Plot of Studentized Residuals",
       x = "Theoretical Quantiles", y = "Studentized Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# QQ plot using ggplot2 without outlier labels
ggplot(qq_data, aes(x = Theoretical, y = Residuals)) +
  geom_point() +  # Plot all points
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +  # Add a reference line (y = x)
  labs(title = "QQ Plot of Studentized Residuals",
       x = "Theoretical Quantiles", y = "Studentized Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

rvf_data <- data.frame(Fitted = fitted(data.model), Residuals = rstud)

# Residuals vs. Fitted Values plot
ggplot(rvf_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_point(data = rvf_data[infl.index, , drop = FALSE], aes(x = Fitted, y = Residuals), color = "red") +
  geom_text(data = rvf_data[infl.index, , drop = FALSE], aes(x = Fitted, y = Residuals, label = rownames(rvf_data)[infl.index]), 
            hjust = 1.3, color = "red", size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  labs(title = "Residuals vs. Fitted Values",
       x = "Fitted Values", y = "Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Residuals vs. Fitted Values plot without outlier labels
ggplot(rvf_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  labs(title = "Residuals vs. Fitted Values",
       x = "Fitted Values", y = "Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

sqrt_std_resid <- sqrt(abs(rstudent(data.model)))

# Create a data frame for plotting
scale_loc_data <- data.frame(
  Fitted = fitted(data.model),
  SqrtStdResid = sqrt_std_resid
)

# Add a column to identify the influential points
scale_loc_data$Influential <- ifelse(1:nrow(scale_loc_data) %in% infl.index, "Yes", "No")

# Create the Scale-Location plot using ggplot2
scale_loc_plot <- ggplot(scale_loc_data, aes(x = Fitted, y = SqrtStdResid, color = Influential)) +
  geom_point() +  # Plot all points
  scale_color_manual(values = c("No" = "black", "Yes" = "red")) +  # Color influential points in red
  geom_smooth(method = "loess", se = FALSE, aes(color = NULL)) +  # Add a loess smoothed line
  labs(title = "Scale-Location Plot",
       x = "Fitted Values", y = "Sqrt of Absolute Standardized Residuals") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide the legend

# Print the plot
print(scale_loc_plot)

print(abs(rstud[infl.index]))
outliers <- which(abs(rstud[infl.index]) > 2)

data = data[,!(names(data) %in% c("Label"))]

plot(data.model, which=4)

# -------------- remove outliers --------------- #
removed_outliers <- data[-c(31, 39, 83, 200, 217, 203, 220), ]
removed_outliers.model <- lm(density ~., data=removed_outliers)
summary(removed_outliers.model)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -20.7148  -6.6085   0.1854   6.4895  20.7061 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1149.95461   54.32657  21.167  < 2e-16 ***
#   age           -0.16185    0.07267  -2.227  0.02691 *  
#   weight         0.38616    0.34410   1.122  0.26294    
# height         0.08338   17.26777   0.005  0.99615    
# neck          98.08663   51.90774   1.890  0.06008 .  
# chest          9.19124   24.19683   0.380  0.70441    
# abdomen     -210.58036   20.50522 -10.270  < 2e-16 ***
#   hip           38.41924   32.27954   1.190  0.23520    
# thigh        -68.76671   32.45698  -2.119  0.03520 *  
#   knee          18.50962   54.40360   0.340  0.73400    
# ankle        -25.40374   60.33785  -0.421  0.67413    
# biceps       -35.34640   37.71797  -0.937  0.34969    
# forearm      -76.10186   45.54304  -1.671  0.09610 .  
# wrist        350.47054  121.10116   2.894  0.00417 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 9.294 on 228 degrees of freedom
# Multiple R-squared:  0.7634,	Adjusted R-squared:  0.7499 
# F-statistic: 56.58 on 13 and 228 DF,  p-value: < 2.2e-16

# --------------

data.no.outliers <- removed_outliers
model.no.outliers <- lm(density ~., data=data.no.outliers)
t.no.outliers <- rstudent(model.no.outliers)

summary(model.no.outliers)

plot(model.no.outliers, which=1:5)

# -------------------- plot residuals vs regressors

# Loop through each regressor
par(mfrow = c(4, 4)) 


for(regressor in names(data.no.outliers)[-which(names(data.no.outliers) == "density")]) {
  plot(data.no.outliers[[regressor]], t.no.outliers, 
       main = paste("Residuals vs", regressor),
       xlab = regressor, 
       ylab = "Residuals")
  abline(h = 0, col = "gray") # Add a horizontal line at 0
}

data.no.outliers.and.no.31 <- removed_outliers[-31, ]
model.no.outliers.and.no.31 <- lm(density ~., data=data.no.outliers.and.no.31)
t.no.outliers.and.no.31 <- rstudent(model.no.outliers.and.no.31)

summary(model.no.outliers.and.no.31)

for(regressor in names(data.no.outliers.and.no.31)[-which(names(data.no.outliers.and.no.31) == "density")]) {
  plot(data.no.outliers.and.no.31[[regressor]], t.no.outliers.and.no.31, 
       main = paste("Residuals vs", regressor),
       xlab = regressor, 
       ylab = "Residuals")
  abline(h = 0, col = "red") # Add a horizontal line at 0
}

avPlots(data.model)

data.no.outliers.and.no.31.171 <- removed_outliers[-c(31, 171) ]
model.no.outliers.and.no.31.171 <- lm(density ~., data=data.no.outliers.and.no.31.171)
t.no.outliers.and.no.31.171 <- rstudent(model.no.outliers.and.no.31.171)

summary(model.no.outliers.and.no.31.171)

avPlots(model.no.outliers.and.no.31.171)


bc <- boxcox(data.model)
lambda <- bc$x[which.max(bc$y)]
data.model.bc <- lm((data$density^lambda - 1) / lambda ~ ., data = data)
plot(rstudent(data.model.bc) ~ fitted(data.model.bc))


regressors <- names(data)[names(data) != "density"]
par(mfrow = c(4, 4))  # Adjust the layout according to the number of regressors
for (regressor in regressors) {
  plot(data[[regressor]], data.model$residuals,
       main=paste("R-student vs.", regressor), 
       xlab=regressor, ylab="R-student")
  abline(h = 0, lty = "dashed", col = "gray")
}














