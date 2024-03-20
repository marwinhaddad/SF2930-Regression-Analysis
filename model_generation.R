setwd("C:/Users/marwi/Documents/RProjects/Project1")

library(caret)
library(ggplot2)
library(MASS)
library(car)
library(corrplot)
library(boot)
library(glmnet)
library(olsrr)

set.seed(NULL)

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
data <- data[-c(31, 39, 83, 200, 217, 203, 220), ]
fit <- lm(density ~., data=data)
summary(fit)

# VIF values
vif_values <- vif(fit)
print(vif_values)

# condition_indicies
condition_number <- kappa(fit)
print(condition_number)

X <- model.matrix(fit)
WtW <- cor(X[,-1])
WtW
#many strong correlations - chest vs weight, weight vs abdomen, wrist vs neck etc

solve(WtW)  # = (W^T W)^{-1}. Compare diagonal values with VIFs

### Condition number, cutoff value is 100 (MPV p.298)
model.eigen <- eigen(WtW)
condition_number <- max(model.eigen$values) / min(model.eigen$values)
### There is multicollinearity, not servesince condition number < 1000
conditional_indicies <- max(model.eigen$values) / model.eigen$values

# ------------- ALL POSSIBLE REGRESSIONS --------------- #
k <- ols_step_all_possible(fit)

original_data <- convert2SI(read.csv("bodyfatmen.csv"))
m0_all <- ols_step_all_possible(lm(density ~., data = original_data))
m0.original <- m0_all$result[which.max(m0_all$result$n), ]

m1.min_aic <- k$result[which.min(k$result$aic), ]
m2.min_bic <- k$result[which.min(k$result$sbc), ]
# m3.min_cp <- k$result[which.min(k$result$cp), ]
m3.min_adjR2 <- k$result[which.max(k$result$adjr), ]
m4.min_rmse <- k$result[which.min(k$result$rmse), ]

fit.all.rmse <- lm(density ~., data = data)
summary(fit.all.rmse)

fit.all.bic <- lm(density ~ weight + abdomen + wrist, data = data)
summary(fit.all.bic)

fit.all.cp <- lm(density ~ age + weight + neck + abdomen + thigh + forearm + wrist, data = data)
summary(fit.all.cp)

fit.all.adjr2 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, data = data)
summary(fit.all.adjr2)

# ------------- FORWARD SELECTION ----------------#
k.forward.p <- ols_step_forward_p(fit)
k.forward.adj_r2 <- ols_step_forward_adj_r2(fit)
k.forward.aic <- ols_step_forward_aic(fit)
k.forward.sbc <- ols_step_forward_sbc(fit)

fit.forward.aic <- lm(density ~ abdomen + weight + wrist + biceps, data = data)
print(ols_mallows_cp(fit.forward.aic, fit))

# ----------------------- LASSO ----------------------#
y <- as.matrix(data[c(1)])
X <- as.matrix(data[c(-1)])

fit.lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 10)

coef.1se <- coef(fit.lasso, s = "lambda.1se")
coef.min <- coef(fit.lasso, s = "lambda.min")

fit.lasso.lm.1se <- lm(density ~ age + height + abdomen + wrist, data = data)
summary(fit.lasso.lm.1se)
print(coef.1se)

fit.lasso.lm.min <- lm(density ~ age + height + neck + abdomen + thigh + biceps + forearm + wrist, data = data)
summary(fit.lasso.lm.min)
print(coef.min)

# ---------------- Cross validation ----------------- #
train_control <- trainControl(method="repeatedcv", number = 10, repeats = 100)

data.original <- read.csv("bodyfatmen.csv")
data.original <- convert2SI(data.original)

model.0 <- train(density ~., data=data.original, method="lm", trControl=train_control)
model.0

model.1 <- train(density ~ age + weight + neck + abdomen + thigh + forearm + wrist, data = data, method="lm", trControl=train_control)
model.1

model.2 <- train(density ~ weight + abdomen + wrist, data=data, method="lm", trControl=train_control)
model.2

model.3 <- train(density ~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, data=data, method="lm", trControl=train_control)
model.3

model.4 <- train(density ~., data=data, method="lm", trControl=train_control)
model.4

model.5 <- train(density ~ abdomen + weight + wrist + biceps, data=data, method="lm", trControl=train_control)
model.5

model.6 <- train(density ~ age + height + abdomen + wrist, data=data, method="lm", trControl=train_control)
model.6

model.7 <- train(density ~ age + height + neck + abdomen + thigh + biceps + forearm + wrist, data=data, method="lm", trControl=train_control)
model.7

# ------------------ BOTSTRAPPING ------------------- #


models <- list(
  model.0,
  model.1,
  model.2,
  model.3,
  model.4,
  model.5,
  model.6,
  model.7
)

m <- 0
for (model in models) {
  cat(sprintf("model.%d\n", m))
  print(model$call$form)
  
  
  bootstrap_coefficients <- function(data, indicies) {
    boot_sample <- data[indicies, ]
    fit <- lm(model$call$form, data = boot_sample)
    return(coef(fit))
  }
  
  boot_results <- boot(data, bootstrap_coefficients, R = 4000)
  
  for (i in 1:ncol(boot_results$t)) {
    ci <- boot.ci(boot_results, type = "bca", index = i)
    lower_bound <- ci$bca[4]
    upper_bound <- ci$bca[5]
    if (i == 1) {
      name <- "intercept"
    } else {
      name <- model$coefnames[i-1]
    }
    if (lower_bound <= 0 && upper_bound >= 0) {
      cat(sprintf("CI for %s: [%f, %f]*\n", name, lower_bound, upper_bound))
    } else {
      cat(sprintf("CI for %s: [%f, %f]\n", name, lower_bound, upper_bound))
    }
  }
  cat("\n")
  m <- m + 1
}


# ------------------ DETERMINE PRESS -------------------- #

m <- 0
for (model in models) {
  cat(sprintf("model.%d\n", m))
  print(model$call$form)
  
  PRESS <- 0
  for (i in 1:nrow(data)) {
    data_loo <- data[-i, ]
    model_loo <- update(model, data=data_loo)
    pred_loo <- predict(model_loo, newdata = data[i, , drop=FALSE])
    PRESS <- PRESS + (data[i, "density"] - pred_loo)^2
  }
  
  cat(sprintf("PRESS for model.%d: %f\n\n", m, PRESS))
  m <- m + 1
}

# ------------------ MODEL SUMMERIES





