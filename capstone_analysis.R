library(readr)
library(dplyr)
library(tidyverse)
library(glmnet)     # Lasso regression
library(caret)      # Data partitioning
library(ggplot2)    # Visualization
library(broom)      # Model output formatting

# ===================== Step 1: Load and Preprocess Data =====================

# Read the dataset
pnas <- read.csv("/Users/klieeu777/CDS492/latex/pnas/pnas.csv")

# Rename columns for clarity
pnas <- pnas %>%
  rename(
    id = subject_id,
    race = demo_race,
    gender = demo_gender,
    firstgen = demo_firstgen,
    total_sleep = TotalSleepTime,
    data_coverage_nights_frac = frac_nights_with_data,
    term_units_std = Zterm_units_ZofZ,
    bedtime_variability = bedtime_mssd
  )

# Remove rows with missing values and inconsistent data
pnas <- na.omit(pnas)
pnas <- pnas %>% filter(firstgen != 2)

# ===================== Step 2: Create Model Matrix =====================

# Convert data to matrix format for glmnet
x <- model.matrix(term_gpa ~ ., data = pnas)[, -1]  # Drop intercept column
y <- pnas$term_gpa

# ===================== Step 3: Split into Training and Test Sets =====================

set.seed(42)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
x_train <- x[train_idx, ]
y_train <- y[train_idx]
x_test <- x[-train_idx, ]
y_test <- y[-train_idx]

# ===================== Step 4: Fit Lasso Model with Cross-Validation =====================

set.seed(42)
cv_lasso <- cv.glmnet(x_train, y_train, alpha = 1, standardize = TRUE)
best_lambda <- cv_lasso$lambda.min
cat("Best lambda:", best_lambda, "\n")

lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)

# ===================== Step 5: Evaluate on Test Set =====================

y_pred <- predict(lasso_model, s = best_lambda, newx = x_test)

rmse <- sqrt(mean((y_test - y_pred)^2))
r2 <- 1 - sum((y_test - y_pred)^2) / sum((y_test - mean(y_test))^2)

cat("Test RMSE:", round(rmse, 4), "\n")
cat("Test R²:", round(r2, 4), "\n")

# ===================== Step 6: Analyze Non-zero Coefficients =====================

coef_pnas <- tidy(lasso_model) %>%
  filter(term != "(Intercept)" & estimate != 0) %>%
  arrange(desc(abs(estimate)))

print(coef_pnas)

# ===================== Step 7: Visualize Coefficients =====================

ggplot(coef_pnas, aes(x = reorder(term, estimate), y = estimate)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Non-zero Coefficients from Lasso Model",
       x = "Predictor",
       y = "Coefficient Estimate") +
  theme_minimal()

# ===================== Step 8: Visualize Actual vs. Predicted =====================

plot_df <- data.frame(
  Actual = y_test,
  Predicted = as.numeric(y_pred)
)

p <- ggplot(plot_df, aes(x = Actual, y = Predicted)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen", size = 0.8) +
  labs(title = "Actual vs. Predicted Term GPA",
       x = "Actual GPA",
       y = "Predicted GPA") +
  coord_fixed(ratio = 1) +
  theme_minimal()

p
