library(quantmod)
library(goftest)
library(WeightedPortTest)
library(moments)

assets <- c("BZ=F", "EURUSD=X", "CADUSD=X", "GC=F", "SI=F")
getSymbols(assets, src = "yahoo", from = "2012-01-01", to="2025-12-31")
# Only taking adjusted closing 
CADUSD_adjclose <- `CADUSD=X`$`CADUSD=X.Adjusted`
EURUSD_adjclose <- `EURUSD=X`$`EURUSD=X.Adjusted`
BZ_adjclose <- `BZ=F`$`BZ=F.Adjusted`
Gold_adjclose <- `GC=F`$`GC=F.Adjusted`
Silver_adjclose <- `SI=F`$`SI=F.Adjusted`

# Performing outer join
merged_dataset <- merge(BZ_adjclose,EURUSD_adjclose,
                         CADUSD_adjclose,Gold_adjclose, Silver_adjclose)
# Dropping rows containing NA
merged_dataset <- na.omit(merged_dataset)

# Renaming the columns to make it easier to read
colnames(merged_dataset) <- c("Brent", "EURUSD", "CADUSD","Gold", "Silver")
# Now for analysis we will take difference in log-terms
returns <- diff(log(merged_dataset))
returns <- na.omit(returns)

# Converting to data frame so sapply loops column-by-column
returns_df <- as.data.frame(returns)

summary_table <- data.frame(
  Asset    = colnames(returns_df),
  Mean     = sapply(returns_df, mean, na.rm = TRUE),
  Median   = sapply(returns_df, median, na.rm = TRUE),
  Maximum  = sapply(returns_df, max, na.rm = TRUE),
  Minimum  = sapply(returns_df, min, na.rm = TRUE),
  Std_Dev  = sapply(returns_df, sd, na.rm = TRUE),
  Skewness = sapply(returns_df, function(x) skewness(x, na.rm = TRUE)),
  Excess_Kurtosis = sapply(returns_df, function(x) kurtosis(x, na.rm = TRUE) - 3)
)

rownames(summary_table) <- NULL
summary_table[, -1] <- round(summary_table[, -1], 5)

print(summary_table)


png("cad_eur_levels.png", width = 2400, height = 1400, res = 200)

dates <- index(merged_dataset)

plot(dates, merged_dataset$EURUSD, type = "l", col = "#2563eb",
     ylim = range(merged_dataset$EURUSD, merged_dataset$CADUSD, na.rm = TRUE),
     xlab = "", ylab = "Exchange Rate",
     main = "EUR/USD and CAD/USD Exchange Rates (2012-2025)")
lines(dates, merged_dataset$CADUSD, col = "#dc2626")
legend("topright", legend = c("EUR/USD", "CAD/USD"),
       col = c("#2563eb", "#dc2626"), lwd = 1.2, bty = "n")
dev.off()

png("gold_levels.png", width = 2400, height = 1400, res = 200)

plot(dates, merged_dataset$Gold, type = "l", col = "#d4af37",
     xlab = "", ylab = "Price (USD)",
     main = "Gold Futures (2012-2025)")

dev.off()

png("brent_silver_levels.png", width = 2400, height = 1400, res = 200)

plot(dates, merged_dataset$Brent, type = "l", col = "black",
     ylim = range(merged_dataset$Brent, merged_dataset$Silver, na.rm = TRUE),
     xlab = "", ylab = "Price (USD)",
     main = "Brent and Silver Futures (2012-2025)")
lines(dates, merged_dataset$Silver, col = "#94a3b8", lwd = 1.2)
legend("topright", legend = c("Brent", "Silver"),
       col = c("black", "#94a3b8"), lwd = 1.2, bty = "n")

dev.off()


# We will now fit time series models, 

library(rugarch)

# Based on literature search, we will take the following configurations
ar <- c(0, 1, 2, 3, 4)   
ma <- c(0, 1, 2)            
garch <- c("sGARCH", "eGARCH", "gjrGARCH")
dist_of_error <- c("std","sstd", "ged","sged")

# Creating a data frame for containing the optimal time series models
optimal_time_series <- data.frame(
  Asset_Name = character(), Garch_Model  = character(), 
  AR_lags = integer(), MA_lags = integer(),
  Distribution_of_Error   = character(), 
  AIC    = numeric(), Q20_p  = numeric(), Q20_sq_p = numeric(),
  KS_p   = numeric(), BERK_p = numeric(), AD_p = numeric()
)

# Vine copula needs i.i.d uniform margins
# So, I am setting up a matrix for that before running the loop
uniform_margins <- matrix(NA, nrow = nrow(returns), ncol = ncol(returns))
colnames(uniform_margins) <- colnames(returns)

for (i in colnames(returns)) {
  cat("Searching for best time series for asset:", i)
  
  best_aic <- Inf
  optimal_time_series_row <- NULL
  uniform_margins_candidate <- NULL
  
  # Doing a grid search,
  for (g in garch) {
      for (d in dist_of_error) {
        for (p in ar) {
          for (q in ma) {
  
            spec <- ugarchspec(
              variance.model = list(model = g, garchOrder = c(1,1)),
              mean.model     = list(armaOrder = c(p, q), include.mean = TRUE),
              distribution.model = d )
    
            fit <- tryCatch(ugarchfit(spec = spec, data = returns[, i], solver = "hybrid"), error = function(e) NULL)
            
            if (!is.null(fit) && fit@fit$convergence == 0) {
              
              std_res <- as.numeric(residuals(fit, standardize = TRUE))
              
              Ljung_mean <- Weighted.Box.test(std_res, lag = 20, type = "Ljung-Box", fitdf = p+q)$p.value
              Ljung_variance <- Weighted.Box.test(std_res, lag = 20, type = "Ljung-Box", fitdf = 2, sqrd.res = TRUE)$p.value
              
              if (Ljung_mean > 0.05 && Ljung_variance > 0.05) {
                # Transforming to Probability Integral Function
                uniform_vals <- as.numeric(pit(fit))
            
                # K-S test
                KS_pvalue <- ks.test(uniform_vals, "punif")$p.value
                # Anderson-Darling Test
                AD_pvalue <- ad.test(uniform_vals, "punif")$p.value
                
                if (KS_pvalue > 0.05 && AD_pvalue > 0.05) {
                  
                  # Berkowitz Test
                  z_vals <- qnorm(uniform_vals)
                  berk_output <- BerkowitzTest(data = z_vals, lags = 1)
                  Berk_pvalue <- berk_output$LRp
                  
                  if (Berk_pvalue > 0.05) {
                    current_aic <- infocriteria(fit)[1]
                  
                    if (current_aic < best_aic) {
                      best_aic   <- current_aic
                      uniform_margins_candidate <- uniform_vals
                      optimal_time_series_row <- data.frame(
                        Asset_Name = i, 
                        Garch_Model = g,
                        AR_lags = p, 
                        MA_lags = q, 
                        Distribution_of_Error = d, 
                        AIC      = round(current_aic, 4), 
                        Q20_p    = round(Ljung_mean, 4),
                        Q20_sq_p = round(Ljung_variance, 4),
                        KS_p     = round(KS_pvalue, 4), 
                        BERK_p   = round(Berk_pvalue, 4),
                        AD_p    = round(AD_pvalue, 4))
                    }
                  }
                }
            }
          }
        }
      }
    }
  }
    
    if (!is.null(optimal_time_series_row)) {
      optimal_time_series  <- rbind(optimal_time_series, optimal_time_series_row)
      uniform_margins[, i] <- uniform_margins_candidate
    } else {
      cat("No time series model cleared all 4 diagnostic checks for:", i)
}}

print(optimal_time_series) # Saw the optimal output

# Now for latex I will get an output,

# Brent
spec_brent <- ugarchspec(
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  distribution.model = "sstd"
)

# EURUSD
spec_eurusd <- ugarchspec(
  mean.model = list(armaOrder = c(4, 2), include.mean = TRUE),
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  distribution.model = "sstd"
)

# CADUSD
spec_cadusd <- ugarchspec(
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  distribution.model = "ged"
)

# Gold
spec_gold <- ugarchspec(
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  distribution.model = "std"
)

# Silver
spec_silver <- ugarchspec(
  mean.model = list(armaOrder = c(4, 2), include.mean = TRUE),
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  distribution.model = "ged"
)

fit_brent  <- ugarchfit(spec = spec_brent,  data = returns[, "Brent"])
fit_eurusd <- ugarchfit(spec = spec_eurusd, data = returns[, "EURUSD"])
fit_cadusd <- ugarchfit(spec = spec_cadusd, data = returns[, "CADUSD"])
fit_gold   <- ugarchfit(spec = spec_gold,   data = returns[, "Gold"])
fit_silver <- ugarchfit(spec = spec_silver, data = returns[, "Silver"])

get_values <- function(fit_object) {
  coef_matrix <- fit_object@fit$matcoef
  
  results <- data.frame(
    Parameter = rownames(coef_matrix),
    Estimate  = round(coef_matrix[, 1], 6),
    Std_Error = round(coef_matrix[, 2], 6),
    P_Value   = round(coef_matrix[, 4], 4)
  )

  return(results)
}

# Run the function on each of your fitted models:
print("BRENT")
print(get_values(fit_brent))

print("EURUSD")
print(get_values(fit_eurusd))

print("CADUSD")
print(get_values(fit_cadusd))

print("GOLD")
print(get_values(fit_gold))

print("SILVER")
print(get_values(fit_silver))



library(VineCopula)
R_Vine <- RVineStructureSelect(
  data = uniform_margins, 
  familyset = c(2, 3, 4, 6, 7, 8, 9, 10,      # Base families
                13, 14, 16, 17, 18, 19, 20,   # 180-degree rotated
                23, 24, 26, 27, 28, 29, 30,   # 90-degree rotated
                33, 34, 36, 37, 38, 39, 40),  # 270-degree rotated
  type = 0,
  trunclevel = 2,
  progress = FALSE
)

# Shows edges where a copula exists
edge_index <- which(R_Vine$family > 0)
# Calculates what copula families are there at each edge
edge_families <- R_Vine$family[edge_index]
# Stores copula family number with 2 parameters
two_param_families <- c(2, 7, 8, 9, 10, 17, 18, 19, 20, 27, 28, 29, 30, 37, 38, 39, 40)
# Shows edges where a two parameter copula exists
edge_index_two_param <- which(edge_families %in% two_param_families)
  
# Before running optimiser, we will enforce bounds on copula family parameters,
# Some of these bounds are in accordance with the bounds of BiCopEst function

par1_bounds <- function(a) {
  # Gaussian is 1 & Student-t is 2
  if (a %in% c(1, 2))      return(c(-0.999, 0.999))
  # Clayton is (3, 13) & Rotated Clayton is (23, 33)
  if (a %in% c(3, 13))     return(c(0.001, 17))       
  if (a %in% c(23, 33))    return(c(-17, -0.001))     
  # Gumbel is (4, 14) & Rotated Gumbel is (24, 34)
  if (a %in% c(4, 14))     return(c(1.001, 15))       
  if (a %in% c(24, 34))    return(c(-15, -1.001))   
  # Joe is (6, 16) & Rotated Joe is (26, 36)
  if (a %in% c(6, 16))     return(c(1.001, 15))       
  if (a %in% c(26, 36))    return(c(-15, -1.001))     
  # BB1 is (7, 17) & Rotated BB1 is (27, 37)
  if (a %in% c(7, 17))     return(c(0.001, 5))       
  if (a %in% c(27, 37))    return(c(-5, -0.001))     
  # BB6 is (8, 18) & Rotated BB6 is (28, 38)
  if (a %in% c(8, 18))     return(c(1.001, 6))       
  if (a %in% c(28, 38))    return(c(-6, -1.001))     
  # BB7 is (9, 19) & Rotated BB7 is (29, 39)
  if (a %in% c(9, 19))     return(c(1.001, 5))       
  if (a %in% c(29, 39))    return(c(-5, -1.001))     
  # BB8 is (10, 20) & Rotated BB8 is (30, 40)
  if (a %in% c(10, 20))    return(c(1.001, 6))       
  if (a %in% c(30, 40))    return(c(-6, -1.001))
}

par2_bounds <- function(a) {
  if (a == 2)              return(c(2.001, 30))       
  if (a %in% c(7, 17))     return(c(1.001, 6))       
  if (a %in% c(27, 37))    return(c(-6, -1.001))     
  if (a %in% c(8, 18))     return(c(1.001, 6))       
  if (a %in% c(28, 38))    return(c(-6, -1.001))     
  if (a %in% c(9, 19))     return(c(0.001, 6))       
  if (a %in% c(29, 39))    return(c(-6, -0.001))     
  if (a %in% c(10, 20))    return(c(0.001, 0.999))    
  if (a %in% c(30, 40))    return(c(-0.999, -0.001))   
}

# Creating bounds based on fitted R vine
parameter1_bounds <- t(sapply(edge_families, par1_bounds))            
parameter2_bounds <- t(sapply(edge_families[edge_index_two_param], par2_bounds))

number_of_parameter2 <- length(edge_index_two_param)
number_of_edges <- length(edge_index)

# First parameter1_bounds and then parameter2_bounds in a 1D vector for optimiser
overall_lowerbound <- c(parameter1_bounds[, 1], if (number_of_parameter2 > 0) parameter2_bounds[, 1] else numeric(0))
overall_upperbound <- c(parameter1_bounds[, 2], if (number_of_parameter2 > 0) parameter2_bounds[, 2] else numeric(0))

EM_Rvine <- function(theta, R_Vine, edge_index , edge_index_two_param , number_of_edges, number_of_parameter2) {
      NewR_Vine <- R_Vine
      NewR_Vine$par[edge_index] <- theta[1:number_of_edges]
      if (number_of_parameter2 > 0) {
        b <- (number_of_edges + 1):(number_of_edges + number_of_parameter2)
        NewR_Vine$par2[edge_index[edge_index_two_param]] <- theta[b]
      }
      NewR_Vine} # Fits R-vine with new parameters theta but existing copula families and 

Calculation_of_NegLogLik <- function(theta, weights, uniform_margins, R_Vine,
                                      edge_index, edge_index_two_param,
                                      number_of_edges, number_of_parameter2) {
  
  RVM <- EM_Rvine(theta, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2)
  ll <-RVineLogLik(uniform_margins, RVM, separate = TRUE, calculate.V = FALSE)$loglik
  return(-sum(weights * ll))} # Gives back the negative of weighted log likelihood of fitted R-vine
  
EM_Optimiser <- function(uniform_margins, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2,
                         overall_lowerbound, overall_upperbound, P11, P22, theta_R1, theta_R2,
                         tol = 1e-4, max_iter=100) {
  
  T <- nrow(uniform_margins)
  ll_prev <- -Inf
  
  for (n in 1:max_iter) {
    # This fits R_Vine to Regime 1 (R1) and Regime 2 (R2)
    NewR_Vine_R1 <- EM_Rvine(theta_R1, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2)
    NewR_Vine_R2 <- EM_Rvine(theta_R2, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2)
    
    # Calculates the likelihood of R vines of Regime 1 and 2
    like_R1 <- exp(RVineLogLik(uniform_margins, NewR_Vine_R1, separate = TRUE)$loglik)
    like_R2 <- exp(RVineLogLik(uniform_margins, NewR_Vine_R2, separate = TRUE)$loglik)
    
    # Stationary probabilities for Regime 1 and Regime 2
    P_R1 <- (1 - P22) / (2 - P11 - P22) 
    P_R2 <- 1 - P_R1
    
    # Predicted Probabilities
    Predicted_R1 <- numeric(T) 
    Predicted_R2 <- numeric(T)
    # Filtered Probabilities
    Filtered_R1 <- numeric(T)
    Filtered_R2 <- numeric(T) 

    Accumulated_loglike  <- 0
    
    for (t in 1:T) {
      # In t=1 we use stationary probabilities for predicted probabilities
      if (t == 1) {
        temporary_R1 <- P_R1
        temporary_R2 <- P_R2
        
      } else { # This is Prediction Step for t>1
        
        temporary_R1 <- P11 * Filtered_R1[t - 1] + (1 - P22) * Filtered_R2[t - 1]
        temporary_R2 <- (1 - P11) * Filtered_R1[t - 1] + P22 * Filtered_R2[t - 1]
      }
      
      # temporary variables change every iteration but predicted_R stores
      Predicted_R1[t] <- temporary_R1
      Predicted_R2[t] <- temporary_R2
      
      # Denominator_filtering shows the denominator in filtering of S_t equation
      Denominator_filtering <- temporary_R1 * like_R1[t] + temporary_R2 * like_R2[t]
      Accumulated_loglike <- Accumulated_loglike + log(Denominator_filtering) 
      
      # Filtering step
      Filtered_R1[t] <- temporary_R1 * like_R1[t] / Denominator_filtering
      Filtered_R2[t] <- temporary_R2 * like_R2[t] / Denominator_filtering
      }
    
    # Now we will implement backward smoothing algorithm
    Smooth_R1 <- numeric(T)
    Smooth_R2 <- numeric(T)
    
    # Initialising the smoothing algorithm in t=T
    Smooth_R1[T] <- Filtered_R1[T] 
    Smooth_R2[T] <- Filtered_R2[T]
    
    # To calculate Kim's smoother, we will break the formula into 4 distinct parts
    
    part11 <- numeric(T - 1)
    part12 <- numeric(T - 1)
    part21 <- numeric(T - 1)
    part22 <- numeric(T - 1)
    
    for (t in (T - 1):1) {
      part11[t] <- (Filtered_R1[t] * P11       * Smooth_R1[t + 1]) / Predicted_R1[t + 1]
      part12[t] <- (Filtered_R1[t] * (1 - P11) * Smooth_R2[t + 1]) / Predicted_R2[t + 1]
      part21[t] <- (Filtered_R2[t] * (1 - P22) * Smooth_R1[t + 1]) / Predicted_R1[t + 1]
      part22[t] <- (Filtered_R2[t] * P22       * Smooth_R2[t + 1]) / Predicted_R2[t + 1]
      
      Smooth_R1[t] <- part11[t] + part12[t]
      Smooth_R2[t] <- part21[t] + part22[t]
    }
    
    # Now before we go to the M step, if E step doesn't lead to much improvement in 
    # log likelihood between 2 consecutive iterations we can terminate the process
    # It means convergence 
    
    if (abs(Accumulated_loglike - ll_prev) < tol) {
      cat("EM converged.\n")
      break
    }
    
    ll_prev <- Accumulated_loglike
    
    # Now we will proceed to M step,
    # We will optimise P11 and P22 first
    P11 <- sum(part11) / sum(Smooth_R1[1:(T - 1)])
    P22 <- sum(part22) / sum(Smooth_R2[1:(T - 1)])
    # Now we will optimise theta_R1 and theta_R2 (parameters of copulas in regime 1 and 2)
    fit_R1 <- nlminb(
      start = theta_R1,
      objective = Calculation_of_NegLogLik,
      lower = overall_lowerbound,
      upper = overall_upperbound,
      weights = Smooth_R1,
      uniform_margins = uniform_margins,
      R_Vine = R_Vine,
      edge_index = edge_index,
      edge_index_two_param = edge_index_two_param,
      number_of_edges = number_of_edges,
      number_of_parameter2 = number_of_parameter2,
      control = list(trace = 0, iter.max = 60, eval.max = 100)
    )
    theta_R1 <- fit_R1$par
    
    fit_R2 <- nlminb(
      start = theta_R2,
      objective = Calculation_of_NegLogLik,
      lower = overall_lowerbound,
      upper = overall_upperbound,
      weights = Smooth_R2,
      uniform_margins = uniform_margins,
      R_Vine = R_Vine,
      edge_index = edge_index,
      edge_index_two_param = edge_index_two_param,
      number_of_edges = number_of_edges,
      number_of_parameter2 = number_of_parameter2,
      control = list(trace = 0, iter.max = 60, eval.max = 100)
    )
    theta_R2 <- fit_R2$par
  }
  
  # Now at the end, when EM is done running, it will return this list
  list(
    P11 = P11, P22 = P22,
    theta_R1 = theta_R1, theta_R2 = theta_R2,
    loglikelihood_value = Accumulated_loglike, iterations = n,
    Predicted_R1 = Predicted_R1, Predicted_R2 = Predicted_R2,
    Filtered_R1 = Filtered_R1, Filtered_R2 = Filtered_R2,
    Smooth_R1 = Smooth_R1, Smooth_R2 = Smooth_R2,
    NewR_Vine_R1 = EM_Rvine(theta_R1, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2),
    NewR_Vine_R2 = EM_Rvine(theta_R2, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2)
  )
}
  
# Configuring multiple starts 

number_of_starts <- 10
best_em  <- NULL
best_log_likelihood  <- -Inf

set.seed(42)

print("Running 10 iterations of EM")

library(readr)
fsi_data <- read_csv("fsi.csv")
fsi_data <- fsi_data[, c("Date", "OFR FSI")]

returns_dates <- as.Date(index(returns))
fsi_data$Date <- as.Date(fsi_data$Date)

# Matching common dates between returns and FSI
match_index <- match(returns_dates, fsi_data$Date)
valid_mask  <- !is.na(match_index)

uniform_margins_aligned <- uniform_margins[valid_mask, ]
fsi_aligned <- fsi_data[match_index[valid_mask], ]

calm_days <- which(fsi_aligned$`OFR FSI` <= 0)
crisis_days <- which(fsi_aligned$`OFR FSI` > 0)

uniform_margins_calm   <- uniform_margins_aligned[calm_days, ]
uniform_margins_crisis <- uniform_margins_aligned[crisis_days, ]

fit_calm <- RVineSeqEst(uniform_margins_calm, R_Vine)
fit_crisis <- RVineSeqEst(uniform_margins_crisis, R_Vine)

# After fitting, we will extract parameters
extract <- function(rvm, edge_index, edge_index_two_param) {
  par1 <- rvm$par[edge_index]
  par2 <- rvm$par2[edge_index[edge_index_two_param]]
  return(c(par1, par2))
}

theta_R1_from_fsi <- extract(fit_calm, edge_index, edge_index_two_param)
theta_R2_from_fsi <- extract(fit_crisis, edge_index, edge_index_two_param)

# Before feeding this to EM algorithm, we need to ensure these are within bounds

theta_R1_from_fsi <- pmin(pmax(theta_R1_from_fsi, overall_lowerbound), overall_upperbound)
theta_R2_from_fsi <- pmin(pmax(theta_R2_from_fsi, overall_lowerbound), overall_upperbound)

states_from_fsi <- ifelse(fsi_aligned$`OFR FSI` > 0, 2, 1)

yesterday_state <- states_from_fsi[-length(states_from_fsi)] 
today_state <- states_from_fsi[-1]

p11_from_fsi <- mean(today_state[yesterday_state == 1] == 1)
p22_from_fsi <- mean(today_state[yesterday_state == 2] == 2)

for (s in 1:number_of_starts) {
  
  if (s == 1) {
    P11 <- p11_from_fsi; P22 <- p22_from_fsi  
    theta_R1 <- theta_R1_from_fsi
    theta_R2 <- theta_R2_from_fsi
  } else {
    
    noise1 <- runif(length(theta_R1), -0.1, 0.1) * (overall_upperbound - overall_lowerbound)
    noise2 <- runif(length(theta_R2), -0.1, 0.1) * (overall_upperbound - overall_lowerbound)
    
    theta_R1 <- pmin(pmax(theta_R1 + noise1, overall_lowerbound), overall_upperbound)
    theta_R2 <- pmin(pmax(theta_R2 + noise2, overall_lowerbound), overall_upperbound)

    P11 <- runif(1, 0.80, 0.95)
    P22 <- runif(1, 0.80, 0.95)
  }
  
  EM_result <- EM_Optimiser(uniform_margins, R_Vine, edge_index, edge_index_two_param, number_of_edges, number_of_parameter2,
                            overall_lowerbound, overall_upperbound,
                            P11, P22, theta_R1, theta_R2,
                            tol = 1e-4, max_iter = 50)
  
  cat("Running iteration", s, "\n")
  cat("Parameter R1 starting values",theta_R1)
  cat("Parameter R2 starting values",theta_R2)
  print(EM_result$loglikelihood_value)
  
  if (!is.null(EM_result) && is.finite(EM_result$loglikelihood_value) && EM_result$loglikelihood_value > best_log_likelihood) {
    best_log_likelihood <- EM_result$loglikelihood_value
    best_em <- EM_result
  }
}
  
cat("Best Log-likelihood", best_log_likelihood)
cat("EM Iterations",best_em$iterations)
print(best_em)
  
# Now we have to compute Kendall's Tau as it was never computed during EM
  
Function_tau <- function(rvm, edge_index) {
  family_name  <- rvm$family[edge_index]
  par1 <- rvm$par[edge_index]
  par2 <- rvm$par2[edge_index]
  mapply(function(f, p1, p2) BiCopPar2Tau(family = f, par = p1, par2 = p2),
         family_name, par1, par2)
}

tau_R1 <- Function_tau(best_em$NewR_Vine_R1, edge_index)
tau_R2 <- Function_tau(best_em$NewR_Vine_R2, edge_index)

print("Regime 1 Tau")
print(round(tau_R1, 4))
print("Regime 2 Tau")
print(round(tau_R2, 4))

# Now we have to see which R_Vine is actually the calm and crisis one
# as during the 

mean_tau_R1 <- mean(abs(tau_R1))
mean_tau_R2 <- mean(abs(tau_R2))

cat("Fitted R1 Tree", round(mean_tau_R1, 4), "\n")
cat("Fitted R2 Tree", round(mean_tau_R2, 4), "\n")

# Final output of the tree structure 
tau_matrix_R1 <- matrix(0, nrow(best_em$NewR_Vine_R1$family), ncol(best_em$NewR_Vine_R1$family))
tau_matrix_R1[edge_index] <- tau_R1
best_em$NewR_Vine_R1$tau <- tau_matrix_R1

tau_matrix_R2 <- matrix(0, nrow(best_em$NewR_Vine_R2$family), ncol(best_em$NewR_Vine_R2$family))
tau_matrix_R2[edge_index] <- tau_R2
best_em$NewR_Vine_R2$tau <- tau_matrix_R2




