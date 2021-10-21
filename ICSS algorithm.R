# ICSS
# ---------------------------------------------------------------------------------------------------------------------------
# ---Description---
# It returns a list containing ICSS's computed results for the input series.
#
# ---Usage---
# ICSS_final(series, x_time = NULL, critical_value = 1.358, demean = F, write = F)
#
#---Arguments---
# series: a univariate series for detecting structure breaks.
# x_time: input a customized series for the real date, default setting is NULL.
# critical_value: Critical value of the statistic is optional. The default critical value is set to be 1.358 (95%) according to Table 1 in Inclan & Tiao (1994).
#                 You can input 1.628 (99%), 1.358 (95%, default setting), 1.224 (90%), and 1.019 (75%). For other available critical values in detailed, please refer to Table 1 in Inclan & Tiao (1994).
# demean: filter the return series with a constant mean value, default setting is F.
# write: output the ICSS's results into a ".csv" file named "ICSS_result.csv", default setting is F.
# ---------------------------------------------------------------------------------------------------------------------------

# ICSS_plot
# ---------------------------------------------------------------------------------------------------------------------------
# ---Description---
# It plots the return series with the mean and standard deviation of periods classified by ICSS's breakpoints.
#
# ---Usage---
# ICSS_plot(ICSS_object, return_extra_interval = 2, sd_times = 1, period_mean = F, double_y_axes = F, x_time = NULL, col_series = "grey", col_sd = "black", lty_sd = 1, col_mean = "black", lty_mean = 1, x_name = NULL, y_left_name = NULL, y_right_name = NULL, title = "Returns with standard deviations")
#
#---Arguments---
# ICSS_object: the object containing the result obtained from the ICSS_final function.
# return_extra_interval: specify the extra y-axis interval for the return plotting, default setting is 2.
# sd_times: the times of standard deviation for plotting, default setting is 1.
# period_mean: plot the periods' mean, default setting is F (a horizontal line of 0). 
# double_y_axes: plot with double y-axis, default setting is F. If the left y-axis is displayed by standard deviation (returns), the right y-axis is displayed by returns (standard deviations).
# x_time: input a customized series for x-axis if you want to plot the real date.
# col_series: the color of the return series, default setting is "grey".
# col_sd: the color of the standard deviation's plotting, default setting is "black".
# lty_sd: the line type of the standard deviation's plotting, default setting is 1 (solid line).
# col_mean: the color of the mean's plotting, default setting is "black".
# lty_mean: the line type of the mean's plotting, default setting is 1 (solid line).
# x_name: x-axis's label, default setting is NULL.
# y_left_name: left y-axis's label, default setting is NULL.
# y_right_name: right y-axis's label, default setting is NULL.
# title: title of the plot.
# ---------------------------------------------------------------------------------------------------------------------------

# ICSS_dummy
# ---------------------------------------------------------------------------------------------------------------------------
# ---Description---
# It returns a matrix of dummy variables.
#
# ---Usage---
# ICSS_dummy(ICSS_object, write = F)
#
#---Arguments---
# ICSS_object: the object containing the result obtained from the ICSS_final function.
# write: output the dummy variable into a ".csv" file named "ICSS_dummy.csv", default setting is F.
# ---------------------------------------------------------------------------------------------------------------------------



# (2) Source code
# ---------------------------------------------------------------------------------------------------------------------------
# ---1."CCSS" function---
# It returns the position of a structural break for a typical series "Y".
CCSS <- function(Y, critical_value){ # Refer to the formula from Inclan & Tiao (1994).
  ct <- 0
  l <- length(Y)
  ck <- matrix(0,l,1)
  dk <- matrix(0,l,1)
  sy <- var(Y)
  a <- matrix(0,l,1)
  
  for(i in 1:l){
    a[i] <- (Y[i] - (sum(Y[1:i])/i))/sqrt((1+1/i)*sy)
  }
  
  for(i in 1:l){
    ct <- ct + a[i]^2
  }
  
  ck[1] <- a[1]^2
  
  for(i in 2:l){
    ck[i] <- ck[i-1]+a[i]^2
  }
  
  for(i in 1:(l-1)){
    dk[i] <- (ck[i]/ct)-(i/l)
  }
  
  dk[l] <- 0
  
  maxk <- max(sqrt(l/2)*abs(dk))
  kstar <- which.max(sqrt(1/2)*abs(dk)) 
  
  if(maxk <= critical_value){
    kstar <- -1
  }
  
  return(kstar) # "kstar" here will return the position of a structural break. If there is no structural break detected, "-1" will be returned.
}



# ---2."ICSS" function---
# It returns all the positions of structural break for "series" by self-calling iteration.
ICSS <- function(Y, series, t1, t2, critical_value){ # "Y" is the full sample series. "series" is the series for detecting structure breaks. 
  k <- CCSS(series, critical_value)                  # "Y" is fixed and "series" is variable. Please note that "Y" is the same to "series" before the first iteration.
  t3 <- t1+k-1 # "t3" is the position of a structural break returned from "ICSS" function that called by each iteration.
  
  if(k != -1){
    store_t3 <<- rbind(store_t3, t3)  # "store_t3" is a pre-defined global variable used to store "t3".
    ICSS(Y,Y[t1:t3],t1,t3,critical_value)
    ICSS(Y,Y[t3:t2],t3,t2,critical_value)
  }
}



# ---3."ICSS_final"---
# It is the main function that calls "CCSS" and "ICSS" functions. It returns a numeric result containing breakpoints positions in ascending order.
ICSS_final <- function(series, x_time = NULL, critical_value = 1.358, demean = F, write = F){  # The default critical value is set to be 1.358 (95%) according to Table 1 in Inclan & Tiao (1994).
  if(ncol(as.matrix(series)) != 1){
    stop("Please input the univariate series.") # Only support the univariate case.
  }
  
  if(is.null(x_time) == F){
    if(length(as.matrix(series)) != length(as.matrix(x_time))){
      stop("Please make sure that the series and x_time maintain the same length.")
    }
  }
  
  if(demean == T){
    series <- as.matrix(series)-mean(as.matrix(series))
  }
  
  store_t3 <<- NULL # Initialize "store_t3" (defined as a global variable) as NULL to prepare for storage for each calling of "ICSS_final" function.
  
  Y <- as.matrix(series)
  series <- as.matrix(series)
  
  t1 <- 1
  t2 <- length(series)
  
  ICSS(Y, series, t1, t2, critical_value)
  
  breakpoint <- sort(as.numeric(store_t3)) # Sort the identified structural breakpoints positions in ascending order.
  num_breakpoint <- length(breakpoint)
  
  
  period_store <- list()
  period_store <- c(period_store, list(series[1:(breakpoint[1]-1)]))
  for(i in c(1:(length(breakpoint)-1))){
    period_store <- c(period_store, list(series[breakpoint[i]:(breakpoint[i+1]-1)]))
  }
  period_store <- c(period_store, list(series[breakpoint[length(breakpoint)]:length(series)]))
  

  observation_store <- NULL
  sd_store <- NULL
  mean_store <- NULL
  for(i in c(1:length(period_store))){
    observation_store <- c(observation_store, length(period_store[[i]]))
    sd_store <- c(sd_store, sd(period_store[[i]]))
    mean_store <- c(mean_store, mean(period_store[[i]]))
  }
  

  if(is.null(x_time) == F){
    x_time <- as.character(x_time)
    periodname_store <- NULL
    periodname_store <- c(periodname_store, paste("period_1", " (1 - ", (breakpoint[1]-1), ")", " (", x_time[1], " - ", x_time[breakpoint[1]-1], ")", sep = ""))
    for(i in c(1:(length(breakpoint)-1))){
      periodname_store <- c(periodname_store, paste("period_", 1+i, " (", breakpoint[i], " - ", (breakpoint[i+1]-1), ")", " (", x_time[breakpoint[i]], " - ", x_time[breakpoint[i+1]-1], ")", sep = ""))
    }
    periodname_store <- c(periodname_store, paste("period_", length(breakpoint)+1, " (", breakpoint[length(breakpoint)], " - ", length(series), ")", " (", x_time[breakpoint[length(breakpoint)]], " - ", x_time[length(series)], ")", sep = ""))
  }
  else{
    periodname_store <- NULL
    periodname_store <- c(periodname_store, paste("period_1", " (1 - ", (breakpoint[1]-1), ")", sep = ""))
    for(i in c(1:(length(breakpoint)-1))){
      periodname_store <- c(periodname_store, paste("period_", 1+i, " (", breakpoint[i], " - ", (breakpoint[i+1]-1), ")", sep = ""))
    }
    periodname_store <- c(periodname_store, paste("period_", length(breakpoint)+1, " (", breakpoint[length(breakpoint)], " - ", length(series), ")", sep = ""))
  }
  
  
  result <- t(rbind(observation_store, sprintf("%0.4f",mean_store), sprintf("%0.4f",sd_store)))
  rownames(result) <- periodname_store
  colnames(result) <- c("Period's observations", "Mean", "Standard deviation")
  
  
  output <- list(series = series, result = noquote(result), num_breakpoint = num_breakpoint, breakpoint = breakpoint, num_period = num_breakpoint+1, period_observation = observation_store, period_mean = mean_store, period_sd = sd_store)
  
  cat(paste("The total observations of the input series is:", length(series)))
  cat("\n")
  cat(paste("The number of breakpoints detected by ICSS is:", length(breakpoint)))
  cat("\n")
  cat("The position of breakpoints detected by ICSS is: ")
  cat(breakpoint)
  cat("\n")
  cat(paste("The number of total periods is:", length(breakpoint)+1))
  cat("\n")
  cat("The observation and respective standard deviation of the period are:")
  cat("\n")
  cat("\n")
  print(noquote(result))
  cat("\n")
  
  if(write == T){
    write.csv(x = noquote(result), file = "ICSS_result.csv")
  }
  
  return(output)
}



# ---4."ICSS_plot"---
# It plots the return series with the mean and standard deviation of periods classified by ICSS's breakpoints.
ICSS_plot <- function(ICSS_object, return_extra_interval = 2, sd_times = 1, period_mean = F, double_y_axes = F, x_time = NULL, col_series = "grey", col_sd = "black", lty_sd = 1, col_mean = "black", lty_mean = 1, x_name = NULL, y_left_name = NULL, y_right_name = NULL, title = "Returns with the mean and standard deviation"){
  
  if(is.null(x_time) == F){
    if(length(as.matrix(ICSS_object$series)) != length(as.matrix(x_time))){
      stop("Please make sure that the series and x_time maintain the same length.")
    }
  }
  
  y_sd <- NULL
  for(i in c(1:length(ICSS_object$period_sd))){
    y_sd <- c(y_sd, rep(ICSS_object$period_sd[i],ICSS_object$period_observation[i]))
  }
  
  if(period_mean == F){
    h_mean <- 0
  }
  else if(period_mean == T){
    h_mean <- NULL
    for(i in c(1:length(ICSS_object$period_mean))){
      h_mean <- c(h_mean, rep(ICSS_object$period_mean[i],ICSS_object$period_observation[i]))
    }
  }
  
  range <- max(abs(ICSS_object$series))+return_extra_interval
  y_range <- c(-range,range)
  
  if(is.null(x_time) == T){
    x_series <- c(1:length(ICSS_object$series))
  }
  else{
    x_series <- x_time
  }

  xname <- "Date"
  yname <- "Return"
  yname2 <- "Standard deviation"
  
  if(is.null(y_left_name) == F){
    yname <- y_left_name
  }
  
  if(is.null(y_right_name) == F){
    yname2 <- y_right_name
  }
  
  if(is.null(x_name) == F){
    xname <- x_name
  }
  
  if(double_y_axes == T){
    par(mar=c(5,4,4,4.3)+0.1)
  }
  
  plot(x = x_series, y = ICSS_object$series, type = "l", col = col_series, xlab = xname, ylab = yname, ylim = y_range, main = title)
  
  if(period_mean == F){
    abline(h = h_mean, lty = lty_mean, col = col_mean)
  }
  else if(period_mean == T){
    lines(x = x_series, y = h_mean, lty = lty_mean, col = col_mean)
  }
  
  par(new=TRUE)
  plot(x = x_series, y = y_sd*sd_times, type = "l", lty = lty_sd, col = col_sd, xlab = "", ylab = "", ylim = y_range)
  lines(x = x_series, y = -y_sd*sd_times, lty = lty_sd, col = col_sd)
  
  if(double_y_axes == T){
    axis(4)
    mtext(yname2, side = 4, line = 3)
  }
  
  if(double_y_axes == T){
    par(mar=c(5,4,4,2)+0.1)
  }
}



# ---5."ICSS_dummy"---
# It is a function returning a matrix of dummy variables.
ICSS_dummy <- function(ICSS_object, write = F){
  breakpoint <-  ICSS_object$breakpoint
  period_observation <- ICSS_object$period_observation
  
  period_begin <- c(1, breakpoint)
  period_end <- period_begin-1+period_observation
  
  dummy <- matrix(0,length(ICSS_object$series), length(period_begin))
  dummy_name <- NULL
  for(i in 1:ncol(dummy)){
    dummy[period_begin[i]:period_end[i],i] <- 1
    dummy_name <- c(dummy_name, paste("dummy variable", i))
  }
  
  colnames(dummy) <- dummy_name
  
  if(write == T){
    write.csv(x = dummy, file = "ICSS_dummy.csv")
  }
  
  return(dummy)
}
# ---------------------------------------------------------------------------------------------------------------------------



# (3) Example
# ---------------------------------------------------------------------------------------------------------------------------
library(rugarch)
data(dmbp) # Example data from "rugarch" package.




ICSS_1 <- ICSS_final(dmbp[,1], write = T) # (default critical value = 1.358 (95%))
                               # "dmbp[,1]" is a log return series of the bilateral Deutschemark/British pound rate constructed from the corresponding U.S. dollar rates (The Bollerslev-Ghysel benchmark dataset)
ICSS_1$result # ICSS result (with default critical value = 1.358 (95%))
ICSS_1$num_breakpoint # number of breakpoints (with default critical value = 1.358 (95%))
ICSS_1$breakpoint # breakpoints (with default critical value = 1.358 (95%))
ICSS_1$num_period # number of periods (with default critical value = 1.358 (95%))
ICSS_1$period_observation # periods' observations (with default critical value = 1.358 (95%))
ICSS_1$period_sd # periods' standard deviation (with default critical value = 1.358 (95%))
ICSS_1$period_mean # periods' mean (with default critical value = 1.358 (95%))
ICSS_1_dummy <- ICSS_dummy(ICSS_1, write = T) # dummy variable (with default critical value = 1.358 (95%))
ICSS_plot(ICSS_1, return_extra_interval = 2, sd_times = 3, period_mean = T, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.358 (95%))") # single axis (left y-axis), plotting periods' mean, with plus and minus 3 times standard deviation
ICSS_plot(ICSS_1, return_extra_interval = 2, sd_times = 1, period_mean = F, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.358 (95%))") # single axis (left y-axis), a horizontal line of 0 instead of plotting periods' mean, with plus and minus 1 times standard deviation
ICSS_plot(ICSS_1, return_extra_interval = 2, sd_times = 3, period_mean = T, double_y_axes = T, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.358 (95%))") # double y-axis (left and right y-axis), plotting periods' mean, with plus and minus 3 times standard deviation

plot(x = c(1:length(ICSS_1$series)), y = ICSS_1$series, type = "l", xlab = "Date", ylab = "Percentage", main = "Deutschemark/British pound exchange rate returns (critical value = 1.358 (95%))")
abline(v = c(1:length(ICSS_1$series))[ICSS_1$breakpoint], lty = 1, col = "red")




ICSS_2 <- ICSS_final(dmbp[,1], critical_value = 1.628) # (critical value = 1.628 (99%))

ICSS_2$result # ICSS result (with critical value = 1.628 (99%))
ICSS_2$num_breakpoint # number of breakpoints (with critical value = 1.628 (99%))
ICSS_2$breakpoint # breakpoints (with critical value = 1.628 (99%))
ICSS_2$num_period # number of periods (with critical value = 1.628 (99%))
ICSS_2$period_observation # periods' observations (with critical value = 1.628 (99%))
ICSS_2$period_sd # periods' standard deviation (with critical value = 1.628 (99%))
ICSS_2$period_mean # periods' mean (with critical value = 1.628 (99%))
ICSS_2_dummy <- ICSS_dummy(ICSS_2) # dummy variable (with critical value = 1.628 (99%))
ICSS_plot(ICSS_2, return_extra_interval = 2, sd_times = 3, period_mean = T, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.628 (99%))") # single axis (left y-axis), plotting periods' mean, with plus and minus 3 times standard deviation
ICSS_plot(ICSS_2, return_extra_interval = 2, sd_times = 1, period_mean = F, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.628 (99%))") # single axis (left y-axis), a horizontal line of 0 instead of plotting periods' mean, with plus and minus 1 times standard deviation
ICSS_plot(ICSS_2, return_extra_interval = 2, sd_times = 3, period_mean = T, double_y_axes = T, col_sd = "red", lty_sd = 2, title = "Deutschemark/British pound exchange rate returns with standard deviations (critical value = 1.628 (99%))") # double y-axis (left and right y-axis), plotting periods' mean, with plus and minus 3 times standard deviation

plot(x = c(1:length(ICSS_2$series)), y = ICSS_2$series, type = "l", xlab = "Date", ylab = "Percentage", main = "Deutschemark/British pound exchange rate returns (critical value = 1.628 (99%))")
abline(v = c(1:length(ICSS_2$series))[ICSS_2$breakpoint], lty = 1, col = "red")
# ---------------------------------------------------------------------------------------------------------------------------