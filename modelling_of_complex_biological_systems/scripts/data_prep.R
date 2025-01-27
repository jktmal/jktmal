library(tidyverse)
library(lubridate)

setwd("/home/jasiek/Dropbox/mcbs/project/")

# moving average function
ma <- function(x,n){#n = number to average over on each tail
  L <- length(x)
  res = x
  for(k in 1:L){
    i1 <- max(c(1,k-n)); i2 <- min(c(L,k+n));
    res[k] = mean(x[i1:i2])
  }
  return(res)
}

DF <- read_csv("https://covidtracking.com/data/download/all-states-history.csv")
DF <- DF %>% mutate(day = yday(ymd(date))+366*(year(date)==2021))

states <- c("OH","IN")
names <- c("Ohio", "Indiana")
pop <- c(11790587,6785644) # data from https://www2.census.gov/programs-surveys/popest/tables/2020-2021/state/totals/NST-EST2021-POP.xlsx

#entry parameters

start_date <- "2020-01-01"
end_date <- "2021-03-07"


ifr_min <- 0
ifr_max <- 0.03

p <- 0 # it was meant to be a cutoff for fraction of population tested

len_per_temp <- 6 # length of periods for cases/tests 7-1
cutoff <- 40 # longest time to death

#add_days <- max(c(0,cutoff - (first_death - first_day)))
add_days <- 0
cut <- 0  # for predicting deaths

#time to death distribution (truncated negative binomial)
alpha <- 21 
beta <- 1.1 
denom <- pnbinom(cutoff, size = alpha, prob=1/(beta+1))
ttd <- dnbinom(0:cutoff,size=alpha, prob = 1/(beta+1))/denom # normalization to 1 of truncated distribution


# main loop for extracting state data
for(s in 1:length(states)){
  State <- states[s]
  name <- names[s]
  N = pop[s] #population of state
  df <- DF %>% filter(state == State) %>% # selection of required variables
    select(state,date,day,positive,positiveIncrease,
           totalTestResults, totalTestResultsIncrease,
           death,deathIncrease)
  df[is.na(df)] <- 0
  
 
  ### data cleaning
  
  cum_deaths <- c(df$death,0) # death column contains already cumulated count
  for (i in 2:length(cum_deaths)) {
    cum_deaths[i] <- min(cum_deaths[1:i], na.rm = T)
  }
  df$deathIncrease <- rev(diff(rev(cum_deaths))) # lagged difference, so we treat death column as base for deathIncrease calculations (prob for consistency)
  
  for(i in nrow(df):1){ # fixing situations when positive/total results on specific day t are negative by scaling rest of n-t results
    if(df$positiveIncrease[i] < 0){
      cc_temp <- df$positiveIncrease[i]
      df$positiveIncrease[i] <- 0
      df$positiveIncrease[i:nrow(df)] = df$positiveIncrease[i:nrow(df)]*
        (1+cc_temp/sum(df$positiveIncrease[i:nrow(df)]))
    }
    if(df$totalTestResultsIncrease[i] < 0){
      test_temp <- df$totalTestResultsIncrease[i]
      df$totalTestResultsIncrease[i] <- 0
      df$totalTestResultsIncrease[i:nrow(df)] =
        df$totalTestResultsIncrease[i:nrow(df)]*
        (1+test_temp/sum(df$totalTestResultsIncrease[i:nrow(df)]))
    }
  }
  
  # moving average [note that we don't use that, but simply ]
  df$positiveIncrease <- ma(df$positiveIncrease,0) 
  df$totalTestResultsIncrease <- ma(df$totalTestResultsIncrease,0)
  # cumulative summing
  df$positive <- rev(cumsum(rev(df$positiveIncrease)))
  df$totalTestResults <- rev(cumsum(rev(df$totalTestResultsIncrease)))
  
  df <- df %>% filter(date >= start_date, date <= end_date)
  
  ### setup initial dates
  first_death <- min(df$day[df$deathIncrease > 0])
  first_day <- min(df$day)
  last_day <- max(df$day)
  new_first_day <- first_day - add_days
  
  if(add_days > 0){d_days <- c(new_first_day:(first_day-1),rev(df$day))}
  if(add_days==0){d_days <- rev(df$day)}
  
  ### for predicting deaths - probably dispensable
  new_death_day <- last_day - cut # yday(end_date)-cut
  d_days <- d_days[d_days <= new_death_day]
  n_d_days <- length(d_days)
  d <- c(rep(0,add_days),rev(df$deathIncrease[df$day <= new_death_day]))
  
  ### cases/tests
  
  # test_days <- rev(df$day)
  
  tests_pre <- rev(df$totalTestResultsIncrease[df$totalTestResults/N <= p])
  cc_pre <- rev(df$positiveIncrease[df$totalTestResults/N <= p])
  cc_pre_total <- sum(cc_pre)
  tests_pre_total <- sum(tests_pre)
  
  cc_post <- rev(df$positiveIncrease[df$totalTestResults/N > p])
  tests_post <- rev(df$totalTestResultsIncrease[df$totalTestResults/N > p])
  test_days <- rev(df$day[df$totalTestResults/N > p])
  
  cum_tests <- rev(df$totalTestResults)
  cum_cc <- rev(df$positive)
  tests <- rev(df$totalTestResultsIncrease)
  cc <- rev(df$positiveIncrease)
  
  n_test_days <- length(test_days)
  
  new_first_day <- min(c(min(d_days),min(test_days)))
  shifted_d_days <- d_days - new_first_day + 1
  shifted_test_days <- test_days - new_first_day + 1
  
  n_days <- max(c(max(shifted_d_days),max(shifted_test_days)))
  days <- new_first_day:last_day
  
  dates = as.Date(days,origin="2019-12-31")
  
  # segment case/test data into periods
  len_per <- len_per_temp  # length of periods in days 6?
  cum_tests_per <- c(0,0)
  
  #while(is.element(0, cum_tests_per[-1])){
  while(prod(cum_tests_per[-1] > 0)==0){
    len_per <- len_per+1
    cc_per <- split(cc_post, ceiling(seq_along(cc_post)/len_per)) # cases in each period
    tests_per <- split(tests_post, ceiling(seq_along(tests_post)/len_per)) # tests in each period
    cum_cc_per <- unlist(lapply(cc_per,sum),use.names = FALSE) # sum of cases in each period
    cum_tests_per <- c(tests_pre_total, unlist(lapply(tests_per,sum),use.names = FALSE)) # sum of tests in each period
    per_ends <- c(shifted_test_days[1]-1,
                    shifted_test_days[cumsum(unlist(lapply(cc_per,length)))]) # (shifted) end date of each period
    n_per <- length(per_ends)-1
  }
  # save data from the states
  ctp_data <- list(state = State, name=name, N=N,
                   cutoff = cutoff, cut = cut, 
                   cc = cc, cum_cc = cum_cc,
                   tests = tests, cum_tests = cum_tests,
                   d = d, days = days,
                   n_days = n_days, d_days = d_days,
                   n_d_days = n_d_days, shifted_d_days = shifted_d_days,
                   test_days = test_days, n_test_days=n_test_days,
                   shifted_test_days=shifted_test_days, ttd = ttd,
                   start_date = start_date,end_date = end_date, 
                   dates = dates, len_per = len_per, n_per = n_per, 
                   cc_per = cc_per,cum_cc_per = cum_cc_per, 
                   tests_per = tests_per,cum_tests_per = cum_tests_per, 
                   per_ends = per_ends, tests_pre = tests_pre,
                   cc_pre = cc_pre, cc_pre_total = cc_pre_total,
                   tests_pre_total = tests_pre_total, p=p,
                   cc_post = cc_post, tests_post = tests_post,
                   ifr_min=ifr_min, ifr_max=ifr_max)
  save(ctp_data, file = paste0("data/",State,"_ctp.RData"))
}